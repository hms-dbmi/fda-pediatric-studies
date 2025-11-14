# ======================================================================
# SETUP: Packages loading
# ======================================================================
# Core tidyverse and utility packages
library(dplyr)         # Data manipulation
library(tidyr)         # Data tidying (pivot_longer, etc.)
library(purrr)         # Functional programming (map, etc.)
library(stringr)       # String operations
library(janitor)       # Cleaning column names
library(lubridate)     # Date handling
library(forcats)       # Factor ordering
library(tibble)

# Statistical modeling
library(lmtest)        # Statistical tests (e.g., coeftest)
library(sandwich)      # Robust / clustered SEs

# Table and summary formatting
library(gtsummary)     # Table formatting for publications
library(modelsummary)  # Compact model tables

# Visualization
library(ggplot2)       # Plotting

# ===================================================
# 1) Read, format the data and set up reference_date
# ===================================================
Sheets <- read.delim(file   = "../completed-fda-pediatric-studies/data/fda_drug_approvals_2011_2023.txt", 
                     header = TRUE, 
                     sep    = "\t")

# Format dates and numeric data as by default all is character
Sheets <- Sheets %>%
  dplyr::mutate(Approval.Date = as.Date(Approval.Date), 
                Study.report.submission.due.date = as.Date(Study.report.submission.due.date),
                Submission.date = as.Date(Submission.date), 
                original_duration_years = as.numeric( original_duration_years))

# Administrative censoring date used in survival analyses
reference_date <- as.Date("2024-12-31")

##############
## Figure 1 ##
##############


# -------------------------------------------------------------------
# Overall KM (not stratified): time from approval to study submission
# or censoring at reference_date
# -------------------------------------------------------------------
data <- Sheets %>%
  dplyr::select(Drug.name, Approval.Date, Study.classification, Submission.date) %>%
  # Event indicator: 1 if submission date present, 0 otherwise
  dplyr::mutate(Submitted = ifelse(!is.na(Submission.date), 1, 0)) %>%
  # Update the data to handle cases with NA 'Submission date'
  dplyr::mutate(
    # If 'Submission date' is not NA, calculate time to submission date
    # If 'Submission date' is NA, calculate time to December 31, 2024
    time_to_event = ifelse(
      is.na(Submission.date), 
      as.numeric(difftime(reference_date, Approval.Date, units = "days")),
      as.numeric(difftime(Submission.date, Approval.Date, units = "days"))
    ),
    event = Submitted  # Event indicator (1 for submitted, 0 for not submitted)
  ) %>%
  # Convert time to years for interpretability
  dplyr::mutate( time_to_event_years = time_to_event / 365.25 )

# Survival object for KM
surv_object <- Surv(time = data$time_to_event_years, event = data$event)

# KM fit without stratification (overall curve)
km_fit <- survfit(surv_object ~ 1)  # '~ 1' means no stratification, overall survival

# -------------------------------------------------------------------
# Build two datasets for expected vs observed completion:
#   - data_real: observed submission times
#   - data_due: expected completion dates (PMR due dates)
# Then combine and fit stratified KM curves by Source
# -------------------------------------------------------------------
data_real <- Sheets %>%
  dplyr::select(Drug.name, Approval.Date, Study.ID, Study.classification, 
                Submission.date) %>%
  dplyr::mutate(Submitted = ifelse(!is.na(Submission.date), 1, 0)) %>%
  dplyr::mutate(
    # Time to event or censoring as above
    time_to_event = ifelse(
      is.na(Submission.date), 
      as.numeric(difftime(reference_date, Approval.Date, units = "days")),
      as.numeric(difftime(Submission.date, Approval.Date, units = "days"))
    ),
    event = Submitted  
  ) %>%
  # For plotting, use a generic name for the observed date column
  dplyr::rename( Due.date.or.submission.date  = Submission.date) %>%
  # Label this source as "data_real" (observed)
  dplyr::mutate(Source = factor("data_real", levels = c("data_due", "data_real")))


data_due <- Sheets %>%
  dplyr::select(Drug.name, Approval.Date, Study.ID, Study.classification, 
                Submission.date, Study.report.submission.due.date) %>%
  # Submitted flag not needed here; we treat due date as an "event" for expected timing
  dplyr::mutate(Submitted = NA) %>%
  # Keep only rows with a defined study report submission due date
  dplyr::filter(!is.na(Study.report.submission.due.date)) %>%
  dplyr::mutate(
    # Expected time: from approval to due date
    time_to_event = as.numeric(difftime(Study.report.submission.due.date, Approval.Date, units = "days")),
    event = 1  # Event indicator (1 since the 'Study report submission due date' is available)
  ) %>%
  # Rename for consistency with data_real
  dplyr::rename( Due.date.or.submission.date = Study.report.submission.due.date ) %>%
  # Label this source as "data_due" (expected)
  dplyr::mutate(Source = factor("data_due", levels = c("data_due", "data_real")))

# Combine expected + observed study-level datasets
combined_data <- bind_rows(data_real, data_due)

# Convert time to years
combined_data$time_to_event_years <- combined_data$time_to_event / 365.25

# Survival object for combined data
surv_object_combined <- Surv(time = combined_data$time_to_event_years, event = combined_data$event)

# KM fit stratified by source (expected vs observed)
km_fit_combined <- survfit(surv_object_combined ~ Source, data = combined_data)

# -------------------------------------------------------------------
# Plot KM curves and compute cumulative incidence at 3, 5, 10 years
# -------------------------------------------------------------------

# Global p-value comparing KM curves (log-rank)
pval <- surv_pvalue(km_fit_combined)$pval
formatted_pval <- ifelse(pval < 0.001, "p < 0.001", sprintf("p = %.3f", pval))

# KM plot: cumulative incidence of study completion
ggsurvplot(km_fit_combined, 
           data = combined_data, 
           xlab = "Time from FDA Approval, y", 
           ylab = "Cumulative Incidence of Study Completion", 
           title = "Time to Study Submission of Pediatric Studies",
           risk.table = TRUE, 
           legend.labs = c("Expected completion date", "Observed completion date"),
           fun = "event",
           palette = c("#d35400", "#1f77b4"),
           pval = formatted_pval)  # Adding p-value for comparing the curves

# Cumulative incidence at 3 years
km_summary_3 <- summary(km_fit_combined, times = 3)
km_summary_at_3_years_by_source <- data.frame(
  Source = km_summary_3$strata,
  Cumulative_Incidence_at_3_Years = 1 - km_summary_3$surv[km_summary_3$time == 3]
)
print(km_summary_at_3_years_by_source)

# Cumulative incidence at 5 years
km_summary_5 <- summary(km_fit_combined, times = 5)
km_summary_at_5_years_by_source <- data.frame(
  Source = km_summary_5$strata,
  Cumulative_Incidence_at_5_Years = 1 - km_summary_5$surv[km_summary_5$time == 5]
)
print(km_summary_at_5_years_by_source)

# Cumulative incidence at 10 years
km_summary_10 <- summary(km_fit_combined, times = 10)
km_summary_at_10_years_by_source <- data.frame(
  Source = km_summary_10$strata,
  Cumulative_Incidence_at_10_Years = 1 - km_summary_10$surv[km_summary_10$time == 10]
)
print(km_summary_at_10_years_by_source)

##############################
## Paired analysis function ##
##############################

# Helper function for formatting p-values
fmt_p <- function(p) ifelse(p < 0.001, "<0.001", 
                            formatC(p, format = "f", digits = 3))

# -------------------------------------------------------------------
# paired_results():
#   - Takes a dataset with expected vs observed times per unit
#   - Computes:
#       * McNemar test at specified timepoints (3, 5, 10 years by default)
#       * Wilcoxon signed-rank test for median delay (observed - expected)
# -------------------------------------------------------------------
paired_results <- function(dat, cutoffs = c(3,5,10)) {
  
  dat <- dat %>%
    mutate(
      # Total follow-up available per study/drug (for eligibility at each cutoff)
      fup_y = as.numeric(difftime(reference_date, Approval.Date, units = "days"))/365.25,
      # Difference between observed and expected times (in years)
      diff_y = observed_time - expected_time
    )
  
  # ----- McNemar at each cutoff (restrict to sufficient follow-up) -----
  one_cutoff <- function(K) {
    # Only include rows with at least K years of follow-up
    elig <- dat %>% filter(fup_y >= K)
    if (nrow(elig) == 0) return(tibble())
    
    # Expected completion by time K (1 if expected_time <= K)
    expK <- as.integer(!is.na(elig$expected_time) & elig$expected_time <= K)
    
    # Observed completion by time K (1 if event occurred by K)
    obsK <- as.integer(!is.na(elig$observed_time) & elig$observed_event == 1 & elig$observed_time <= K)
    
    # 2x2 table of expected vs observed completion
    tab <- table(expK, obsK)
    mc  <- mcnemar.test(tab, correct = TRUE)
    
    tibble(
      Timepoint   = paste0(K, " years"),
      Expected    = round(mean(expK) * 100, 1),
      Observed    = round(mean(obsK) * 100, 1),
      P_value     = fmt_p(mc$p.value),
      Expected_n  = sum(expK), Observed_n = sum(obsK),
      N_eligible  = nrow(elig)
    )
  }
  
  # Apply McNemar at each cutoff
  mc_tbl <- bind_rows(lapply(cutoffs, one_cutoff))
  
  # ----- Wilcoxon signed-rank on paired differences (uses rows with both times) -----
  both <- dat %>% filter(!is.na(expected_time) & !is.na(observed_time))
  if (nrow(both) > 0) {
    wilx <- wilcox.test(both$observed_time, both$expected_time,
                        paired = TRUE, exact = FALSE, conf.int = TRUE)
    med  <- median(both$diff_y, na.rm = TRUE)
    ci   <- wilx$conf.int
    
    wilx_row <- tibble(
      Timepoint  = paste0("Median delay (observed − expected) [HL 95% CI: ",
                          round(ci[1], 2), ", ", round(ci[2], 2), "]"),
      Expected   = NA_real_,
      Observed   = round(med, 2),   # in years
      P_value    = fmt_p(wilx$p.value),
      Expected_n = NA_integer_, Observed_n = NA_integer_, N_eligible = nrow(both)
    )
    
    mc_tbl <- bind_rows(mc_tbl, wilx_row)
  }
  
  mc_tbl %>% select(Timepoint, Expected, Observed, P_value, Expected_n, Observed_n, N_eligible)
}

# Build paired dataset for expected vs observed study submission times
paired_f1 <- combined_data %>%
  transmute(
    Approval.Date,
    reference_date,
    expected_time = as.numeric(difftime(Due.date.or.submission.date, Approval.Date, units = "days"))/365.25,
    observed_time = as.numeric(difftime(Submission.date,Approval.Date, units = "days"))/365.25,
    observed_event = as.integer(!is.na(Submission.date))
  )

# Paired comparison table for Figure 1 at 3, 5, and 10 years
res_tbl_f1 <- paired_results(paired_f1, cutoffs = c(3,5,10))
res_tbl_f1

##############
## Figure 2 ##
##############
# Fix specific drug naming inconsistency for one PMR ID
Sheets$Drug.name[Sheets$Study.ID == "PMR 3292-2"] <- "Akynzeo2018"

# -------------------------------------------------------------------
# Observed labeling addition (data_lab):
#   - One row per drug
#   - Uses "All (non-released) PREA labeling additions complete (date)"
# -------------------------------------------------------------------
data_lab <- Sheets %>%
  dplyr::select(Drug.name, Approval.Date, Indication, Study.ID, 
                 Submission.date, All..non.released..PREA.labeling.additions.complete..date.) %>% 
  dplyr::distinct(Drug.name, .keep_all = TRUE) %>%
  dplyr::mutate( Final.Label.Added = ifelse(!is.na(All..non.released..PREA.labeling.additions.complete..date.), 1, 0)) %>%
  dplyr::mutate(
    time_to_addition = ifelse(
      is.na(All..non.released..PREA.labeling.additions.complete..date.), 
      as.numeric(difftime(reference_date, Approval.Date, units = "days")),
      as.numeric(difftime(All..non.released..PREA.labeling.additions.complete..date., Approval.Date, units = "days"))
    ),
    event = Final.Label.Added  # Event indicator (1 for final label added, 0 for not added)
  ) %>%
  dplyr::mutate(time_to_addition_years = time_to_addition / 365.25, 
                Source = factor("data_lab", levels = c("data_last_ind", "data_lab"))) %>%
  dplyr::select( Drug.name, Study.ID, Approval.Date, Indication, time_to_addition, event, time_to_addition_years, Source)

# Survival object for observed label additions
surv_object_lab <- Surv(time  = data_lab$time_to_addition_years, 
                        event = data_lab$event)

km_fit_lab <- survfit(surv_object_lab ~ 1)  # '~ 1' = no stratification, overall survival

# -------------------------------------------------------------------
# Expected labeling completion timing based on last due date per drug
# (data_last_ind)
# -------------------------------------------------------------------
data_last_ind <- Sheets %>%
  dplyr::select(Drug.name, Approval.Date, Indication, Study.ID, 
                Study.report.submission.due.date,
                All..non.released..PREA.labeling.additions.complete..date.) %>%
  dplyr::group_by( Drug.name ) %>%
  # Add a new column with the latest date within each group
  dplyr::mutate(Last_due_date_per_indication = max(Study.report.submission.due.date, na.rm = TRUE)) %>%
  # Ungroup to remove the grouping structure
  dplyr::ungroup() %>%
  dplyr::distinct( Drug.name, .keep_all = TRUE) %>%
  dplyr::group_by( Drug.name ) %>%
  dplyr::slice_head(n = 1) %>%
  dplyr::ungroup() %>%
  dplyr::mutate( Final.Label.Added = ifelse(!is.na(Last_due_date_per_indication), 1, 0),
                 time_to_addition  = as.numeric(as.Date(Last_due_date_per_indication) - as.Date(Approval.Date)),
    event = Final.Label.Added
  ) %>%
  dplyr::mutate( time_to_addition_years = time_to_addition / 365.25, 
                 Source = factor("data_last_ind", levels = c("data_last_ind", "data_lab"))) %>%
  dplyr::select( Drug.name, Study.ID, Approval.Date, Indication, time_to_addition, event, time_to_addition_years, Source)

# Combine expected (last_due_date) and observed (data_lab)
combined_data_labels <- bind_rows(data_last_ind, data_lab)

# Survival object for combined label data
surv_object_combined_labels <- Surv(time  = combined_data_labels$time_to_addition_years, 
                                    event = combined_data_labels$event)

# KM fit stratified by Source (expected vs observed labeling additions)
km_fit_combined_labels <- survfit(surv_object_combined_labels ~ Source, data = combined_data_labels)

# Plot and cumulative incidence numbers for Figure 2

# p-value for difference between expected vs observed label addition curves
pval_labels <- surv_pvalue(km_fit_combined_labels)$pval_labels
formatted_pval_labels <- ifelse(pval < 0.001, "p < 0.001", sprintf("p = %.3f", pval))


ggsurvplot(km_fit_combined_labels, 
           data = combined_data_labels, 
           xlab = "Time from FDA Approval, y", 
           ylab = "Cumulative Incidence of Label Additions", 
           title = "Time to Label Addition: Actual vs. Due Date",
           risk.table = TRUE, 
           legend.labs = c("Expected completion date", "Observed labeling addition date"),
           fun = "event",
           palette = c("#d35400", "#1f77b4"),
           pval = formatted_pval_labels)  # Adding p-value for comparing the curves

# -------------------------------------------------------------------
# Paired analysis for Figure 2 (expected vs observed labeling addition)
# -------------------------------------------------------------------
paired_f2 <- inner_join(
  data_last_ind %>%
    transmute( Drug.name,
               Approval.Date = as.Date(Approval.Date),
              expected_time = time_to_addition_years),
  data_lab %>%
    transmute( Drug.name,
              observed_time = time_to_addition_years,
              observed_event = as.integer(event)),
  by = "Drug.name"
) %>%
  mutate(reference_date = reference_date)

# Paired comparison table for Figure 2 at 3, 5, and 10 years
res_tbl_f2 <- paired_results(paired_f2, cutoffs = c(3,5,10))
res_tbl_f2

sessionInfo()

