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


# Statistical modeling
library(lmtest)        # Statistical tests (e.g., coeftest)
library(sandwich)      # Robust / clustered SEs (e.g., vcovCL)
library(survival)      # Survival analysis (Cox models, Surv objects, coxph, cox.zph)
library(broom)         # Tidy model outputs (tidy(), glance(), augment() for regression/survival models)


# Table and summary formatting
library(gtsummary)     # Table formatting for publications
library(modelsummary)  # Compact model tables

# Visualization
library(ggplot2)       # Plotting

# =============================
# 1) Read and format the data
# =============================
Sheets <- read.delim(file   = "../fda-pediatric-studies/data/fda_drug_approvals_2011_2023.txt", 
                     header = TRUE, 
                     sep    = "\t")

label_data <- read.delim(file   = "../fda-pediatric-studies/data/pediatric_labeling_additions.txt", 
                        header = TRUE, 
                        sep    = "\t")


# =================================
# 2) Subset: PMRs due by 2024-12-31
# =================================

# Restrict to PMRs with a due date on or before 2024-12-31
december_due <- Sheets %>%
  dplyr::filter(
    !is.na(Study.report.submission.due.date) & 
      Study.report.submission.due.date <= as.Date("2024-12-31"))

dec_due_count <- nrow(december_due)
nrow(december_due)

# Follow-up time (years) from approval to the censoring date (2024-12-31)
december_due <- december_due %>%
  dplyr::mutate(
    follow_up = as.numeric(difftime(as.Date("2024-12-31"), `Approval.Date`, units = "days")) / 365.25
  )

# Harmonize study identifier across original and replacement PMRs
december_due <- december_due %>%
  dplyr::mutate(
    Study.ID = coalesce(
      Study.ID,
      New.PMR..2,
      New.PMR.replacing.superseded.PMR
    )
  )

# Merge PMR data with labeling-change data
merged_df <- merge(
  december_due,
  label_data,
  by = "Study.ID",
  all.x = TRUE
)


# ==============================
# 3) Labeling change variables
# ==============================
merged_df <- merged_df %>%
  dplyr::mutate(
    # Label change date: blank/"N" -> NA, parse m/d/y
    label_change_date =
      Resulting.label.change.from.PMR..N.if.not..date.if.Y. %>%
      str_trim() %>%                 # remove stray spaces
      na_if("N") %>%                 # turn "N" into NA
      mdy(quiet = TRUE),             # parse m/d/y or m/d/Y -> Date
    # Flag: did any label change occur?
    label_change_flag = !is.na(label_change_date), 
    # Time from drug approval to labeling addition (years)
    time_to_label_addition_by_study = as.numeric(difftime(label_change_date, `Approval.Date`, units = "days")) / 365.25)  %>%
  dplyr::mutate(
    Study.classification = factor(Study.classification),
    ATC_category         = factor(ATC_category),
    BLA_category         = factor(ifelse(BLA == "N", "Non-Biologic",
                                         ifelse(is.na(BLA), NA, "Biologic"))),
    approval_year        = year(Approval.Date),
    approval_year_c      = approval_year - mean(approval_year, na.rm = TRUE)
  ) 


# ==================================
# 4) Collapse to drug-level dataset
# ==================================
# Columns needed at drug level
merged_df_subset <- merged_df %>%
      dplyr::select( Drug.name, label_change_date, label_change_flag, time_to_label_addition_by_study,
                     Neonate, Infant, Early.childhood,Late.childhood, Adolescent, Unspecified, 
                     ATC_category, BLA_category, approval_year,approval_year_c, 
                     PREA.labeling.change.type.x, PREA.labeling.change.type.y) %>%
      unique()

# Helper functions for collapsing study-level info to drug-level

# earliest non-NA; if all NA -> NA (for Dates/numerics)
min_non_na <- function(x) if (all(is.na(x))) NA else suppressWarnings(min(x, na.rm = TRUE))

# first non-NA value; NA if all NA (for scalar text/nums)
first_non_na <- function(x) { i <- which(!is.na(x)); if (length(i)) x[min(i)] else NA }

# any "yes" in a vector (Y/Yes/1/TRUE)
any_yes <- function(x) {
  if (is.logical(x)) return(any(x, na.rm = TRUE))
  if (is.numeric(x)) return(any(x == 1, na.rm = TRUE))        # tighten to ==1; change to x != 0 if preferred
  xr <- toupper(trimws(as.character(x)))
  any(xr %in% c("Y","YES","T","TRUE","1"), na.rm = TRUE)
}

# collapse distinct values into comma-separated string; NA if none
collapse_or_na <- function(x) {
  u <- sort(unique(na.omit(x)))
  if (length(u) == 0) NA_character_ else paste(u, collapse = ", ")
}

age_cols <- c( "Neonate", "Infant", "Early.childhood", 
               "Late.childhood", "Adolescent", "Unspecified" )

# Collapse to one row per drug
dedup <- merged_df_subset %>%
  dplyr::group_by(Drug.name) %>%
  dplyr::summarise(
    # Age groups: Y if any row says yes
    across(all_of(age_cols), ~ if_else(any_yes(.), "Y", "N")),
    
    # First non-NA across type_x OR type_y; stays NA if none exist
    prea_labeling_change_first = {
      v <- dplyr::coalesce(PREA.labeling.change.type.x, PREA.labeling.change.type.y)
      nn <- which(!is.na(v))
      if (length(nn) > 0) v[[nn[1]]] else NA_character_
    },
    
    # Label change flag: Y if any row has Y
    label_change_flag = if_else(any_yes(label_change_flag), "Y", "N"),
    
    # Earliest label change date / time to label addition
    label_change_date = min_non_na(label_change_date),
    time_to_label_addition_by_study = min_non_na(time_to_label_addition_by_study),
    
    # ATC / BLA categories and original change-type columns
    ATC_category = collapse_or_na(ATC_category),
    BLA_category = collapse_or_na(BLA_category),
    PREA.labeling.change.type.x = collapse_or_na(PREA.labeling.change.type.x),
    PREA.labeling.change.type.y = collapse_or_na(PREA.labeling.change.type.y),
    
    # Approval year info
    approval_year   = first_non_na(approval_year),
    approval_year_c = first_non_na(approval_year_c),
    
    .groups = "drop"
  )

# Total distinct drugs (for optional % of all drugs)
n_total <- dplyr::n_distinct(dedup$Drug.name)

# ======================================================================
# PART A — Type of labeling addition
# ======================================================================
# Denominator per labeling-change type
type_denoms <- dedup %>%                              
  dplyr::filter(!is.na(prea_labeling_change_first)) %>%
  dplyr::group_by(prea_labeling_change_first) %>%
  dplyr::summarise(n_stratum = n_distinct(Drug.name), .groups = "drop")

# Among those with a labeling addition, summarize counts & time to addition
type_summary <- dedup %>%
  dplyr::filter(label_change_flag == "Y", !is.na(prea_labeling_change_first)) %>%
  dplyr::group_by(prea_labeling_change_first) %>%
  dplyr::summarise(
    n_with_addition = n_distinct(Drug.name),
    time_mean = mean(time_to_label_addition_by_study, na.rm = TRUE),
    time_sd   = sd(time_to_label_addition_by_study, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::left_join(type_denoms, by = "prea_labeling_change_first") %>%
  dplyr::mutate(
    pct_within_stratum = 100 * n_with_addition / n_stratum
  ) %>%
  dplyr::arrange(desc(n_with_addition))

type_summary

# One-way ANOVA: time to label addition by type of labeling change
anova_df <- dedup %>%
  dplyr::filter(label_change_flag == "Y",
         !is.na(prea_labeling_change_first),
         !is.na(time_to_label_addition_by_study)) %>%
  select(Drug.name, prea_labeling_change_first, time_to_label_addition_by_study)

anova_result <- aov(time_to_label_addition_by_study ~ prea_labeling_change_first, data = anova_df)
summary(anova_result)

# Post-hoc pairwise comparisons (Tukey)
tukey_result <- TukeyHSD(anova_result)

tukey_formatted <- broom::tidy(tukey_result) %>%
  dplyr::mutate(
    comparison   = contrast,
    mean_diff    = estimate,
    p_value_fmt  = ifelse(adj.p.value < 0.001, "<0.001",
                          sprintf("%.3f", adj.p.value)),
    ci           = sprintf("%.2f to %.2f", conf.low, conf.high),
    report       = sprintf(
      "%s: mean difference = %.2f years (95%% CI %.2f–%.2f), p = %s",
      comparison, mean_diff, conf.low, conf.high, p_value_fmt
    )
  ) %>%
  select(comparison, mean_diff, conf.low, conf.high, adj.p.value, report)

tukey_formatted

# ======================================================================
# PART B — Pediatric age groups
# ======================================================================
age_long <- dedup %>%
  select(Drug.name, label_change_flag, time_to_label_addition_by_study, all_of(age_cols)) %>%
  pivot_longer(cols = all_of(age_cols),
               names_to = "age_group", values_to = "in_group_flag")

# Denominator per age group (how many drugs in each group)
age_denoms <- age_long %>%
  dplyr::filter(in_group_flag == "Y") %>%
  dplyr::group_by(age_group) %>%
  dplyr::summarise(n_stratum = n_distinct(Drug.name), .groups = "drop")

# Among those, how many had a labeling addition + time stats
age_summary <- age_long %>%
  dplyr::filter(in_group_flag == "Y", label_change_flag == "Y") %>%
  dplyr::group_by(age_group) %>%
  dplyr::summarise(
    n_with_addition = n_distinct(Drug.name),
    time_mean = round(mean(time_to_label_addition_by_study, na.rm = TRUE),1),
    time_sd   = round(sd(time_to_label_addition_by_study, na.rm = TRUE),1),
    .groups = "drop"
  ) %>%
  dplyr::left_join(age_denoms, by = "age_group") %>%
  dplyr::mutate(
    pct_within_stratum = 100 * n_with_addition / n_stratum
  ) %>%
  dplyr::arrange(desc(n_with_addition))

age_summary

# ======================================================================
# PART C — Therapeutic areas (ATC)
# ======================================================================
# If a drug has multiple ATC categories in a single cell (comma-separated),
# split them into multiple rows so each category is counted correctly.
atc_denoms <- dedup %>%
  dplyr::filter(!is.na(ATC_category), ATC_category != "") %>%
  separate_rows(ATC_category, sep = "\\s*,\\s*") %>%
  dplyr::group_by(ATC_category) %>%
  dplyr::summarise(n_stratum = n_distinct(Drug.name), .groups = "drop")

atc_summary <- dedup %>%
  dplyr::filter(label_change_flag == "Y", !is.na(ATC_category), ATC_category != "") %>%
  separate_rows(ATC_category, sep = "\\s*,\\s*") %>%
  dplyr::group_by(ATC_category) %>%
  dplyr::summarise(
    n_with_addition = n_distinct(Drug.name),
    time_mean = round(mean(time_to_label_addition_by_study, na.rm = TRUE),1),
    time_sd   = round(sd(time_to_label_addition_by_study, na.rm = TRUE),1),
    .groups = "drop"
  ) %>%
  dplyr::left_join(atc_denoms, by = "ATC_category") %>%
  dplyr::mutate(
    pct_within_stratum = 100 * n_with_addition / n_stratum
  ) %>%
  dplyr::arrange(desc(n_with_addition))

atc_summary

# ======================================================================
# PART D — Regression models (logistic & Cox)
# ======================================================================

# Helper: map Y/N/TRUE/FALSE/1/0 to 0/1
to01 <- function(x) {
  y <- toupper(trimws(as.character(x)))
  as.integer(y %in% c("Y","YES","T","TRUE","1"))
}

# Data for regression
data_tbl3 <- dedup %>%
  dplyr::mutate(
    # outcome (1 = has labeling addition)
    y = to01(label_change_flag),
    
    # Covariates
    BLA_category = fct_relevel(as.factor(BLA_category), "Non-Biologic"),
    ATC_category = as.factor(ATC_category),
    ATC_category = fct_infreq(ATC_category) %>% fct_relevel(levels(.)[1]),
    
    # ensure approval_year_c defined if missing
    approval_year_c = if (!"approval_year_c" %in% names(.)) {
      approval_year - mean(approval_year, na.rm = TRUE)
    } else approval_year_c
  ) %>%
  # age flags to 0/1
  dplyr::mutate(across(all_of(age_cols), to01)) %>%
  # drop rows with missing essentials
  tidyr::drop_na(y, BLA_category, ATC_category, approval_year_c)

# ----------------------------- #
# Logistic model: presence of labeling addition
# ----------------------------- #
f_glm <- as.formula(
  paste("y ~ ATC_category + BLA_category + approval_year_c +", 
        paste(age_cols, collapse = " + "))
)

fit_tbl3_glm <- glm(f_glm, family = binomial, data = data_tbl3)

tbl3_or <- tidy(fit_tbl3_glm, conf.int = TRUE, exponentiate = TRUE) %>%
  dplyr::filter(term != "(Intercept)") %>%
  transmute(
    term,
    OR = estimate,
    CI_low = conf.low,
    CI_high = conf.high,
    p.value,
    `Adjusted OR (95% CI)` = sprintf("%.2f (%.2f–%.2f)", OR, CI_low, CI_high)
  )

tbl3_or

# ----------------------------- #
# Cox model: time to labeling addition
# ----------------------------- #
cox_df <- data_tbl3 %>%
  dplyr::mutate(
    event = y,
    followup_years = max(time_to_label_addition_by_study, na.rm = TRUE),
    time = ifelse(event == 1L, time_to_label_addition_by_study, followup_years)
  ) %>%
  tidyr::drop_na(time)

f_cox <- as.formula(
  paste("Surv(time, event) ~ ATC_category + BLA_category + approval_year_c +", 
        paste(age_cols, collapse = " + "))
)


fit_tbl3_cox <- coxph(f_cox, data = cox_df)

# Test proportional hazards assumption
ph_test <- cox.zph(fit_tbl3_cox)
print(ph_test)

tbl3_hr <- tidy(fit_tbl3_cox, exponentiate = TRUE, conf.int = TRUE) %>%
  dplyr::filter(term != "(Intercept)") %>%
  transmute(
    term,
    HR = estimate,
    CI_low = conf.low,
    CI_high = conf.high,
    p.value,
    `Adjusted HR (95% CI)` = sprintf("%.2f (%.2f–%.2f)", HR, CI_low, CI_high)
  )

tbl3_hr

# ====================================================
# Cox model with time-varying effect of approval year
# (to address PH violation for approval_year_c)
# ====================================================
fit_tbl3_cox_adj <- coxph(
  Surv(time, event) ~ ATC_category + BLA_category + 
    tt(approval_year_c) + # time-varying effect
    Neonate + Infant + Early.childhood + Late.childhood +
    Adolescent + Unspecified,
  data = cox_df,
  tt = function(x, t, ...) x * log(t)
)

tbl3_hr_adj <- tidy(fit_tbl3_cox_adj, exponentiate = TRUE, conf.int = TRUE) %>%
  dplyr::filter(term != "(Intercept)") %>%
  transmute(
    term,
    HR = estimate,
    CI_low = conf.low,
    CI_high = conf.high,
    p.value,
    `Adjusted HR (95% CI)` = sprintf("%.2f (%.2f–%.2f)", HR, CI_low, CI_high)
  )

tbl3_hr_adj

sessionInfo()
