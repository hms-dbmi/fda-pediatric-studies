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
library(sandwich)      # Robust / clustered SEs

# Table and summary formatting
library(gtsummary)     # Table formatting for publications
library(modelsummary)  # Compact model tables

# Visualization
library(ggplot2)       # Plotting

# =============================
# 1) Read and format the data
# =============================
Sheets <- read.delim(file   = "../fda-pediatric-studies/data/t1.txt", 
                     header = TRUE, 
                     sep    = "\t")

# Format dates and numeric data as by default all is character
Sheets <- Sheets %>%
  dplyr::mutate(Approval.Date = as.Date(Approval.Date), 
                Study.report.submission.due.date = as.Date(Study.report.submission.due.date),
                Submission.date = as.Date(Submission.date), 
                original_duration_years = as.numeric( original_duration_years))

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

summary(december_due$follow_up)
sum(december_due$Submission.date <= december_due$Study.report.submission.due.date, na.rm = TRUE)

# Flags for "completed" and "completed by due date"
december_due <- december_due %>%
  dplyr::mutate(
    completed = ifelse(Status.of.current.PMR %in% c("Fulfilled", "Submitted"), "Y", "N"),
    completed_by_due_date = ifelse(
      !is.na(Submission.date) &
        !is.na(Study.report.submission.due.date) &
        Submission.date <= Study.report.submission.due.date,
      "Y", "N"
    )
  )

# Subsets
on_time   <- december_due %>%
  filter(completed_by_due_date == "Y")

completed <- december_due %>%
  filter(completed == "Y")

# ==========================
# 3) Descriptive counts / % 
# ==========================

## 3.1 By Study classification
by_study_type <- bind_rows(
  on_time      %>% transmute(class = Study.classification, set = "on_time"),
  completed    %>% transmute(class = Study.classification, set = "completed"),
  december_due %>% transmute(class = Study.classification, set = "december_due")
) %>%
  dplyr::filter(!is.na(class)) %>%                                # drop NAs like your == checks did
  dplyr::count(set, class) %>%                                    # counts within each set
  tidyr::pivot_wider(names_from = set, values_from = n, values_fill = 0) %>%
  dplyr::mutate(
    pct_on_time_of_due   = 100 * if_else(december_due > 0, on_time   / december_due, NA_real_),
    pct_completed_of_due = 100 * if_else(december_due > 0, completed / december_due, NA_real_)
  ) %>%
  dplyr::arrange(desc(december_due)) %>%
  dplyr::rename(Study.classification = class)

by_study_type

## 3.2 By age groups
age_cols <- c( "Neonate", "Infant", "Early.childhood", 
               "Late.childhood", "Adolescent", "Unspecified" )

# Helper: summarise one age-flag column across the three sets
summarize_one_age_col <- function(col) {
  on  <- sum(on_time[[col]]      == "Y", na.rm = TRUE)
  cm  <- sum(completed[[col]]    == "Y", na.rm = TRUE)
  due <- sum(december_due[[col]] == "Y", na.rm = TRUE)
  
  tibble::tibble(
    age_group_col = col,
    on_time       = on,
    completed     = cm,
    december_due  = due,
    pct_on_time_of_due   = ifelse(due > 0, 100 * on / due, NA_real_),
    pct_completed_of_due = ifelse(due > 0, 100 * cm / due, NA_real_)
  )
}

# By age summary
by_age_flags_all <- purrr::map_dfr(age_cols, summarize_one_age_col)
by_age_flags_all


## 3.3 By therapeutic area (ATC_category)
by_therapeutic_area <- bind_rows(
  on_time      %>% transmute(class = ATC_category, set = "on_time"),
  completed    %>% transmute(class = ATC_category, set = "completed"),
  december_due %>% transmute(class = ATC_category, set = "december_due")
) %>%
  dplyr::filter(!is.na(class)) %>%                                # drop NAs like your == checks did
  dplyr::count(set, class) %>%                                    # counts within each set
  tidyr::pivot_wider(names_from = set, values_from = n, values_fill = 0) %>%
  dplyr::mutate(
    pct_on_time_of_due   = 100 * if_else(december_due > 0, on_time   / december_due, NA_real_),
    pct_completed_of_due = 100 * if_else(december_due > 0, completed / december_due, NA_real_)
  ) %>%
  dplyr::arrange(desc(december_due)) %>%
  dplyr::rename(atc_category = class)

by_therapeutic_area

# =====================================================
# 4) Regression analysis (logistic with clustered SEs)
# =====================================================

# Base data frame for regression
data_for_glm_base <- december_due %>%
  mutate(
    # Make sure Approval.Date is Date (some pipelines add " UTC")
    Approval.Date = ymd(sub("\\s*UTC$", "", as.character(Approval.Date))),
    
    ATC_category         = factor(ATC_category),
    Study.classification = factor(Study.classification),
    Neonate              = factor(Neonate),
    Infant               = factor(Infant),
    Early.childhood      = factor(Early.childhood),
    Late.childhood       = factor(Late.childhood),
    Adolescent           = factor(Adolescent),
    Unspecified          = factor(ifelse(Unspecified == "N", "N",
                                         ifelse(is.na(Unspecified), NA, "Y"))),
    BLA_category         = factor(ifelse(BLA == "N", "Non-Biologic",
                                         ifelse(is.na(BLA), NA, "Biologic"))),
    approval_year        = year(Approval.Date),
    approval_year_c      = approval_year - mean(approval_year, na.rm = TRUE)
  )

# ===============================================
# 4.1 Helper: logistic model with clustered SEs
# ===============================================
#' Fit clustered logistic model and return ORs with 95% CI
#' @param df data.frame with predictors and Drug.name
#' @param outcome_col string; name of binary outcome column in df
#' @param drug_col string; name of clustering variable (drug identifier)
#' @param age_flag_cols character vector; names of age-flag columns

fit_glm_for <- function(df, 
                        outcome_col,
                        drug_col = "Drug.name",
                        age_flag_cols = age_cols) {
  
  ysym <- rlang::sym(outcome_col)
  
  # 1) build analysis frame: binary y + factors + 0/1 age flags
  df2 <- df %>%
    dplyr::transmute(
      y = dplyr::case_when(
        !!ysym %in% c(TRUE, 1, "1", "Y", "Yes", "yes") ~ 1L,
        !!ysym %in% c(FALSE, 0, "0", "N", "No", "no")  ~ 0L,
        TRUE ~ NA_integer_
      ),
      Study.classification, ATC_category, BLA_category, approval_year_c,
      dplyr::across(
        dplyr::all_of(age_flag_cols),
        ~ dplyr::case_when(
          . %in% c(TRUE, 1, "1", "Y", "Yes", "yes") ~ 1L,
          . %in% c(FALSE, 0, "0", "N", "No", "no", "I") ~ 0L,  # "I" treated as 0
          TRUE ~ 0L  # default to 0; change to NA_integer_ if you prefer to drop
        )
      )
    ) %>%
    tidyr::drop_na(y, Study.classification, ATC_category, BLA_category, approval_year_c) %>%
    dplyr::mutate(
      Study.classification = fct_infreq(Study.classification) %>% fct_relevel(levels(.)[1]),
      ATC_category         = fct_infreq(ATC_category)         %>% fct_relevel(levels(.)[1]),
      BLA_category         = fct_relevel(BLA_category, "Non-Biologic")
    )
  
  # 2) build formula: add all age flags as additive 0/1 covariates
  age_terms <- paste(age_flag_cols, collapse = " + ")
  fml <- as.formula(paste(
    "y ~ Study.classification + ATC_category + BLA_category + approval_year_c +", age_terms
  ))
  print(fml)
  
  # 3) fit
  m <- glm(fml, data = df2, family = binomial)
  
  # 4) cluster-robust vcov by drug
  V_cl <- sandwich::vcovCL(m, cluster = df2[[drug_col]])
  
  coefs <- stats::coef(m)
  se    <- sqrt(diag(V_cl))
  z     <- coefs / se
  pval  <- 2 * stats::pnorm(abs(z), lower.tail = FALSE)
  
  # 95% CI on log-odds, then exponentiate
  crit  <- stats::qnorm(0.975)
  lo    <- coefs - crit * se
  hi    <- coefs + crit * se
  
  tib_or <- tibble::tibble(
    term      = names(coefs),
    estimate  = exp(coefs),
    conf.low  = exp(lo),
    conf.high = exp(hi),
    p.value   = pval
  )
  
  list(
    model        = m,
    vcov_cluster = V_cl,
    or           = tib_or
  )
}

# ====================================================
# 4.2 Fit models: completed vs completed by due date
# ====================================================
res_completed           <- fit_glm_for(data_for_glm_base, "completed")
res_completed_by_duedate <- fit_glm_for(data_for_glm_base, "completed_by_due_date")

# ==================================
# 5) Format OR tables for reporting
# ==================================
fmt_num <- function(x, digits = 2) formatC(x, format = "f", digits = digits)
fmt_p <- function(p) ifelse(p < 0.001, "<0.001", formatC(round(p, 3), format = "f", digits = 3))

# Completed (any time)
or_tbl_completed <- res_completed$or %>%
  filter(term != "(Intercept)") %>%
  mutate(
    OR_CI = paste0(fmt_num(estimate, 2),
                   " (", fmt_num(conf.low, 2), "–", fmt_num(conf.high, 2), ")"),
    p_formatted = fmt_p(p.value)
  ) %>%
  select(term, OR_CI, p_formatted)

or_tbl_completed

# Completed by due date
or_tbl_completed_duedate <-  res_completed_by_duedate$or %>%
  filter(term != "(Intercept)") %>%
  mutate(
    OR_CI = paste0(fmt_num(estimate, 2),
                   " (", fmt_num(conf.low, 2), "–", fmt_num(conf.high, 2), ")"),
    p_formatted = fmt_p(p.value)
  ) %>%
  select(term, OR_CI, p_formatted)

or_tbl_completed_duedate


sessionInfo()