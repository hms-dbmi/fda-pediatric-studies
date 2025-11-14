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
Sheets <- read.delim(file   = "../fda-pediatric-studies/data/fda_drug_approvals_2011_2023.txt", 
                     header = TRUE, 
                     sep    = "\t")

# Format dates and numeric data as by default all is character
Sheets <- Sheets %>%
  dplyr::mutate(Approval.Date = as.Date(Approval.Date), 
                Study.report.submission.due.date = as.Date(Study.report.submission.due.date),
                Submission.date = as.Date(Submission.date), 
                original_duration_years = as.numeric( original_duration_years))

# =============================
# 2) Descriptive summaries
# =============================
total_PMRs <- sum(!is.na(Sheets$Drug.name)) 
total_PMRs

sum(!is.na(unique(Sheets$Drug.name))) 
sum(!is.na(Sheets$Indication))

# Study classification: N and %
Sheets %>%
  dplyr::group_by(Study.classification) %>%
  dplyr::summarize( n = n(), 
                    perc = round( n/total_PMRs*100, 2))
# Age groups
Sheets %>%
  tidyr::pivot_longer(
    cols = c(9:14), 
    names_to = "age_group", 
    values_to = "flag"
  ) %>%
  dplyr::filter(flag == "Y") %>%  
  dplyr::group_by(age_group) %>%
  dplyr::summarise(n = n(), .groups = "drop") %>%
  dplyr::mutate( perc = paste0(round(100 * n /total_PMRs, 1), "%"))

# ATC categories
Sheets %>%
  dplyr::group_by(ATC_category) %>%
  dplyr::summarize( n = n(), 
                    perc = round( n/total_PMRs*100, 2))

# Overall mean/median/SD for duration
Sheets %>%
  dplyr::summarise( mean(original_duration_years, na.rm = TRUE),
                    median(original_duration_years, na.rm = TRUE),
                    sd_duration = sd(original_duration_years, na.rm = TRUE))

# By BLA (Non-Biologic vs Biologic)
Sheets %>%
  filter(BLA == "N") %>%
  summarise(mean(original_duration_years, na.rm = TRUE),
            median(original_duration_years, na.rm = TRUE),
            sd_duration = sd(original_duration_years, na.rm = TRUE))

Sheets %>%
  filter(!BLA == "N") %>%
  summarise(mean(original_duration_years, na.rm = TRUE),
            median(original_duration_years, na.rm = TRUE),
            sd_duration = sd(original_duration_years, na.rm = TRUE))

Sheets$BLA_category <- ifelse(Sheets$BLA == "N", "Non-Biologic", "Biologic")

# Mean (95% CI) by Study classification
Sheets %>%
  dplyr::mutate( `Study.classification` = as.factor(`Study.classification`)) %>%
  dplyr::group_by( `Study.classification` ) %>%
  dplyr::summarize(
    n = sum(!is.na(original_duration_years)),
    sd_duration = sd(original_duration_years, na.rm = TRUE),
    mean_val = mean(original_duration_years, na.rm = TRUE),
    se = sd(original_duration_years, na.rm = TRUE) / sqrt(n),
    t_crit = qt(0.975, df = n - 1),  # 95% CI critical value
    ci_low = mean_val - t_crit * se,
    ci_high = mean_val + t_crit * se,
    mean_ci = sprintf("%.1f [%.1f, %.1f]", mean_val, ci_low, ci_high)
  ) %>%
  dplyr::select(`Study.classification`, n, mean_ci)

# Mean (95% CI) by ATC category
Sheets %>%
  dplyr::mutate( `ATC_category` = as.factor(`ATC_category`)) %>%
  dplyr::group_by( `ATC_category` ) %>%
  dplyr::summarize(
    n = sum(!is.na(original_duration_years)),
    mean_val = mean(original_duration_years, na.rm = TRUE),
    se = sd(original_duration_years, na.rm = TRUE) / sqrt(n),
    t_crit = qt(0.975, df = n - 1),  # 95% CI critical value
    ci_low = mean_val - t_crit * se,
    ci_high = mean_val + t_crit * se,
    mean_ci = sprintf("%.1f [%.1f, %.1f]", mean_val, ci_low, ci_high)
  ) %>%
  dplyr::select(ATC_category, n, mean_ci)

# =============================
# 3) Linear Model
# =============================
# Age related columns
age_cols <- c( "Neonate", "Infant", "Early.childhood", 
               "Late.childhood", "Adolescent", "Unspecified" )

# Build modeling frame
df_lin <- Sheets %>%
  # Convert the age columns to numeric indicators 0/1
  dplyr::mutate(across(all_of(age_cols),
                ~ case_when(. %in% c("Y","Yes",1,TRUE) ~ 1L,
                            . %in% c("N","No",0,FALSE) ~ 0L,
                            TRUE ~ NA_integer_))) %>%
  dplyr::transmute(
    y  = as.numeric(original_duration_years),
    Study.classification = factor(Study.classification),  # efficacy / safety / PK
    ATC_category         = factor(ATC_category),
    approval_year        = year(`Approval.Date`),
    BLA_category = factor(BLA_category),
    !!!syms(age_cols)   # the six 0/1 age indicators
  ) %>%
  tidyr::drop_na(y, Study.classification, ATC_category, BLA_category, approval_year, !!!syms(age_cols)) %>%
  mutate(
    # Reference levels: most frequent for study/ATC, Non-Biologic for BLA
    Study.classification = fct_relevel(fct_infreq(Study.classification), levels(fct_infreq(Study.classification))[1]),
    ATC_category         = fct_relevel(fct_infreq(ATC_category),         levels(fct_infreq(ATC_category))[1]),
    BLA_category         = fct_relevel(BLA_category, "Non-Biologic"),
    approval_year_c      = approval_year - mean(approval_year)
  )

# Model formula
form_lin <- reformulate(
  termlabels = c("Study.classification", age_cols, "ATC_category", "BLA_category", "approval_year_c"),
  response   = "y"
)
form_lin

# Fit the multivariable linear model 
m_lin <- lm(form_lin, data = df_lin)
summary(m_lin)

# Coefficient-level effects w/ 95% CI (standard SEs)
coef_lin <- broom::tidy(m_lin, conf.int = TRUE)
coef_lin

# Global p-values (Type II tests)
p_global_lin <- car::Anova(m_lin, type = 2) %>% broom::tidy() %>%
  dplyr::mutate(
    pval_sci = formatC(p.value, format = "e", digits = 2),
    pval_sig = ifelse(p.value < 0.05, "yes", "no"),
    pval_formatted = ifelse(p.value < 0.001, "<0.001",
                            formatC(round(p.value, 3), format = "f", digits = 3))
  )

p_global_lin

# ===================================
# 4) Cluster-robust SEs by Drug name
# ===================================
cov_drug <- vcovCL(m_lin, cluster = df_lin$Drug.name)

# Coefficient test with clustered SEs
coeftest(m_lin, vcov = cov_drug)

# CIs and p-values from clustered SEs (2 decimals)
modelsummary(
  m_lin,
  vcov = cov_drug,
  estimate  = "{estimate}", 
  statistic = "[{conf.low}, {conf.high}]  p={p.value}",
  fmt = 2  )

sessionInfo()