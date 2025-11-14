# Pediatric Studies and Labeling Additions Required by the US FDA for Novel Drug Approvals, 2011–2023

This repository contains the R scripts and analytical workflow used in the study *“Pediatric Studies and Labeling Additions Required by the US FDA for Novel Drug Approvals, 2011–2023.”* The code reproduces all main analyses presented in the manuscript, including generation of **Tables 1–3** and the **supplementary materials**. The repository includes data preparation, statistical analyses (ANOVA, logistic regression, and Cox proportional hazards models), and figure/table formatting steps used for publication. All analyses were conducted in R, and the repository is organized for reproducibility and transparency.

---

## Datasets and Authorization

This repository includes the curated dataset compiled for the analysis, derived from the **FDA’s Postmarketing Requirements and Commitments Database**.  
No patient-level or confidential data are included.  
For additional details about dataset generation, please refer to the *Methods* section of the manuscript.

---

## Dependencies

All analyses were performed in **R (≥4.2.0)**.  
The following R packages are required to reproduce the results:

```r
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
```

## Repository Organization
```
├── /data/                 # Publicly available curated datasets
│   ├── fda_drug_approvals_2011_2023.txt
│   ├── pediatric_labeling_additions.txt
│   └── 
├── /code/                 # R scripts for all analyses
│   ├── table1_linear_regression.R
│   ├── table2_logistic_models.R
│   ├── table3_cox_models.R
│   └── supplementary_analysis.R
└── README.md
```

## Citation
If you use this dataset or code, please cite:
[ CITATION AND LINK TO THE PUBLICATION TO BE ADDED ]


## License
Licensed under the [Apache License, Version 2.0](http://www.apache.org/licenses/LICENSE-2.0).  


