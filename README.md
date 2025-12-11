# Pediatric Studies and Labeling Additions Required by the US FDA for Novel Drugs Approved from 2011-2023: A retrospective cohort study

This repository contains the R scripts and analytical workflow used in the study *“Pediatric Studies and Labeling Additions Required by the US FDA for Novel Drugs Approved from 2011-2023: A retrospective cohort study.”* The code reproduces all main analyses presented in the manuscript, including generation of **Tables 1–3, Figures 1-2** and the **supplementary table**. The repository includes data preparation, statistical analyses (ANOVA, logistic regression, and Cox proportional hazards models), and figure/table formatting steps used for publication. All analyses were conducted in R, and the repository is organized for reproducibility and transparency.

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
├── /data/                 
│   ├── fda_drug_approvals_2011_2023.txt
│   ├── pediatric_labeling_additions.txt
│   └── 
├── /code/                
│   ├── Table1.R
│   ├── Table2.R
│   ├── Table3.R
│   ├── Figure1andFigure2.R
│   └── SupplementaryTable.R
└── README.md
```

## Citation
If you use this dataset or code, please cite:
[ CITATION AND LINK TO THE PUBLICATION TO BE ADDED ]


## License
Licensed under the [Apache License, Version 2.0](http://www.apache.org/licenses/LICENSE-2.0).  


