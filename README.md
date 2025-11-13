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
# Core packages
library(tidyverse)
library(broom)
library(forcats)
library(lubridate)

# Statistical modeling
library(car)
library(sandwich)
library(lmtest)
library(survival)

# Table and summary formatting
library(modelsummary)
library(gt)
library(janitor)

# For visualization and data cleaning
library(ggplot2)
library(readxl)
```

## Repository Organization
```
├── /data/                 # Publicly available curated datasets
│   ├── fda_drug_approvals_2011_2023.csv
│   ├── pediatric_labeling_additions.csv
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


