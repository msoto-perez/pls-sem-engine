# PLS-SEM Engine: A Transparent Reflective Engine in Base R 📊

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19120651.svg)](https://doi.org/10.5281/zenodo.19120651)

**Version:** 1.1.0 (2026-03-24)

**plssemengine** provides a transparent, modular, and reproducible implementation of Partial Least Squares Structural Equation Modeling (PLS-SEM), specifically designed for **reflective measurement models (Mode A)**[cite: 4, 118].

## 🌟 Purpose and Philosophy

The software prioritizes:
* [cite_start]**Algorithmic Transparency:** Implementation via pure base R matrix operations[cite: 14, 35, 112].
* [cite_start]**Explicit Analytical Control:** No hidden heuristics or automatic re-specifications[cite: 7, 34, 119].
* [cite_start]**Modular Architecture:** Clear separation between estimation, inference, and prediction[cite: 5, 31, 45].

## ⚙️ Computational Workflow

[cite_start]The engine follows a standardized PLS-SEM pipeline[cite: 46, 56]:

1. [cite_start]**Data Standardization:** Handles centered and standardized scales[cite: 48, 58].
2. [cite_start]**Iterative Mode A Estimation:** Factorial and Centroid weighting schemes[cite: 62, 63].
3. [cite_start]**Measurement Evaluation:** Loadings, CR, AVE, and HTMT[cite: 6, 89].
4. [cite_start]**Structural Estimation:** Path coefficients via OLS on latent scores[cite: 6, 69].
5. [cite_start]**Inference:** Bootstrap resampling for structural significance[cite: 6, 70].
6. [cite_start]**Prediction:** Strict k-fold cross-validation (PLSpredict)[cite: 6, 71, 78].

> [cite_start]**Note:** Deterministic sign alignment is available but disabled by default to preserve bootstrap distribution integrity[cite: 65, 66].

## 🚀 What's New in V1.1.0
* Added **f-squared ($f^2$)** effect size calculation.
* Implemented native **plotting functions** (`plot_model_results`).
* [cite_start]Fixed composite variance to 1 (Wold, 1982)[cite: 63].
* [cite_start]Added support for multiple inner weighting schemes[cite: 62].

## 🛠️ Minimal Example
```r
# Load the engine
source("R/pls_engine.R")

# 1. Generate data
set.seed(123)
data <- data.frame(
  SQ1=rnorm(100), SQ2=rnorm(100), SQ3=rnorm(100),
  CS1=rnorm(100), CS2=rnorm(100), CS3=rnorm(100),
  CL1=rnorm(100), CL2=rnorm(100), CL3=rnorm(100)
)

# 2. Define Models
mm <- list(
  Quality = c("SQ1", "SQ2", "SQ3"),
  Satisfaction = c("CS1", "CS2", "CS3"),
  Loyalty = c("CL1", "CL2", "CL3")
)

sm <- list(
  Satisfaction ~ Quality,
  Loyalty ~ Satisfaction + Quality
)

# 3. Run Analysis
model <- pls_sem(data=data, measurement_model=mm, structural_model=sm)

# 4. Results & Plots
model$tables$table4  # Structural Paths & f2
plot_model_results(model)

```



\## 📖 Citation

If you use this software, please cite:



\*\*Manuscript:\*\*

\*\*Soto-Perez, M. (2026).\*\* \*A transparent PLS-SEM engine for reflective measurement models in R\*. SoftwareX. (Under review) \[cite\_start]\[cite: 1].



\*\*Software Archive:\*\*

\*\*Soto-Perez, M. (2026).\*\* \*plssemengine (Version 1.1.0)\*. Zenodo. \[cite\_start]\[https://doi.org/10.5281/zenodo.19120651](https://doi.org/10.5281/zenodo.19120651).



\## ✉️ Contact

\*\*Dr. \[cite\_start]M. Soto-Perez\*\* Email: \[msoto@up.edu.mx](mailto:msoto@up.edu.mx)   

Universidad Panamericana, Mexico



