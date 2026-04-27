# PLSsemEngine: A Transparent PLS-SEM Engine in Base R 📊

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19120651.svg)](https://doi.org/10.5281/zenodo.19120651)

[cite_start]**Version:** 1.2.0 (2026-04-25) 

[cite_start]**PLSsemEngine** provides a transparent, modular, and reproducible implementation of Partial Least Squares Structural Equation Modeling (PLS-SEM), specifically designed for **composite-based Mode A estimation** of reflective models[cite: 4, 5].

---

## 🌟 Purpose and Philosophy

The software prioritizes:

* [cite_start]**Algorithmic Transparency:** Implementation via pure base R matrix operations to ensure long-term stability[cite: 31, 44, 46].
* [cite_start]**Methodological Bridges:** Native integration with `lavaan` for CB-SEM/CFA cross-validation[cite: 34, 35].
* [cite_start]**Explicit Analytical Control:** No hidden heuristics or automatic re-specifications; interpretive support is optional and researcher-led[cite: 8, 28, 133].
* [cite_start]**Modular Architecture:** Clear separation between estimation, inference, and prediction components[cite: 6, 52, 59].

---

## ⚙️ Computational Workflow

[cite_start]The engine follows a standardized and inspectable PLS-SEM pipeline[cite: 60]:

1. [cite_start]**Data Standardization:** Handles mean-centered and standardized scales[cite: 55, 67].
2. [cite_start]**Iterative Mode A Estimation:** Factorial weighting scheme by default[cite: 69, 72].
3. [cite_start]**Measurement Evaluation:** Loadings, CR, AVE, HTMT, and **HTMT2**[cite: 7, 102, 128].
4. [cite_start]**Structural Estimation:** Path coefficients via OLS on latent scores with $f^2$ effect sizes[cite: 60, 205].
5. [cite_start]**Inference:** Non-parametric percentile bootstrap for structural significance[cite: 81, 82].
6. [cite_start]**Prediction:** Strict k-fold cross-validation following the `PLSpredict` protocol[cite: 88, 89].
7. [cite_start]**Model Fit:** Assessment via SRMR, $d_{ULS}$, and $d_G$[cite: 292].

> [cite_start]**Note:** Deterministic sign alignment is implemented to ensure stability across resamples and eliminate sign indeterminacy[cite: 78, 107].

---

## 🚀 What's New in V1.2.0 (Response to Reviewers)

* [cite_start]**Methodological Bridge:** Added `export_lavaan_syntax()` to translate PLS specifications for `lavaan`[cite: 286, 287].
* [cite_start]**Interpretive Layer:** Added `interpret_model()` for diagnostic guidance based on established literature without forcing mechanical decisions[cite: 288, 289].
* [cite_start]**Advanced Metrics:** Implemented **HTMT2** for congeneric models and global fit indices (SRMR, $d_{ULS}$, $d_G$)[cite: 292].
* [cite_start]**Professional Packaging:** The software is now a fully versioned R package installable via `devtools`[cite: 284].

---

## 🛠️ Minimal Example

```r
# Install and load the engine
# devtools::install_github("msoto-perez/PLSsemEngine")
library(PLSsemEngine)

# 1. Generate data
set.seed(123)
data <- data.frame(
  SQ1=rnorm(100), SQ2=rnorm(100), SQ3=rnorm(100),
  CS1=rnorm(100), CS2=rnorm(100), CS3=rnorm(100),
  CL1=rnorm(100), CL2=rnorm(100), CL3=rnorm(100)
)

# 2. Define Models using Native R structures
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

# 4. Methodological Bridge & Interpretation
export_lavaan_syntax(mm, sm)
interpret_model(model)

# 5. View Results
print(model$tables$table4)  # Structural Paths
```

📖 Citation
If you use this software, please cite:

Manuscript:
Soto-Perez, M. (2026). A transparent PLSsemEngine for composite-based Mode A estimation of reflective models in R. SoftwareX. (Under review) .

Software Archive:
Soto-Perez, M. (2026). PLSsemEngine (Version 1.2.0). Zenodo. https://doi.org/10.5281/zenodo.19120651.

✉️ Contact
Dr. M. Soto-Perez Email: msoto@up.edu.mx 

Universidad Panamericana, Mexico.

