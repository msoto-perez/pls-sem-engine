# plssemengine: A Transparent PLS-SEM Engine in Base R

**Version:** 1.1.0 (2026-03-19)

**Repository:** https://github.com/msoto-perez/pls-sem-engine

**DOI (Zenodo):** https://doi.org/10.5281/zenodo.19120651

**License:** MIT

---

## 1. Purpose
`plssemengine` provides a transparent, modular, and reproducible implementation of Partial Least Squares Structural Equation Modeling (PLS-SEM). It is specifically designed for reflective measurement models (Mode A).

The software prioritizes:
* **Algorithmic transparency** via pure base R matrix operations.
* **Explicit analytical control** without hidden heuristics.
* **Modular separation** between estimation, inference, and prediction.

---

## 2. Computational Workflow
The engine implements a standardized PLS-SEM pipeline:
* Data standardization.
* Iterative Mode A estimation (Factorial/Centroid schemes).
* Measurement model evaluation (Loadings, CR, AVE, HTMT).
* Structural model estimation via OLS.
* Bootstrap inference for path coefficients.
* k-fold cross-validation (PLSpredict).

**Note:** Deterministic sign alignment is available but disabled by default to preserve bootstrap distribution integrity.

---

## 3. What's New in v1.1.0
* **Effect Size ($f^2$):** Automatically calculated and integrated into Table 4.
* **Visualization:** Native plotting functions `plot_model_results()` and `plot_structural_model()`.
* **Standardized Variance:** Composite variance is fixed to 1 following Wold (1982).
* **Inner Schemes:** Support for multiple inner weighting schemes.

---

## 4. Minimal Example (Copy-paste ready)

```r
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
---

## 5. Citation
If you use this software in your research, please cite the accompanying paper:

Soto-Perez, M. (2026). A transparent PLS-SEM engine for reflective measurement models in R. *SoftwareX*. (Manuscript under review).

**Software archive:**
Soto-Perez, M. (2026). *plssemengine (Version 1.1.0)*. Zenodo. https://doi.org/10.5281/zenodo.19120651

---

## 6. Contact
**Dr. M. Soto-Perez** 

Email: msoto@up.edu.mx  
Universidad Panamericana, Mexico

