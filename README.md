# pls-sem-engine

A transparent PLS-SEM engine for reflective measurement models (Mode A)
implemented in base R.

**Current version:** v1.0.4  
**Repository:** https://github.com/msoto-perez/pls-sem-engine  
**DOI (Zenodo archive):** https://doi.org/10.5281/zenodo.18794941 
**License:** MIT

------------------------------------------------------------------------

## 1. Purpose

`pls-sem-engine` provides a transparent and reproducible implementation
of Partial Least Squares Structural Equation Modeling (PLS-SEM)
restricted to reflective measurement models (Mode A).

The software was developed to prioritize:

-   Algorithmic transparency\
-   Explicit analytical control\
-   Reproducibility\
-   Clear separation between estimation, inference, and prediction

The engine avoids automated model re-specification, heuristic threshold
classifications, and interpretative flags.

------------------------------------------------------------------------

## 2. Computational Workflow

The engine follows a modular workflow:

1.  Data standardization\
2.  Iterative Mode A latent variable score estimation\
3.  Measurement model evaluation\
4.  Structural model estimation via ordinary least squares\
5.  Bootstrap inference (observation-level resampling)\
6.  k-fold cross-validation for predictive assessment

All computations are implemented using base R matrix operations.

Deterministic sign alignment is applied after convergence to ensure
reproducible latent variable orientation.

------------------------------------------------------------------------

## 3. Scope

### Supported

-   Reflective measurement models (Mode A)\
-   Bootstrap confidence intervals for structural coefficients\
-   HTMT discriminant validity assessment\
-   Full collinearity variance inflation factors (VIF)\
-   Predictive evaluation via k-fold cross-validation\
-   Paper-ready tabular outputs

### Not Supported (by design)

-   Formative measurement models (Mode B)\
-   Automatic mediation or moderation testing\
-   Automatic model re-specification\
-   Multigroup analysis\
-   Measurement invariance testing

These restrictions are intentional to preserve methodological clarity
and researcher control.

------------------------------------------------------------------------

## 4. Installation

Clone the repository:

    git clone https://github.com/msoto-perez/pls-sem-engine

The engine requires only base R and has no external dependencies.

------------------------------------------------------------------------

## 5. Minimal Example

To test the engine immediately, you can run this self-contained example with simulated data:

```r
source("pls_engine_v1.0.4.R")

# 1. Generate dummy data (n=100)
set.seed(123)
data <- data.frame(
  SQ1 = rnorm(100), SQ2 = rnorm(100), SQ3 = rnorm(100),
  CS1 = rnorm(100), CS2 = rnorm(100), CS3 = rnorm(100),
  CL1 = rnorm(100), CL2 = rnorm(100), CL3 = rnorm(100)
)

# 2. Define Measurement Model (List of blocks)
measurement_model <- list(
  Service_Quality       = c("SQ1", "SQ2", "SQ3"),
  Customer_Satisfaction = c("CS1", "CS2", "CS3"),
  Customer_Loyalty      = c("CL1", "CL2", "CL3")
)

# 3. Define Structural Model (Formula syntax)
structural_model <- list(
  Customer_Satisfaction ~ Service_Quality,
  Customer_Loyalty      ~ Customer_Satisfaction + Service_Quality
)

# 4. Run PLS-SEM Analysis
model <- pls_sem(
  data = data,
  measurement_model = measurement_model,
  structural_model = structural_model,
  nboot = 500,
  k = 5
)

# 5. Inspect Results
print(model$path_coefficients)
print(model$reliability)

------------------------------------------------------------------------

## 6. Numerical Validation

The engine was numerically validated against the established `plspm`
implementation.

Absolute differences in path coefficients and R² values were
consistently below 0.003, confirming computational consistency.

------------------------------------------------------------------------

## 7. Reproducibility

-   Version controlled via Git\
-   Archived via Zenodo DOI\
-   Deterministic sign alignment implemented\
-   Explicit seed control for bootstrap and cross-validation

This ensures traceability between published results and the exact
software version used.

------------------------------------------------------------------------

## 8. Citation

If you use this software, please cite:

Soto-Pérez, M. (2026).  
pls-sem-engine (Version 1.0.4) [Software].  
Zenodo. https://doi.org/10.5281/zenodo.18794941

------------------------------------------------------------------------

## 9. Contact

msoto@up.edu.mx
