PLSSEMENGINE: A TRANSPARENT PLS-SEM ENGINE IN BASE R



Version: 1.1.0 (2026-03-24)

Repository: https://github.com/msoto-perez/pls-sem-engine

DOI (Zenodo): https://doi.org/10.5281/zenodo.19120651

License: MIT



\--------------------------------------------------

1\. PURPOSE

\--------------------------------------------------

plssemengine provides a transparent, modular, and reproducible

implementation of Partial Least Squares Structural Equation

Modeling (PLS-SEM). It is specifically designed for

reflective measurement models (Mode A).



The software prioritizes:

\- Algorithmic transparency via pure base R matrix operations.

\- Explicit analytical control without hidden heuristics.

\- Modular separation between estimation, inference, and prediction.



\--------------------------------------------------

2\. COMPUTATIONAL WORKFLOW

\--------------------------------------------------

The engine implements a standardized PLS-SEM pipeline:

\- Data standardization.

\- Iterative Mode A estimation (Factorial/Centroid schemes).

\- Measurement model evaluation (Loadings, CR, AVE, HTMT).

\- Structural model estimation via OLS.

\- Bootstrap inference for path coefficients.

\- k-fold cross-validation (PLSpredict).



Note: Deterministic sign alignment is available but disabled by

default to preserve bootstrap distribution integrity.



\--------------------------------------------------

3\. WHAT'S NEW IN V1.1.0

\--------------------------------------------------

\- Added f-squared (f2) effect size calculation.

\- Implemented native plotting functions (plot\_model\_results).

\- Fixed composite variance to 1 (Wold, 1982).

\- Added support for multiple inner weighting schemes.



\--------------------------------------------------

4\. MINIMAL EXAMPLE (Copy-paste ready)

\--------------------------------------------------



source("pls\_engine.R")



\# 1. Generate data

set.seed(123)

data <- data.frame(

&#x20; SQ1=rnorm(100), SQ2=rnorm(100), SQ3=rnorm(100),

&#x20; CS1=rnorm(100), CS2=rnorm(100), CS3=rnorm(100),

&#x20; CL1=rnorm(100), CL2=rnorm(100), CL3=rnorm(100)

)



\# 2. Define Models

mm <- list(

&#x20; Quality = c("SQ1", "SQ2", "SQ3"),

&#x20; Satisfaction = c("CS1", "CS2", "CS3"),

&#x20; Loyalty = c("CL1", "CL2", "CL3")

)



sm <- list(

&#x20; Satisfaction \~ Quality,

&#x20; Loyalty \~ Satisfaction + Quality

)



\# 3. Run Analysis

model <- pls\_sem(data=data, measurement\_model=mm, structural\_model=sm)



\# 4. Results \& Plots

model$tables$table4  # Structural Paths \& f2

plot\_model\_results(model)



\--------------------------------------------------

5\. CITATION

\--------------------------------------------------

If you use this software, please cite:



Soto-Perez, M. (2026). A transparent PLS-SEM engine for

reflective measurement models in R. SoftwareX.

(Manuscript under review).



Software archive:

Soto-Perez, M. (2026). plssemengine (Version 1.1.0).

Zenodo. https://doi.org/10.5281/zenodo.19120651



\--------------------------------------------------

6\. CONTACT

\--------------------------------------------------

Dr. M. Soto-Perez

Email: msoto@up.edu.mx

Universidad Panamericana, Mexico

