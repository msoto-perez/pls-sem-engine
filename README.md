\# PLS-SEM Engine: A Transparent Reflective Engine in Base R 📊



\[!\[DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19120651.svg)](https://doi.org/10.5281/zenodo.19120651)

\[cite\_start]\*\*Version:\*\* 1.1.0 (2026-03-24) 



\[cite\_start]\*\*plssemengine\*\* provides a transparent, modular, and reproducible implementation of Partial Least Squares Structural Equation Modeling (PLS-SEM), specifically designed for \*\*reflective measurement models (Mode A)\*\*. \[cite: 4, 118]



\## 🌟 Purpose and Philosophy

The software prioritizes:

\* \[cite\_start]\*\*Algorithmic Transparency:\*\* Implementation via pure base R matrix operations. \[cite: 14, 35, 112]

\* \[cite\_start]\*\*Explicit Analytical Control:\*\* No hidden heuristics or automatic re-specifications. \[cite: 7, 34, 119]

\* \[cite\_start]\*\*Modular Architecture:\*\* Clear separation between estimation, inference, and prediction. \[cite: 5, 31, 45]



\## ⚙️ Computational Workflow

\[cite\_start]The engine follows a standardized PLS-SEM pipeline\[cite: 46, 56]:

1\. \[cite\_start]\*\*Data Standardization:\*\* Handles centered and standardized scales. \[cite: 48, 58]

2\. \[cite\_start]\*\*Iterative Mode A Estimation:\*\* Factorial and Centroid weighting schemes. \[cite: 62, 63]

3\. \[cite\_start]\*\*Measurement Evaluation:\*\* Loadings, CR, AVE, and HTMT. \[cite: 6, 89]

4\. \[cite\_start]\*\*Structural Estimation:\*\* Path coefficients via OLS on latent scores. \[cite: 6, 69]

5\. \[cite\_start]\*\*Inference:\*\* Bootstrap resampling for structural significance. \[cite: 6, 70]

6\. \[cite\_start]\*\*Prediction:\*\* Strict k-fold cross-validation (PLSpredict). \[cite: 6, 71, 78]



> \[cite\_start]\*\*Note:\*\* Deterministic sign alignment is available but disabled by default to preserve bootstrap distribution integrity. \[cite: 65, 66]



\## 🚀 What's New in V1.1.0

\* \[cite\_start]Added \*\*f-squared ($f^2$)\*\* effect size calculation. 

\* \[cite\_start]Implemented native \*\*plotting functions\*\* (`plot\_model\_results`). 

\* \[cite\_start]Fixed composite variance to 1 (Wold, 1982). \[cite: 63]

\* \[cite\_start]Added support for multiple inner weighting schemes. \[cite: 62]



\## 🛠️ Minimal Example

```r

\# Load the engine

source("R/pls\_engine.R")



\# 1. Generate synthetic data

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



