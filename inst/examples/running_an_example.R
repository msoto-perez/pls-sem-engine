#################################################################
# ILLUSTRATIVE EXAMPLE: Full PLS-SEM Workflow
# Package: PLSsemEngine (v1.2.0)
# Purpose: Demonstrates transparency, modularity, and diagnostics
#################################################################

devtools::install_github("msoto-perez/PLSsemEngine", force = TRUE)
library(PLSsemEngine)

# =============================================================
# 1. DATA GENERATION (Reproducible Synthetic Dataset)
# =============================================================
# Simulates a Service Marketing framework (N = 300) [cite: 152, 153]
set.seed(123)
n <- 300

# Theoretical latent constructs
Service_Quality <- rnorm(n)
Customer_Satisfaction <- 0.6 * Service_Quality + rnorm(n, sd = 0.6)
Customer_Loyalty <- 0.55 * Customer_Satisfaction + 0.25 * Service_Quality + rnorm(n, sd = 0.6)

# Transformation function: Latent scores to 7-point Likert items [cite: 153]
latent_to_item <- function(latent, loading) {
  x <- loading * latent + rnorm(length(latent), sd = sqrt(1 - loading^2))
  x <- scale(x)
  as.numeric(cut(x, breaks = quantile(x, probs = seq(0, 1, length.out = 8)), 
                 labels = 1:7, include.lowest = TRUE))
}

# Construct the dataframe SQ1-SQ3, CS1-CS3, CL1-CL3 [cite: 153]
simulated_data <- data.frame(
  SQ1 = latent_to_item(Service_Quality, 0.82), SQ2 = latent_to_item(Service_Quality, 0.78), SQ3 = latent_to_item(Service_Quality, 0.74),
  CS1 = latent_to_item(Customer_Satisfaction, 0.80), CS2 = latent_to_item(Customer_Satisfaction, 0.76), CS3 = latent_to_item(Customer_Satisfaction, 0.72),
  CL1 = latent_to_item(Customer_Loyalty, 0.81), CL2 = latent_to_item(Customer_Loyalty, 0.77), CL3 = latent_to_item(Customer_Loyalty, 0.73)
)

# =============================================================
# 2. MODEL SPECIFICATION (Native R Interface)
# =============================================================
# Defining the reflective measurement model [cite: 161, 172]
measurement_model <- list(
  Service_Quality = c("SQ1", "SQ2", "SQ3"),
  Customer_Satisfaction = c("CS1", "CS2", "CS3"),
  Customer_Loyalty = c("CL1", "CL2", "CL3")
)

# Defining the structural paths using R formulas [cite: 166, 172]
structural_model <- list(
  Customer_Satisfaction ~ Service_Quality,
  Customer_Loyalty ~ Customer_Satisfaction + Service_Quality
)

# =============================================================
# 3. MODEL EXECUTION
# =============================================================
# High-level wrapper coordinates estimation, bootstrap, and prediction [cite: 179]
# nboot = 500 is used for stable inference as per SoftwareX standards [cite: 84]
model <- pls_sem(
  data = simulated_data,
  measurement_model = measurement_model,
  structural_model = structural_model,
  nboot = 500,
  k = 5
)

# =============================================================
# 4. COMPREHENSIVE OUTPUTS (Aligned with Manuscript Tables)
# =============================================================

# 4.1 Measurement Model Assessment (Table 1) [cite: 58]
# Includes Loadings, Composite Reliability (CR), AVE, and R2 [cite: 7]
print(model$measurement_model)

# 4.2 Discriminant Validity (Table 3) [cite: 99, 133]
# Includes both HTMT and HTMT2 for congeneric models [cite: 292]
print(model$discriminant_validity$HTMT)
print(model$discriminant_validity$HTMT2)

# 4.3 Structural Model and Inference (Table 4) [cite: 58]
# Includes Path Coefficients, Bootstrap CIs, and f2 Effect Sizes [cite: 7]
print(model$structural_model)

# 4.4 Predictive Relevance (Table 5) [cite: 88, 92]
# PLSpredict results: RMSE, MAE, and Q2_predict [cite: 7]
print(model$predictive_relevance)

# =============================================================
# 5. ADVANCED DIAGNOSTICS & BRIDGES (Reviewer Requests)
# =============================================================

# 5.1 Common Method Bias (CMB) [cite: 99, 292]
# Full Collinearity VIF following Kock (2015)
print(model$diagnostics$common_method_bias)

# 5.2 Global Model Fit [cite: 292]
# Indices: SRMR, d_ULS, and d_G
print(model$diagnostics$global_fit)

# 5.3 Indirect Effects (Mediation Analysis) [cite: 85]
# Computed as the product of corresponding paths
PLSsemEngine:::get_indirect_effects(model)

# 5.4 Methodological Bridge: CB-SEM Syntax [cite: 34, 287]
# Automatically translates PLS specification to lavaan code
export_lavaan_syntax(measurement_model, structural_model)

# 5.5 Optional Structured Interpretive Layer [cite: 36, 288]
# Contextualizes metrics against established literature thresholds
interpret_model(model)

# =============================================================
# 6. VISUALIZATION [cite: 283]
# =============================================================
# 6.1 Plot conceptual structural model
plot_structural_model(structural_model)

# 6.2 Plot final model with results (Beta values and R2)
plot_model_results(model, show_r2 = TRUE)

# 6.3 Export results for external reporting
export_scores(model, "latent_scores_v120.csv") [cite: 100, 107]
