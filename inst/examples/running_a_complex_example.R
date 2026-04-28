################################################################################
# Complex Example: Extended Technology Acceptance Model (TAM)
# Version: 1.2.0
# Description: Demonstrates scalability with a 5-construct model, 
#              automated circular plotting, and mediation analysis.
################################################################################

# Load the package
devtools::install_github("msoto-perez/PLSsemEngine", force = TRUE)
library(PLSsemEngine)

# =====================
# 1. SIMULATED DATA (N = 400)
# =====================
set.seed(456)
n <- 400

# Core structural relationships (Latent Variables)
System_Quality <- rnorm(n)
Info_Quality   <- rnorm(n)

Perceived_Ease_of_Use <- 0.45 * System_Quality + 0.35 * Info_Quality + rnorm(n, sd = 0.5)
Perceived_Usefulness  <- 0.50 * Perceived_Ease_of_Use + 0.40 * Info_Quality + rnorm(n, sd = 0.5)
Intention_to_Use      <- 0.55 * Perceived_Usefulness + 0.30 * Perceived_Ease_of_Use + rnorm(n, sd = 0.4)

# Helper function to generate Likert-scale items (1-7)
latent_to_item <- function(latent, loading) {
  x <- loading * latent + rnorm(length(latent), sd = sqrt(1 - loading^2))
  x <- scale(x)
  as.numeric(cut(
    x,
    breaks = quantile(x, probs = seq(0, 1, length.out = 8)),
    labels = 1:7,
    include.lowest = TRUE
  ))
}

# Generate dataset with item names corresponding to abbreviations
data_tam <- data.frame(
  SQ1 = latent_to_item(System_Quality, 0.85),
  SQ2 = latent_to_item(System_Quality, 0.80),
  SQ3 = latent_to_item(System_Quality, 0.75),
  
  IQ1 = latent_to_item(Info_Quality, 0.82),
  IQ2 = latent_to_item(Info_Quality, 0.78),
  IQ3 = latent_to_item(Info_Quality, 0.76),
  
  PEOU1 = latent_to_item(Perceived_Ease_of_Use, 0.88),
  PEOU2 = latent_to_item(Perceived_Ease_of_Use, 0.84),
  PEOU3 = latent_to_item(Perceived_Ease_of_Use, 0.79),
  
  PU1 = latent_to_item(Perceived_Usefulness, 0.86),
  PU2 = latent_to_item(Perceived_Usefulness, 0.82),
  PU3 = latent_to_item(Perceived_Usefulness, 0.80),
  
  ITU1 = latent_to_item(Intention_to_Use, 0.89),
  ITU2 = latent_to_item(Intention_to_Use, 0.85),
  ITU3 = latent_to_item(Intention_to_Use, 0.81)
)

# =====================
# 2. MODEL SPECIFICATION (Using Abbreviations)
# =====================
# Abbreviations (SYST, INFO, PEOU, PU, ITU) ensure clean circular plots
measurement_model_tam <- list(
  SYST = c("SQ1", "SQ2", "SQ3"),
  INFO = c("IQ1", "IQ2", "IQ3"),
  PEOU = c("PEOU1", "PEOU2", "PEOU3"),
  PU   = c("PU1", "PU2", "PU3"),
  ITU  = c("ITU1", "ITU2", "ITU3")
)

# Define structural paths using the new abbreviations
structural_model_tam <- list(
  PEOU ~ SYST + INFO,
  PU   ~ PEOU + INFO,
  ITU  ~ PU + PEOU
)

# =====================
# 3. MODEL EXECUTION
# =====================
cat("Executing PLS-SEM for Complex TAM Model...\n")
model_tam <- pls_sem(
  data = data_tam,
  measurement_model = measurement_model_tam,
  structural_model = structural_model_tam,
  k = 5,
  nboot = 200 # Faster for examples; use 500+ for final research
)

# =====================
# 4. COMPREHENSIVE RESULTS
# =====================
# Display main results with descriptive names (v1.2.0 standard)
print(model_tam$measurement_model)      # Loadings, Reliability, and R2
print(model_tam$discriminant_validity)  # HTMT and HTMT2 matrices
print(model_tam$structural_model)       # Path Coefficients and f2
print(model_tam$predictive_relevance)   # PLSpredict (RMSE, Q2_predict)

# Multicollinearity & Global Fit
print(model_tam$diagnostics$common_method_bias)
print(model_tam$diagnostics$global_fit)

# Indirect Effects (Mediation analysis)
get_indirect_effects(model_tam)

# Narrative Interpretation Layer
interpret_model(model_tam)

# =====================
# 5. VISUALIZATION (Automated Circular Layout)
# =====================
# Using abbreviations allows for standard box sizes without overlapping
plot_model_results(model_tam, 
                   show_r2 = FALSE,  
                   box_width = 1.1, 
                   box_height = 0.3)

# Save high-resolution circular plot for external reporting
plot_model_results(model_tam, 
                   layout     = NULL, 
                   save_plot  = TRUE, 
                   file_name  = "TAM_Circular_Results.png")

cat("\nAnalysis complete. Results and plots are ready.\n")
