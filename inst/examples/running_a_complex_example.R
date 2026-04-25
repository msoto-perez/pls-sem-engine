# ==============================================================================
# Complex Example: Extended Technology Acceptance Model (TAM)
# Demonstrates the scalability of the plssemengine with a 5-construct model
# and dynamic plotting capabilities.
# ==============================================================================

source("pls_engine.R")

# =====================
# 1. SIMULATED DATA (N = 400)
# =====================
set.seed(456)
n <- 400

# Core structural relationships
System_Quality <- rnorm(n)
Info_Quality   <- rnorm(n)

Perceived_Ease_of_Use <- 0.45 * System_Quality + 0.35 * Info_Quality + rnorm(n, sd = 0.5)
Perceived_Usefulness  <- 0.50 * Perceived_Ease_of_Use + 0.40 * Info_Quality + rnorm(n, sd = 0.5)
Intention_to_Use      <- 0.55 * Perceived_Usefulness + 0.30 * Perceived_Ease_of_Use + rnorm(n, sd = 0.4)

# Helper function to generate Likert-scale items
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

# Generate dataset
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
# 2. MODEL SPECIFICATION
# =====================
measurement_model_tam <- list(
  System_Quality        = c("SQ1", "SQ2", "SQ3"),
  Info_Quality          = c("IQ1", "IQ2", "IQ3"),
  Perceived_Ease_of_Use = c("PEOU1", "PEOU2", "PEOU3"),
  Perceived_Usefulness  = c("PU1", "PU2", "PU3"),
  Intention_to_Use      = c("ITU1", "ITU2", "ITU3")
)

# Complex structural paths (multiple mediations)
structural_model_tam <- list(
  Perceived_Ease_of_Use ~ System_Quality + Info_Quality,
  Perceived_Usefulness  ~ Perceived_Ease_of_Use + Info_Quality,
  Intention_to_Use      ~ Perceived_Usefulness + Perceived_Ease_of_Use
)

# =====================
# 3. FINAL EXECUTION 
# =====================
cat("Running complex PLS-SEM model...\n")
model_tam <- pls_sem(
  data = data_tam,
  measurement_model = measurement_model_tam,
  structural_model = structural_model_tam,
  k = 5,
  nboot = 200 # Set to 200 for faster example execution, use 500-1000 for papers
)

# =====================
# 4. RESULTS & PLOTS
# =====================
model_tam$tables$table1
model_tam$tables$table2
model_tam$tables$table3
model_tam$tables$table4
model_tam$tables$table5
model_tam$tables$cmb 

# Further diagnostics
get_indirect_effects(model_tam)
get_references()

# Test the dynamic plotting with 5 constructs
plot_structural_model(structural_model_tam)
plot_structural_model(structural_model_tam, save_plot = TRUE, file_name = "conceptual model.png")

plot_model_results(model_tam)
plot_model_results(model_tam, 
                   box_height = 0.5,
                   box_width = 1.3,  
                   cex_node = 0.85,
                   cex_r2 = 0.50,
                   cex_beta = 0.80, 
                   save_plot = TRUE, 
                   file_name = "model_with_results.png")
