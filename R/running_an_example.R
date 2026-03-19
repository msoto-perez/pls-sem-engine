# Example script to reproduce a full PLS-SEM workflow
# using simulated data (N = 300).

source("pls_engine.R")

# =====================
# SIMULATED DATA
# =====================
set.seed(123)
n <- 300

Service_Quality <- rnorm(n)

Customer_Satisfaction <-
  0.6 * Service_Quality +
  rnorm(n, sd = 0.6)

Customer_Loyalty <-
  0.55 * Customer_Satisfaction +
  0.25 * Service_Quality +
  rnorm(n, sd = 0.6)

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

simulated_data <- data.frame(
  SQ1 = latent_to_item(Service_Quality, 0.82),
  SQ2 = latent_to_item(Service_Quality, 0.78),
  SQ3 = latent_to_item(Service_Quality, 0.74),
  CS1 = latent_to_item(Customer_Satisfaction, 0.80),
  CS2 = latent_to_item(Customer_Satisfaction, 0.76),
  CS3 = latent_to_item(Customer_Satisfaction, 0.72),
  CL1 = latent_to_item(Customer_Loyalty, 0.81),
  CL2 = latent_to_item(Customer_Loyalty, 0.77),
  CL3 = latent_to_item(Customer_Loyalty, 0.73)
)

# =====================
# MODEL
# =====================
measurement_model <- list(
  Service_Quality = c("SQ1","SQ2","SQ3"),
  Customer_Satisfaction = c("CS1","CS2","CS3"),
  Customer_Loyalty = c("CL1","CL2","CL3")
)

structural_model <- list(
  Customer_Satisfaction ~ Service_Quality,
  Customer_Loyalty ~ Customer_Satisfaction + Service_Quality
)

# =====================
# FINAL EXECUTION 
# =====================
model <- pls_sem(
  data = simulated_data,
  measurement_model = measurement_model,
  structural_model = structural_model,
  k = 5,
  nboot = 25
)

# Tables paper-ready
model$tables$table1
model$tables$table2
model$tables$table3
model$tables$table4
model$tables$table5
model$tables$cmb 

# Further diagnostics
get_indirect_effects(model)
export_scores(model, "scores.csv")
get_references()
htmt_item_diagnostics(model)

# Warnings
model$warnings

# Warnings
model$warnings

# Plot structural model
plot_structural_model(structural_model)

# Plot structural model (with results)
plot_model_results(model)

