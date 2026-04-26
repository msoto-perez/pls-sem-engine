############################################################
# VALIDATION SCRIPT – Numerical comparison with plspm
# Demonstrates numerical equivalence between the proposed 
# base R PLS-SEM engine and the established plspm package.
############################################################

# =================================================================
# 1. Environment Setup (Load both engines)
# =================================================================
if (!requireNamespace("plspm", quietly = TRUE)) install.packages("plspm")
library(plspm)

devtools::install_github("msoto-perez/pls-sem-engine", force = TRUE)
library(plssemengine)

# =====================
# 2. SIMULATED DATA
# =====================
set.seed(123)
n <- 300

Service_Quality <- rnorm(n)
Customer_Satisfaction <- 0.6 * Service_Quality + rnorm(n, sd = 0.6)
Customer_Loyalty <- 0.55 * Customer_Satisfaction + 0.25 * Service_Quality + rnorm(n, sd = 0.6)

latent_to_item <- function(latent, loading) {
  x <- loading * latent + rnorm(length(latent), sd = sqrt(1 - loading^2))
  x <- scale(x)
  as.numeric(cut(x, breaks = quantile(x, probs = seq(0, 1, length.out = 8)), labels = 1:7, include.lowest = TRUE))
}

simulated_data <- data.frame(
  SQ1 = latent_to_item(Service_Quality, 0.82), SQ2 = latent_to_item(Service_Quality, 0.78), SQ3 = latent_to_item(Service_Quality, 0.74),
  CS1 = latent_to_item(Customer_Satisfaction, 0.80), CS2 = latent_to_item(Customer_Satisfaction, 0.76), CS3 = latent_to_item(Customer_Satisfaction, 0.72),
  CL1 = latent_to_item(Customer_Loyalty, 0.81), CL2 = latent_to_item(Customer_Loyalty, 0.77), CL3 = latent_to_item(Customer_Loyalty, 0.73)
)

# =====================
# 3. plspm EXECUTION
# =====================
path_matrix <- rbind(
  Service_Quality = c(0, 0, 0),
  Customer_Satisfaction = c(1, 0, 0),
  Customer_Loyalty = c(1, 1, 0)
)
colnames(path_matrix) <- rownames(path_matrix)

blocks <- list(c("SQ1","SQ2","SQ3"), c("CS1","CS2","CS3"), c("CL1","CL2","CL3"))
modes <- c("A","A","A")

plspm_model <- plspm(simulated_data, path_matrix, blocks, modes = modes)

# =====================
# 4. PROPOSED ENGINE EXECUTION
# =====================
measurement_model <- list(
  Service_Quality = c("SQ1", "SQ2", "SQ3"),
  Customer_Satisfaction = c("CS1", "CS2", "CS3"),
  Customer_Loyalty = c("CL1", "CL2", "CL3")
)
structural_model <- list(
  Customer_Satisfaction ~ Service_Quality,
  Customer_Loyalty ~ Customer_Satisfaction + Service_Quality
)

# Run core engine (no bootstrap needed for point estimates comparison)
proposed_model <- pls_engine(simulated_data, measurement_model, structural_model)

# =====================
# 5. COMPARISON REPORT
# =====================
# Extract path coefficients
plspm_paths <- plspm_model$path_coefs
prop_paths <- proposed_model$paths

comparison <- data.frame(
  Relationship = c(
    "Service_Quality -> Customer_Satisfaction",
    "Service_Quality -> Customer_Loyalty",
    "Customer_Satisfaction -> Customer_Loyalty",
    "R2 (Customer_Satisfaction)",
    "R2 (Customer_Loyalty)"
  ),
  plspm = c(
    plspm_paths["Customer_Satisfaction", "Service_Quality"],
    plspm_paths["Customer_Loyalty", "Service_Quality"],
    plspm_paths["Customer_Loyalty", "Customer_Satisfaction"],
    plspm_model$inner_model["Customer_Satisfaction", "R2"],
    plspm_model$inner_model["Customer_Loyalty", "R2"]
  ),
  Proposed_Engine = c(
    prop_paths$Customer_Satisfaction["Service_Quality"],
    prop_paths$Customer_Loyalty["Service_Quality"],
    prop_paths$Customer_Loyalty["Customer_Satisfaction"],
    proposed_model$r2["Customer_Satisfaction"],
    proposed_model$r2["Customer_Loyalty"]
  )
)

comparison$Absolute_Difference <- abs(comparison$plspm - comparison$Proposed_Engine)

cat("\n=================================================================\n")
cat(" NUMERICAL COMPARISON: plspm vs Proposed Engine\n")
cat("=================================================================\n\n")
print(comparison, row.names = FALSE, digits = 4)
cat("\n=================================================================\n")

