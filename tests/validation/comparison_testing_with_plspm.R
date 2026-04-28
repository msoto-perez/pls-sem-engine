############################################################
# VALIDATION SCRIPT – Numerical comparison with plspm
# Demonstrates numerical equivalence between PLSsemEngine 
# and the established plspm package for Mode A estimation.
############################################################

# =================================================================
# 1. Environment Setup
# =================================================================
if (!requireNamespace("plspm", quietly = TRUE)) install.packages("plspm")
library(plspm)

devtools::install_github("msoto-perez/PLSsemEngine", force = TRUE)
library(PLSsemEngine)

# =====================
# 2. SIMULATED DATA
# =====================
set.seed(123)
n <- 300

# Create latent variables with defined theoretical paths
Service_Quality <- rnorm(n)
Customer_Satisfaction <- 0.6 * Service_Quality + rnorm(n, sd = 0.6)
Customer_Loyalty <- 0.55 * Customer_Satisfaction + 0.25 * Service_Quality + rnorm(n, sd = 0.6)

# Function to transform latent scores into observed 7-point Likert items
latent_to_item <- function(latent, loading) {
  x <- loading * latent + rnorm(length(latent), sd = sqrt(1 - loading^2))
  x <- scale(x)
  as.numeric(cut(x, breaks = quantile(x, probs = seq(0, 1, length.out = 8)), 
                 labels = 1:7, include.lowest = TRUE))
}

# Construct the data frame
simulated_data <- data.frame(
  SQ1 = latent_to_item(Service_Quality, 0.82), SQ2 = latent_to_item(Service_Quality, 0.78), SQ3 = latent_to_item(Service_Quality, 0.74),
  CS1 = latent_to_item(Customer_Satisfaction, 0.80), CS2 = latent_to_item(Customer_Satisfaction, 0.76), CS3 = latent_to_item(Customer_Satisfaction, 0.72),
  CL1 = latent_to_item(Customer_Loyalty, 0.81), CL2 = latent_to_item(Customer_Loyalty, 0.77), CL3 = latent_to_item(Customer_Loyalty, 0.73)
)

# =====================
# 3. EXECUTION
# =====================

# --- A. plspm Execution ---
path_matrix <- rbind(
  Service_Quality = c(0, 0, 0),
  Customer_Satisfaction = c(1, 0, 0),
  Customer_Loyalty = c(1, 1, 0)
)
colnames(path_matrix) <- rownames(path_matrix)
blocks <- list(c("SQ1","SQ2","SQ3"), c("CS1","CS2","CS3"), c("CL1","CL2","CL3"))
plspm_model <- plspm(simulated_data, path_matrix, blocks, modes = c("A","A","A"))

# --- B. PLSsemEngine Execution ---
measurement_model <- list(
  Service_Quality = c("SQ1", "SQ2", "SQ3"),
  Customer_Satisfaction = c("CS1", "CS2", "CS3"),
  Customer_Loyalty = c("CL1", "CL2", "CL3")
)
structural_model <- list(
  Customer_Satisfaction ~ Service_Quality,
  Customer_Loyalty ~ Customer_Satisfaction + Service_Quality
)

# Run internal engine for direct point estimate comparison
proposed_model <- PLSsemEngine:::pls_engine(simulated_data, measurement_model, structural_model)

# =====================
# 4. DATA EXTRACTION
# =====================

# 1. plspm Extraction
plspm_p1 <- plspm_model$path_coefs["Customer_Satisfaction", "Service_Quality"]
plspm_p2 <- plspm_model$path_coefs["Customer_Loyalty", "Service_Quality"]
plspm_p3 <- plspm_model$path_coefs["Customer_Loyalty", "Customer_Satisfaction"]
plspm_r2_1 <- plspm_model$inner_summary["Customer_Satisfaction", "R2"]
plspm_r2_2 <- plspm_model$inner_summary["Customer_Loyalty", "R2"]

# 2. PLSsemEngine Extraction
# Swapping p2 and p3 indices to match plspm's structural matrix order
prop_p1 <- proposed_model$paths[[1]][1]  # SQ -> CS (Remains the same)
prop_p2 <- proposed_model$paths[[2]][2]  # SQ -> CL (Changed index from 1 to 2)
prop_p3 <- proposed_model$paths[[2]][1]  # CS -> CL (Changed index from 2 to 1)

# =====================
# 5. COMPARISON REPORT
# =====================

comparison <- data.frame(
  Relationship = c(
    "Path: SQ -> CS",
    "Path: SQ -> CL",
    "Path: CS -> CL",
    "R2: Satisfaction",
    "R2: Loyalty"
  ),
  plspm = as.numeric(c(plspm_p1, plspm_p2, plspm_p3, plspm_r2_1, plspm_r2_2)),
  Proposed_Engine = as.numeric(c(prop_p1, prop_p2, prop_p3, proposed_model$r2[1], proposed_model$r2[2]))
)

# Calculate Absolute Difference
comparison$Absolute_Diff <- abs(comparison$plspm - comparison$Proposed_Engine)

# Round for professional display (Columns 2, 3, and 4)
comparison[, 2:4] <- round(comparison[, 2:4], 4)

# Print Output
cat("\n--- CROSS-SOFTWARE NUMERICAL VALIDATION ---\n")
print(comparison)

# Final Summary check
if(max(comparison$Absolute_Diff, na.rm = TRUE) < 0.005) {
  cat("\nValidation SUCCESS: Numerical equivalence confirmed.\n")
} else {
  cat("\nValidation WARNING: Significant numerical differences detected.\n")
}
