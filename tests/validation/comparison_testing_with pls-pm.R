############################################################
# VALIDATION SCRIPT – Numerical comparison with plspm
#
# This script is provided for validation purposes only.
# It is not required to run the pls-sem-engine.
############################################################

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

library(plspm)

# Path matrix
path_matrix <- rbind(
  Service_Quality = c(0, 0, 0),
  Customer_Satisfaction = c(1, 0, 0),
  Customer_Loyalty = c(1, 1, 0)
)

colnames(path_matrix) <- rownames(path_matrix)


# Blocks
blocks <- list(
  c("SQ1","SQ2","SQ3"),
  c("CS1","CS2","CS3"),
  c("CL1","CL2","CL3")
)

# Modes
modes <- c("A","A","A")

# Run plspm
plspm_model <- plspm(
  simulated_data,
  path_matrix,
  blocks,
  modes = modes
)

# Extract results
plspm_loadings <- plspm_model$outer_model[, c("name","loading")]
plspm_paths <- plspm_model$path_coefs
plspm_r2 <- plspm_model$inner_model[, c("name","R2")]

plspm_model$path_coefs

