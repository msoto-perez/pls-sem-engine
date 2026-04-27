#################################################################
# MONTE CARLO BENCHMARKING (ROBUSTNESS VALIDATION)
# Standalone execution script.
# Generates synthetic data, defines the model, and runs the 
# Monte Carlo validation to assess empirical stability and bias.
#################################################################

# =================================================================
# 1. Environment Setup (Load the engine)
# =================================================================

# Install and load the latest version of the engine directly from the repository
devtools::install_github("msoto-perez/pls-sem-engine", force = TRUE)
library(plssemengine)

# =================================================================
# 2. Data Simulation and Model Definition
# =================================================================
set.seed(123)
n_obs <- 300

# Simulate Latent Variables
Service_Quality <- rnorm(n_obs)
Customer_Satisfaction <- 0.6 * Service_Quality + rnorm(n_obs, sd = 0.6)
Customer_Loyalty <- 0.55 * Customer_Satisfaction + 0.25 * Service_Quality + rnorm(n_obs, sd = 0.6)

# Helper to create items
latent_to_item <- function(latent, loading) {
  x <- loading * latent + rnorm(length(latent), sd = sqrt(1 - loading^2))
  x <- scale(x)
  as.numeric(cut(x, breaks = quantile(x, probs = seq(0, 1, length.out = 8)), 
                 labels = 1:7, include.lowest = TRUE))
}

# Create Dataset
simulated_data <- data.frame(
  SQ1 = latent_to_item(Service_Quality, 0.82), SQ2 = latent_to_item(Service_Quality, 0.78), SQ3 = latent_to_item(Service_Quality, 0.74),
  CS1 = latent_to_item(Customer_Satisfaction, 0.80), CS2 = latent_to_item(Customer_Satisfaction, 0.76), CS3 = latent_to_item(Customer_Satisfaction, 0.72),
  CL1 = latent_to_item(Customer_Loyalty, 0.81), CL2 = latent_to_item(Customer_Loyalty, 0.77), CL3 = latent_to_item(Customer_Loyalty, 0.73)
)

# Define Models
measurement_model <- list(
  Service_Quality = c("SQ1", "SQ2", "SQ3"),
  Customer_Satisfaction = c("CS1", "CS2", "CS3"),
  Customer_Loyalty = c("CL1", "CL2", "CL3")
)

structural_model <- list(
  Customer_Satisfaction ~ Service_Quality,
  Customer_Loyalty ~ Customer_Satisfaction + Service_Quality
)

# =================================================================
# 3. Monte Carlo Function Definition
# =================================================================
run_monte_carlo_benchmark <- function(data, measurement_model, structural_model, reps = 100, seed = 123) {
  set.seed(seed)
  
  # Base estimation to get "population" targets
  base_model <- plssemengine:::pls_engine(data, measurement_model, structural_model)
  target_paths <- unlist(base_model$paths)
  
  # Extract empirical correlation matrix to simulate multivariate normal data
  items <- unlist(measurement_model)
  cor_matrix <- cor(data[, items, drop = FALSE], use = "pairwise.complete.obs")
  p <- ncol(cor_matrix)
  
  # Base R Cholesky decomposition for data generation
  chol_matrix <- chol(cor_matrix)
  
  cat(sprintf("\n=================================================================\n"))
  cat(sprintf(" Running Monte Carlo Validation (%d replications)...\n", reps))
  cat(sprintf("=================================================================\n"))
  
  results_list <- list()
  
  for (i in seq_len(reps)) {
    # Generate synthetic data (Pure base R)
    Z <- matrix(rnorm(n_obs * p), nrow = n_obs, ncol = p)
    sim_data <- as.data.frame(Z %*% chol_matrix)
    colnames(sim_data) <- items
    
    # Estimate model on synthetic data
    sim_model <- tryCatch({
      plssemengine:::pls_engine(sim_data, measurement_model, structural_model)
    }, error = function(e) return(NULL))
    
    if (!is.null(sim_model)) {
      results_list[[i]] <- unlist(sim_model$paths)
    }
  }
  
  # Combine results
  res_matrix <- do.call(rbind, results_list)
  
  # Calculate empirical stability metrics
  mean_estimates <- colMeans(res_matrix, na.rm = TRUE)
  empirical_sd <- apply(res_matrix, 2, sd, na.rm = TRUE)
  bias <- mean_estimates - target_paths
  
  benchmark_report <- data.frame(
    Path = names(target_paths),
    Target_Beta = round(target_paths, 3),
    MC_Mean = round(mean_estimates, 3),
    MC_Bias = round(bias, 4),
    Empirical_SD = round(empirical_sd, 4),
    row.names = NULL
  )
  
  print(benchmark_report)
  cat("=================================================================\n\n")
  
  invisible(benchmark_report)
}

# =================================================================
# 4. Execute Benchmark
# =================================================================
# Running 100 replications. 
# (For final validation reporting, you may increase reps to 500 or 1000)
mc_results <- run_monte_carlo_benchmark(
  data = simulated_data, 
  measurement_model = measurement_model, 
  structural_model = structural_model, 
  reps = 100
)

