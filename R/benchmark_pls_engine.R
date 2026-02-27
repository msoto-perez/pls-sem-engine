############################################################
# BENCHMARK – PLS ENGINE vs plspm
############################################################
# Note: pls_sem() includes additional diagnostics beyond bootstrap,
# which may increase execution time relative to plspm.

source("pls_engine.R")
library(plspm)

# =========================
# Helper: simulate data
# =========================
simulate_pls_data <- function(n) {
  
  Service_Quality <- rnorm(n)
  
  Customer_Satisfaction <-
    0.6 * Service_Quality + rnorm(n, sd = 0.6)
  
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
  
  data.frame(
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
}

# =========================
# Model specification
# =========================

measurement_model <- list(
  Service_Quality = c("SQ1","SQ2","SQ3"),
  Customer_Satisfaction = c("CS1","CS2","CS3"),
  Customer_Loyalty = c("CL1","CL2","CL3")
)

structural_model <- list(
  Customer_Satisfaction ~ Service_Quality,
  Customer_Loyalty ~ Customer_Satisfaction + Service_Quality
)

path_matrix <- rbind(
  Service_Quality = c(0, 0, 0),
  Customer_Satisfaction = c(1, 0, 0),
  Customer_Loyalty = c(1, 1, 0)
)

colnames(path_matrix) <- rownames(path_matrix)

blocks <- list(
  c("SQ1","SQ2","SQ3"),
  c("CS1","CS2","CS3"),
  c("CL1","CL2","CL3")
)

modes <- c("A","A","A")

# =========================
# Benchmark settings
# =========================

sample_sizes <- c(100, 500, 1000, 3000)
bootstrap_reps <- 200

set.seed(123)  # Ensures reproducible simulation across runs

results <- list()

for (n in sample_sizes) {
  
  cat("Running benchmark for n =", n, "\n")
  
  data_sim <- simulate_pls_data(n)
  
  # --- pls_engine (no bootstrap)
  t_engine <- system.time({
    pls_engine(
      data_sim,
      measurement_model,
      structural_model
    )
  })[3]
  
  # --- pls_sem with bootstrap
  t_pls_sem <- system.time({
    pls_sem(
      data = data_sim,
      measurement_model = measurement_model,
      structural_model = structural_model,
      nboot = bootstrap_reps,
      k = 5
    )
  })[3]
  
  # --- plspm with bootstrap
  t_plspm <- system.time({
    plspm(
      data_sim,
      path_matrix,
      blocks,
      modes = modes,
      boot.val = TRUE,
      br = bootstrap_reps
    )
  })[3]
  
  results[[length(results) + 1]] <- data.frame(
    Sample_Size = n,
    Engine_NoBoot_sec = t_engine,
    Engine_Boot200_sec = t_pls_sem,
    plspm_Boot200_sec = t_plspm
  )
}

benchmark_table <- do.call(rbind, results)

print(benchmark_table)

write.csv(benchmark_table, "benchmark_results.csv", row.names = FALSE)


