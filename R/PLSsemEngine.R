# Version: 1.2.0
# Date: 2026-04-25
# email: msoto@up.edu.mx

#################################################
# PLSsemEngine: Transparent PLS-SEM Mode A Estimation
#
# Design principles:
# - Methodological parsimony
# - Full transparency of estimation steps
# - No automatic interpretative rules or thresholds
# - Researcher-led model evaluation
#
# This engine is intentionally limited to reflective
# measurement models (Mode A) in its current version.
#################################################

# Core PLS-SEM estimation engine.
# Implements the standard PLS iterative algorithm for
# reflective measurement models (Mode A) using pure R.
#
# The function estimates:
# - Latent variable scores
# - Indicator loadings
# - Structural path coefficients
# - R-squared values
#
# No formative measurement, model re-specification,
# or interpretative flags are implemented by design.

#################################################################
# CORE ENGINE: PLS-SEM Estimation (Mode A)
#################################################################

#' @title Core PLS-SEM Estimation Engine (Mode A)
#'
#' @description Implements the standard Partial Least Squares (PLS) iterative algorithm 
#' restricted to reflective measurement models (Mode A). The engine estimates latent 
#' variable scores, indicator loadings, and structural path coefficients using pure 
#' matrix operations. By design, it avoids automatic model re-specification or 
#' hidden interpretative heuristics.
#'
#' @param data A data frame or matrix containing the observed indicators.
#' @param measurement_model A named list defining the reflective constructs and their corresponding indicators.
#' @param structural_model A list of formulas defining the structural relationships between latent constructs.
#' @param max_iter Integer. Maximum number of iterations for the PLS algorithm. Default is 300.
#' @param tol Numeric. Tolerance threshold for algorithmic convergence. Default is 1e-6.
#' @param sign_correction Logical. If TRUE, applies deterministic sign alignment. Disabled by default to prevent bootstrap truncation.
#' @param inner_scheme Character. Weighting scheme for the inner approximation ("factorial" or "centroid"). Default is "factorial".
#'
#' @return A list containing latent scores, outer weights, indicator loadings, structural paths, R-squared values, and metadata.
#' @export
pls_engine <- function(data,
                       measurement_model,
                       structural_model,
                       max_iter = 300,
                       tol = 1e-6,
                       sign_correction = FALSE,
                       inner_scheme = "factorial") { 
  
  # Standardise raw indicators to zero mean and unit variance for scale invariance
  X_raw <- as.matrix(data)
  X <- scale(X_raw) 
  
  if (any(is.na(X))) {
    stop("Scaling produced NA values. Check for zero-variance indicators in the dataset.")
  }
  
  blocks <- measurement_model
  constructs <- names(blocks)
  
  outer_weights <- vector("list", length(blocks))
  names(outer_weights) <- constructs
  
  # --- 1. Adjacency Matrix for Inner Model
  # Construct the connection network (C) based on the specified structural model.
  # Inner relationships in PLS are bidirectional for the purpose of weight calculation.
  C <- matrix(0, nrow = length(constructs), ncol = length(constructs), 
              dimnames = list(constructs, constructs))
  for (f in structural_model) {
    lhs <- as.character(f[[2]])
    rhs <- all.vars(f[[3]])
    for (r in rhs) {
      C[r, lhs] <- 1
      C[lhs, r] <- 1 
    }
  }
  
  # --- Initialisation of Latent Variable Scores
  # Initial proxy scores are computed as the simple, equally-weighted average 
  # of their assigned empirical indicators, followed by standardisation.
  scores <- sapply(blocks, function(items) {
    rowMeans(X[, items, drop = FALSE], na.rm = TRUE)
  })
  scores <- scale(scores)
  colnames(scores) <- constructs
  
  # --- PLS Iterative Algorithm (Mode A)
  for (iter in seq_len(max_iter)) {
    
    scores_old <- scores
    
    # --- 2. Inner Approximation
    # Constructs communicate based on the chosen inner weighting scheme.
    R <- cor(scores, use = "pairwise.complete.obs")
    if (inner_scheme == "centroid") {
      E <- sign(R) * C
    } else {
      E <- R * C # Factorial scheme as the default and recommended approach
    }
    
    # Compute the inner approximation (Z)
    Z <- scores %*% E
    
    # --- 3. Outer Approximation (Weight Update Mode A)
    for (j in seq_along(blocks)) {
      items <- blocks[[j]]
      
      # For Mode A (reflective), external weights are derived from the bivariate 
      # correlation between the indicators and their corresponding inner proxy.
      w <- apply(X[, items, drop = FALSE], 2, function(x) cor(x, Z[, j], use = "pairwise.complete.obs"))
      
      if (all(is.na(w))) next
      
      # --- 4. Fix Composite Variance to 1 
      # Normalise the newly calculated scores by dividing by their standard deviation 
      # to ensure unit variance, adhering strictly to the classic PLS algorithm.
      unscaled_score <- as.numeric(X[, items, drop = FALSE] %*% w)
      w <- w / sd(unscaled_score, na.rm = TRUE)
      
      outer_weights[[j]] <- w
      scores[, j] <- X[, items, drop = FALSE] %*% w
    }
    
    # Explicitly standardise all latent scores at the end of the iteration
    scores <- scale(scores)
    
    # Convergence criterion: algorithm stops when the maximum absolute difference 
    # between successive latent scores falls below the defined tolerance threshold.
    if (max(abs(scores - scores_old), na.rm = TRUE) < tol) break
  }
  
  # --- Sign Indeterminacy Correction
  # Disabled by default to prevent bootstrap distribution truncation (Rönkkö et al., 2015).
  # If enabled, scores are deterministically multiplied by -1 if their correlation 
  # with their dominant indicator is negative.
  if (sign_correction) {
    for (j in seq_along(blocks)) {
      items <- blocks[[j]]
      temp_loadings <- cor(X[, items, drop=FALSE], scores[, j], use = "pairwise.complete.obs")
      max_ind <- which.max(abs(temp_loadings))
      
      if (!is.na(temp_loadings[max_ind]) && temp_loadings[max_ind] < 0) {
        scores[, j] <- scores[, j] * -1
      }
    }
  }
  
  # --- Indicator Loadings
  # Calculate final reflective loadings as the correlation between items and constructs
  loadings <- matrix(NA, nrow = sum(lengths(blocks)), ncol = length(blocks))
  rownames(loadings) <- unlist(blocks)
  colnames(loadings) <- constructs
  
  for (j in seq_along(blocks)) {
    for (item in blocks[[j]]) {
      loadings[item, j] <- cor(X[, item], scores[, j], use = "pairwise.complete.obs")
    }
  }
  
  # --- Structural Model Estimation
  # Estimate path coefficients and R-squared values using Ordinary Least Squares (OLS)
  paths <- list()
  r2 <- numeric(0)
  
  for (f in structural_model) {
    lhs <- as.character(f[[2]])
    rhs <- all.vars(f[[3]])
    
    fit <- lm(scores[, lhs] ~ scores[, rhs, drop = FALSE])
    paths[[lhs]] <- coef(fit)[-1]
    r2[lhs] <- summary(fit)$r.squared
  }
  
  list(
    scores = scores,
    loadings = loadings,
    paths = paths,
    r2 = r2,
    outer_weights = outer_weights,
    inner_scheme = inner_scheme 
  )
}

#################################################################
# INFERENCE: Bootstrap Structural Model
#################################################################

#' @title Non-Parametric Bootstrap for Structural Paths
#'
#' @description Executes a case-level percentile bootstrap to assess the statistical 
#' significance of structural path coefficients. Crucially, to capture the total 
#' sampling uncertainty, this function performs a full-model re-estimation for every 
#' bootstrap replica, recalculating both outer weights and inner paths.
#'
#' @param data A data frame containing the observed indicators.
#' @param measurement_model A named list defining the reflective constructs.
#' @param structural_model A list of formulas defining the structural paths.
#' @param nboot Integer. Number of bootstrap resamples. Default is 500.
#' @param seed Integer. Random seed for reproducibility. Default is 123.
#' @param sign_correction Logical. If TRUE, applies deterministic sign alignment per resample.
#' @param inner_scheme Character. Weighting scheme for the inner approximation.
#'
#' @return A data frame containing the estimated path coefficients for all bootstrap iterations.
#' @export
bootstrap_paths <- function(data,
                            measurement_model,
                            structural_model,
                            nboot = 500,
                            seed = 123,
                            sign_correction = FALSE,
                            inner_scheme = "factorial") { 
  
  # Initialise random seed for exact reproducibility across execution runs
  set.seed(seed)
  n <- nrow(data)
  results <- list()
  
  for (b in seq_len(nboot)) {
    
    # Generate bootstrap sample via case-level resampling with replacement
    idx <- sample(seq_len(n), replace = TRUE)
    data_b <- data[idx, , drop = FALSE]
    
    # Execute full-model re-estimation to accurately capture total sampling uncertainty.
    # This prevents the underestimation of standard errors associated with partial updates.
    engine_b <- pls_engine(
      data = data_b,
      measurement_model = measurement_model,
      structural_model = structural_model,
      sign_correction = sign_correction,
      inner_scheme = inner_scheme 
    )
    
    # Extract and store path coefficients for the current resample
    for (f in structural_model) {
      lhs <- as.character(f[[2]])
      rhs <- all.vars(f[[3]])
      betas <- engine_b$paths[[lhs]]
      
      for (i in seq_along(rhs)) {
        results[[length(results) + 1]] <- data.frame(
          Path = paste(rhs[i], "->", lhs),
          Beta = betas[i]
        )
      }
    }
  }
  
  do.call(rbind, results)
}

#################################################################
# PREDICTIVE EVALUATION: PLSpredict (Cross-Validation)
#################################################################

#' @title PLSpredict: Out-of-Sample Predictive Evaluation
#'
#' @description Implements the PLSpredict algorithm using k-fold cross-validation. 
#' The function estimates the model exclusively on the training set and explicitly 
#' projects the outer weights onto the test set. Predictive performance (RMSE) is 
#' then compared against a naive linear model (LM) benchmark.
#'
#' @param data A data frame containing the observed indicators.
#' @param measurement_model A named list defining the reflective constructs.
#' @param structural_model A list of formulas defining the structural paths.
#' @param k Integer. Number of folds for cross-validation. Default is 5.
#' @param seed Integer. Random seed for fold generation to ensure reproducibility. Default is 123.
#' @param sign_correction Logical. Applies deterministic sign alignment to training models.
#' @param inner_scheme Character. Weighting scheme for the inner approximation.
#'
#' @return A list containing a predictive evaluation table (RMSE and Q2_predict) and conditional warnings.
#' @export
plspredict <- function(data,
                       measurement_model,
                       structural_model,
                       k = 5,
                       seed = 123,
                       sign_correction = FALSE,
                       inner_scheme = "factorial") { 
  
  # Set seed to guarantee that fold allocation remains consistent
  set.seed(seed)
  n <- nrow(data)
  folds <- sample(rep(1:k, length.out = n))
  
  results <- list()
  
  for (fold in seq_len(k)) {
    
    # Partition data into mutually exclusive training and test sets
    train <- data[folds != fold, , drop = FALSE]
    test  <- data[folds == fold,  , drop = FALSE]
    
    # Estimate the core PLS model exclusively on the training data 
    model <- pls_engine(train, measurement_model, structural_model, 
                        sign_correction = sign_correction,
                        inner_scheme = inner_scheme) 
    
    # Standardise test data using the exact mean and standard deviation from the training set
    X_train <- scale(as.matrix(train))
    X_test  <- scale(as.matrix(test),
                     center = attr(X_train, "scaled:center"),
                     scale  = attr(X_train, "scaled:scale"))
    
    # Allocate matrix for test scores
    test_scores <- matrix(NA,
                          nrow = nrow(X_test),
                          ncol = length(measurement_model))
    
    colnames(test_scores) <- names(measurement_model)
    
    # Project training outer weights onto the standardised test data to compute 
    # out-of-sample latent scores without causing data leakage
    for (j in seq_along(measurement_model)) {
      
      items <- measurement_model[[j]]
      w <- model$outer_weights[[j]]
      
      test_scores[, j] <- X_test[, items, drop = FALSE] %*% w
    }
    
    test_scores <- scale(test_scores)
    
    # Predict endogenous indicators structurally
    for (f in structural_model) {
      
      lhs <- as.character(f[[2]])
      rhs <- all.vars(f[[3]])
      if (length(rhs) == 0) next
      
      y_hat <- as.matrix(test_scores[, rhs, drop = FALSE]) %*%
        model$paths[[lhs]]
      
      inds <- measurement_model[[lhs]]
      
      for (ind in inds) {
        
        lambda <- model$loadings[ind, lhs]
        if (is.na(lambda)) next
        
        # Calculate PLS prediction on a standardised scale
        y_pls_std <- as.numeric(y_hat * lambda)
        
        # Re-scale predictions back to the original indicator metric using training statistics
        mean_train <- mean(train[[ind]], na.rm = TRUE)
        sd_train   <- sd(train[[ind]], na.rm = TRUE)
        
        y_pls <- y_pls_std * sd_train + mean_train
        
        y_obs <- test[[ind]]
        
        rmse_pls <- sqrt(mean((y_obs - y_pls)^2, na.rm = TRUE))
        
        # --- Linear Model (LM) Benchmark
        
        # Estimate the naive LM model using training scores
        df_train <- data.frame(
          y = train[[ind]],
          model$scores[, rhs, drop = FALSE]
        )
        lm_fit <- lm(y ~ ., data = df_train)
        
        # Predict using test scores for the LM benchmark
        df_test <- data.frame(
          test_scores[, rhs, drop = FALSE]
        )
        colnames(df_test) <- colnames(model$scores[, rhs, drop = FALSE])
        
        y_lm <- predict(lm_fit, newdata = df_test)
        
        rmse_lm <- sqrt(mean((y_obs - y_lm)^2, na.rm = TRUE))
        
        results[[length(results) + 1]] <- data.frame(
          Construct = lhs,
          Indicator = ind,
          RMSE_PLS = rmse_pls,
          RMSE_LM = rmse_lm
        )
      }
    }
  }
  
  res <- do.call(rbind, results)
  
  # --- Predictive Relevance Assessment
  # Aggregate RMSE across all folds for comparison
  table5 <- do.call(
    rbind,
    lapply(
      split(seq_len(nrow(res)),
            paste(res$Construct, res$Indicator)),
      function(i) {
        data.frame(
          Construct = res$Construct[i[1]],
          Indicator = res$Indicator[i[1]],
          RMSE_PLS = mean(res$RMSE_PLS[i]),
          RMSE_LM  = mean(res$RMSE_LM[i])
        )
      })
  )
  
  # Calculate Q2_predict: values > 0 indicate the PLS model outperforms the naive LM benchmark
  table5$Q2_predict <- ifelse(
    table5$RMSE_LM == 0,
    NA,
    1 - (table5$RMSE_PLS^2 / table5$RMSE_LM^2)
  )
  
  # Conditional diagnostic warning logic
  warnings <- character(0)
  if (any(table5$Q2_predict < 0, na.rm = TRUE)) {
    warnings <- paste(
      "Negative Q2_predict detected:",
      "PLS-based predictions are outperformed by the linear benchmark,",
      "indicating low predictive relevance (Shmueli et al., 2019)."
    )
    warning(warnings, call. = FALSE)
  }
  
  rownames(table5) <- NULL
  
  list(
    table = table5,
    warnings = warnings
  )
}
#################################################################
# DIAGNOSTICS: Common Method Bias (Full Collinearity VIF)
#################################################################

#' @title Compute Full Collinearity VIF for Common Method Bias
#'
#' @description Assesses potential common method bias using the full collinearity 
#' VIF approach proposed by Kock (2015). Values above 3.3 may indicate 
#' pathological collinearity. In line with the engine's design philosophy, 
#' no automatic classification or model modification is applied.
#'
#' @param scores A matrix or data frame of estimated latent variable scores.
#' @param digits Integer. Number of decimal places for the output. Default is 2.
#'
#' @return A data frame containing the VIF values for each construct.
#' @export
compute_cmb_vif <- function(scores, digits = 2) {
  
  constructs <- colnames(scores)
  vif <- numeric(length(constructs))
  names(vif) <- constructs
  
  # Iterate over each construct to compute its VIF against all others
  for (i in seq_along(constructs)) {
    
    y <- scores[, i]
    x <- scores[, -i, drop = FALSE]
    
    # Fit a linear model predicting the focal construct using the remaining constructs
    fit <- lm(y ~ x)
    r2 <- summary(fit)$r.squared
    
    # Calculate Variance Inflation Factor (VIF)
    vif[i] <- 1 / (1 - r2)
  }
  
  data.frame(
    Construct = constructs,
    VIF = round(vif, digits),
    row.names = NULL
  )
}

#################################################################
# DIAGNOSTICS: Discriminant Validity (HTMT & HTMT2)
#################################################################

#' @title Compute HTMT and HTMT2 Ratios for Discriminant Validity
#'
#' @description Computes the Heterotrait-Monotrait ratio of correlations (HTMT) 
#' following Henseler et al. (2015), and the revised HTMT2 following Roemer 
#' et al. (2021). HTMT2 uses the geometric mean to relax the tau-equivalence 
#' assumption, making it suitable for congeneric models.
#'
#' @param data A data frame containing the observed indicators.
#' @param measurement_model A named list defining the reflective constructs.
#' @param digits Integer. Number of decimal places for the output matrices. Default is 2.
#'
#' @return A list containing the HTMT and HTMT2 symmetric matrices.
#' @export
compute_htmt_metrics <- function(data, measurement_model, digits = 2) {
  
  # Standardise data to ensure scale invariance in correlations
  X <- scale(as.matrix(data))
  constructs <- names(measurement_model)
  
  # Initialise empty matrices for both metrics
  htmt_matrix <- matrix(NA, nrow = length(constructs), ncol = length(constructs), dimnames = list(constructs, constructs))
  htmt2_matrix <- matrix(NA, nrow = length(constructs), ncol = length(constructs), dimnames = list(constructs, constructs))
  
  # Helper function for geometric mean (strictly requires positive values)
  geom_mean <- function(x) {
    x_pos <- x[x > 0 & !is.na(x)]
    if (length(x_pos) == 0) return(NA)
    exp(mean(log(x_pos)))
  }
  
  # Iterate over unique pairs of constructs
  for (i in seq_along(constructs)) {
    for (j in seq_along(constructs)) {
      
      if (i >= j) next # Calculate only the lower triangle to save computation time
      
      items_i <- measurement_model[[constructs[i]]]
      items_j <- measurement_model[[constructs[j]]]
      
      # Extract absolute heterotrait correlations (between constructs)
      cor_ij <- abs(cor(X[, items_i, drop = FALSE], X[, items_j, drop = FALSE], use = "pairwise.complete.obs"))
      
      # Extract absolute monotrait correlations (within constructs)
      cor_ii <- abs(cor(X[, items_i, drop = FALSE], use = "pairwise.complete.obs"))
      cor_jj <- abs(cor(X[, items_j, drop = FALSE], use = "pairwise.complete.obs"))
      
      cor_ii_vec <- cor_ii[lower.tri(cor_ii)]
      cor_jj_vec <- cor_jj[lower.tri(cor_jj)]
      
      # --- HTMT (Arithmetic mean approach)
      mean_hetero_a <- mean(cor_ij, na.rm = TRUE)
      mean_mono_i_a <- mean(cor_ii_vec, na.rm = TRUE)
      mean_mono_j_a <- mean(cor_jj_vec, na.rm = TRUE)
      
      htmt_val <- mean_hetero_a / sqrt(mean_mono_i_a * mean_mono_j_a)
      htmt_matrix[i, j] <- htmt_val
      htmt_matrix[j, i] <- htmt_val
      
      # --- HTMT2 (Geometric mean approach for congeneric models)
      mean_hetero_g <- geom_mean(as.vector(cor_ij))
      mean_mono_i_g <- geom_mean(cor_ii_vec)
      mean_mono_j_g <- geom_mean(cor_jj_vec)
      
      htmt2_val <- mean_hetero_g / sqrt(mean_mono_i_g * mean_mono_j_g)
      htmt2_matrix[i, j] <- htmt2_val
      htmt2_matrix[j, i] <- htmt2_val
    }
  }
  
  list(
    HTMT = round(htmt_matrix, digits),
    HTMT2 = round(htmt2_matrix, digits)
  )
}

#################################################################
# DIAGNOSTICS: Global Model Fit (SRMR, d_ULS, d_G)
#################################################################

#' @title Compute Global Model Fit Indices
#'
#' @description Calculates global fit indices including the Standardised Root Mean 
#' Square Residual (SRMR), exact squared Euclidean distance (d_ULS), and exact 
#' geodesic distance (d_G). Calculations are based on the saturated model-implied 
#' correlation matrix following Henseler et al. (2014).
#'
#' @param engine An object containing the core PLS estimation results.
#' @param data A data frame containing the observed indicators.
#' @param measurement_model A named list defining the reflective constructs.
#' @param digits Integer. Number of decimal places for the output. Default is 3.
#'
#' @return A data frame containing the computed fit metrics.
#' @export
compute_model_fit <- function(engine, data, measurement_model, digits = 3) {
  
  items <- unlist(measurement_model)
  
  # Empirical correlation matrix (S)
  S <- cor(data[, items, drop = FALSE], use = "pairwise.complete.obs")
  
  # Build strict loading matrix (Lambda) assuming 0 for cross-loadings
  constructs <- names(measurement_model)
  Lambda <- matrix(0, nrow = length(items), ncol = length(constructs),
                   dimnames = list(items, constructs))
  
  for (cn in constructs) {
    inds <- measurement_model[[cn]]
    Lambda[inds, cn] <- engine$loadings[inds, cn]
  }
  
  # Construct correlation matrix representing the saturated inner model
  R <- cor(engine$scores, use = "pairwise.complete.obs")
  
  # Model-implied correlation matrix (Sigma_hat)
  Sigma_hat <- Lambda %*% R %*% t(Lambda)
  diag(Sigma_hat) <- 1
  
  # Extract lower triangular residuals (e)
  lower_tri_idx <- lower.tri(S)
  e <- S[lower_tri_idx] - Sigma_hat[lower_tri_idx]
  
  # 1. SRMR (Standardised Root Mean Square Residual)
  srmr <- sqrt(mean(e^2))
  
  # 2. d_ULS (Squared Euclidean distance)
  d_uls <- sum(e^2)
  
  # 3. d_G (Geodesic distance)
  # Calculated as half the sum of squared logs of the eigenvalues of S * Sigma_hat^-1
  d_g <- NA
  tryCatch({
    eig_vals <- eigen(S %*% solve(Sigma_hat), only.values = TRUE)$values
    eig_vals <- Re(eig_vals)
    eig_vals <- eig_vals[eig_vals > 0] # Retain only positive real parts for stability
    d_g <- 0.5 * sum((log(eig_vals))^2)
  }, error = function(e) {
    d_g <<- NA # Handle potential singularity issues gracefully
  })
  
  data.frame(
    Metric = c("SRMR", "d_ULS", "d_G"),
    Value = c(round(srmr, digits), round(d_uls, digits), round(d_g, digits)),
    stringsAsFactors = FALSE
  )
}

#################################################################
# WRAPPER: High-Level Analysis Workflow (Version 1.2.0)
#################################################################

#' @title Estimate a PLS-SEM Model (Mode A)
#'
#' @description This is the main high-level wrapper of the PLSsemEngine. It coordinates 
#' the estimation, assessment, and predictive evaluation stages, returning all results 
#' in a structured, publication-ready S3 object. It deliberately restricts estimation 
#' to reflective blocks (Mode A) and avoids hidden model modifications.
#'
#' @param data A data frame containing the observed indicators.
#' @param measurement_model A named list defining the reflective constructs and their items.
#' @param structural_model A list of formulas defining the inner paths.
#' @param k Integer. Number of folds for PLSpredict cross-validation. Default is 5.
#' @param nboot Integer. Number of bootstrap resamples for inference. Default is 500.
#' @param digits Integer. Number of decimal places for output tables. Default is 2.
#' @param sign_correction Logical. If TRUE, applies deterministic sign alignment.
#' @param inner_scheme Character. Weighting scheme for the inner model ("factorial" or "centroid"). Default is "factorial".
#'
#' @return A list of class 'pls_model' containing structured tables for measurement, structural, and predictive assessments.
#' @importFrom stats coef cor lm predict quantile sd
#' @importFrom graphics arrows rect strheight strwidth text
#' @importFrom grDevices dev.off png
#' @importFrom utils write.csv
#' @export
pls_sem <- function(data,
                    measurement_model,
                    structural_model,
                    k = 5,
                    nboot = 500,
                    digits = 2,
                    sign_correction = FALSE,
                    inner_scheme = "factorial") { 
  
  # =============================================================
  # 1. CORE ESTIMATION
  # =============================================================
  # Execute the iterative PLS algorithm (Mode A) to obtain latent scores,
  # outer weights, and base structural coefficients.
  engine <- pls_engine(
    data = data,
    measurement_model = measurement_model,
    structural_model = structural_model,
    sign_correction = sign_correction,
    inner_scheme = inner_scheme 
  )
  
  # =============================================================
  # 2. MEASUREMENT MODEL ASSESSMENT
  # =============================================================
  # Extract item loadings, assigning each indicator to its dominant construct
  t1_items <- data.frame(
    Construct = apply(engine$loadings, 1, function(x)
      colnames(engine$loadings)[which.max(abs(x))]),
    Item = rownames(engine$loadings),
    Loading = round(
      apply(engine$loadings, 1, function(x) max(x, na.rm = TRUE)),
      digits
    ),
    row.names = NULL
  )
  
  # Calculate Construct Reliability (CR), Average Variance Extracted (AVE), 
  # and R-squared for endogenous constructs
  t1_constructs <- do.call(
    rbind,
    lapply(names(measurement_model), function(cn) {
      items <- measurement_model[[cn]]
      lambda <- engine$loadings[items, cn]
      lambda2 <- lambda^2
      
      # Composite Reliability (CR) calculation based on standardised loadings
      den <- (sum(lambda))^2 + sum(1 - lambda2)
      CR <- ifelse(den == 0, NA, (sum(lambda))^2 / den)
      
      # Average Variance Extracted (AVE) calculation
      AVE <- mean(lambda2)
      
      # Retrieve R-squared if the construct acts as an endogenous variable
      r2_val <- ifelse(cn %in% names(engine$r2), engine$r2[cn], NA)
      
      data.frame(
        Construct = cn,
        `Composite Reliability (CR)` = round(CR, digits),
        AVE = round(AVE, digits),
        R2 = ifelse(is.na(r2_val), NA, round(r2_val, digits)),
        row.names = NULL,
        check.names = FALSE
      )
    })
  )
  
  # Merge item and construct metrics into a consolidated, clean reporting table
  measurement_results <- merge(t1_items, t1_constructs, by = "Construct", all.x = TRUE)
  measurement_results <- measurement_results[order(measurement_results$Construct, measurement_results$Item), ]
  rownames(measurement_results) <- NULL
  
  # =============================================================
  # 3. DISCRIMINANT VALIDITY 
  # =============================================================
  # Compute HTMT (classic) and HTMT2 (congeneric) metrics
  htmt_results <- compute_htmt_metrics(data, measurement_model, digits)
  
  # =============================================================
  # 4. STRUCTURAL MODEL & INFERENCE 
  # =============================================================
  # Execute full-model bootstrap resampling to derive standard errors and CIs
  boot <- bootstrap_paths(
    data = data,
    measurement_model = measurement_model,
    structural_model = structural_model,
    nboot = nboot,
    sign_correction = sign_correction,
    inner_scheme = inner_scheme 
  )
  
  # Compile structural path coefficients, 95% confidence intervals, and f2 effect sizes
  structural_results <- do.call(
    rbind,
    lapply(structural_model, function(f) {
      
      lhs <- as.character(f[[2]])
      rhs <- all.vars(f[[3]])
      
      betas <- engine$paths[[lhs]]
      r2_incl <- engine$r2[lhs] # R-squared with the focal predictor included
      
      do.call(
        rbind,
        lapply(seq_along(rhs), function(i) {
          
          path_name <- paste(rhs[i], "->", lhs)
          bdist <- boot$Beta[boot$Path == path_name]
          
          # f2 Effect Size calculation (incremental predictive contribution)
          rhs_excl <- setdiff(rhs, rhs[i])
          if (length(rhs_excl) == 0) {
            r2_excl <- 0
          } else {
            fit_excl <- lm(engine$scores[, lhs] ~ engine$scores[, rhs_excl, drop = FALSE])
            r2_excl <- summary(fit_excl)$r.squared
          }
          
          # Mathematical safeguard against R2 values >= 1
          f2_val <- ifelse(r2_incl >= 1, NA, (r2_incl - r2_excl) / (1 - r2_incl))
          
          data.frame(
            From = rhs[i],
            To = lhs,
            `Path Coefficient (beta)` = round(betas[i], digits),
            CI_low = round(quantile(bdist, 0.025), digits),
            CI_high = round(quantile(bdist, 0.975), digits),
            f2 = round(f2_val, digits),
            row.names = NULL,
            check.names = FALSE
          )
        })
      )
    })
  )
  
  # =============================================================
  # 5. PREDICTIVE EVALUATION 
  # =============================================================
  # Execute k-fold cross-validation to assess out-of-sample predictive relevance
  pred <- plspredict(
    data = data,
    measurement_model = measurement_model,
    structural_model = structural_model,
    k = k,
    sign_correction = sign_correction,
    inner_scheme = inner_scheme 
  )
  
  # =============================================================
  # 6. COMPLEMENTARY DIAGNOSTICS 
  # =============================================================
  # Full Collinearity VIF for detecting potential Common Method Bias
  cmb_table <- compute_cmb_vif(engine$scores)
  
  # Global Model Fit indices (SRMR, d_ULS, d_G) based on saturated inner matrices
  fit_table <- compute_model_fit(engine, data, measurement_model, digits)
  
  # =============================================================
  # 7. FINAL OBJECT ASSEMBLY
  # =============================================================
  # Organise results into a structured S3 object to satisfy usability and traceability
  model <- list(
    measurement_model = measurement_results,     
    discriminant_validity = htmt_results,        
    structural_model = structural_results,       
    predictive_relevance = pred$table,           
    diagnostics = list(
      common_method_bias = cmb_table,            
      global_fit = fit_table                     
    ),
    latent_scores = as.data.frame(round(engine$scores, digits)),
    raw_engine = engine,                         # Retained for advanced internal diagnostics
    specification = list(
      measurement = measurement_model, 
      structural = structural_model
    ),
    warnings = list(
      measurement = NULL,
      structural = NULL,
      predictive = pred$warnings
    ),
    meta = list(
      n_obs = nrow(data),
      constructs = names(measurement_model),
      inner_scheme = inner_scheme,
      version = "1.2.0"
    )
  )
  
  class(model) <- "pls_model"
  return(model)
}

#################################################################
# UTILITIES: Export Latent Scores
#################################################################

#' @title Export Latent Variable Scores
#'
#' @description Extracts and exports the estimated latent variable scores to a 
#' CSV file. A sequential 'Case' identifier is appended to improve data 
#' traceability and facilitate subsequent external analyses.
#'
#' @param model An object of class 'pls_model' generated by \code{pls_sem}.
#' @param file Character. The file path and name for the exported CSV. Default is "latent_scores.csv".
#'
#' @return Invisibly returns the data frame of scores while writing the CSV file to disk.
#' @export
export_scores <- function(model, file = "latent_scores.csv") {
  
  # Access the internally stored latent scores from the model object
  scores <- as.data.frame(model$latent_scores)
  
  # Append a sequential case identifier to guarantee traceability 
  scores$Case <- seq_len(nrow(scores))
  
  # Reorder columns to ensure the 'Case' identifier appears first
  scores <- scores[, c("Case", setdiff(colnames(scores), "Case"))]
  
  write.csv(scores, file, row.names = FALSE)
  invisible(scores)
}

#################################################################
# DIAGNOSTICS: HTMT Item-Level Diagnostics
#################################################################

#' @title HTMT-Guided Item Diagnostics
#'
#' @description When discriminant validity is compromised at the construct level 
#' (HTMT > threshold), this function provides a granular, item-level diagnostic. 
#' It calculates cross-construct correlations to help researchers isolate specific 
#' problematic indicators that may be inflating the trait ratios.
#'
#' @param model An object of class 'pls_model'.
#' @param threshold Numeric. The HTMT threshold used to flag discriminant validity issues. Default is 0.85.
#' @param digits Integer. Number of decimal places for the output. Default is 2.
#'
#' @return A data frame detailing problematic construct pairs and their item-level cross-correlations, or NULL if no issues are detected.
#' @export
htmt_item_diagnostics <- function(model, threshold = 0.85, digits = 2) {
  
  # Extract the standard arithmetic HTMT matrix from the model object
  htmt <- model$discriminant_validity$HTMT 
  
  scores <- model$latent_scores
  loadings <- model$measurement_model
  
  constructs <- colnames(htmt)
  
  # Identify matrix indices where the HTMT exceeds the conservative threshold
  problems <- which(htmt > threshold, arr.ind = TRUE)
  
  if (nrow(problems) == 0) return(NULL)
  
  out <- list()
  
  for (k in seq_len(nrow(problems))) {
    
    c1 <- constructs[problems[k, 1]]
    c2 <- constructs[problems[k, 2]]
    
    items_c1 <- loadings$Item[loadings$Construct == c1]
    items_c2 <- loadings$Item[loadings$Construct == c2]
    
    # Execute item-construct cross-correlation analysis for the first construct
    for (it in items_c1) {
      # Utilise pairwise.complete.obs to ensure statistical robustness against missing data
      r <- cor(scores[, c2], scores[, c1], use = "pairwise.complete.obs")
      out[[length(out) + 1]] <- data.frame(
        Construct_A = c1,
        Construct_B = c2,
        Item = it,
        `Construct-level correlation` = round(abs(r), digits)
      )
    }
    
    # Execute item-construct cross-correlation analysis for the second construct
    for (it in items_c2) {
      r <- cor(scores[, c1], scores[, c2], use = "pairwise.complete.obs")
      out[[length(out) + 1]] <- data.frame(
        Construct_A = c2,
        Construct_B = c1,
        Item = it,
        `Construct-level correlation` = round(abs(r), digits)
      )
    }
  }
  
  do.call(rbind, out)
}

#################################################################
# INFERENCE: Indirect Effects (Mediation)
#################################################################

#' @title Calculate Indirect Effects
#'
#' @description Computes indirect effects as the mathematical product of the 
#' respective structural path coefficients. Consistent with the engine's design 
#' philosophy, no automated mediation classifications (e.g., "full" or "partial") 
#' are provided, ensuring the theoretical interpretation remains entirely in the 
#' hands of the researcher.
#'
#' @param model An object of class 'pls_model'.
#' @param digits Integer. Number of decimal places for the output. Default is 3.
#'
#' @return A data frame containing the calculated indirect effects, or NULL if no structural mediators exist.
#' @export
get_indirect_effects <- function(model, digits = 3) {
  
  t4 <- model$structural_model
  
  if (is.null(t4)) return(NULL)
  
  # Robust extraction of the beta coefficient column to ensure version compatibility
  beta_col <- grep("Path", colnames(t4), value = TRUE)[1]
  
  out <- list()
  
  # Identify genuine mediators (constructs acting as both antecedent and outcome)
  mediators <- intersect(unique(t4$To), unique(t4$From))
  
  for (m in mediators) {
    
    incoming <- t4[t4$To == m, , drop = FALSE]
    outgoing <- t4[t4$From == m, , drop = FALSE]
    
    if (nrow(incoming) == 0 || nrow(outgoing) == 0) next
    
    # Calculate the product of coefficients for all incoming and outgoing paths
    for (i in seq_len(nrow(incoming))) {
      for (j in seq_len(nrow(outgoing))) {
        
        out[[length(out) + 1]] <- data.frame(
          From = incoming$From[i],
          Via  = m,
          To   = outgoing$To[j],
          `Indirect Effect` = round(
            incoming[[beta_col]][i] * outgoing[[beta_col]][j],
            digits
          ),
          stringsAsFactors = FALSE
        )
      }
    }
  }
  
  if (length(out) == 0) return(NULL)
  do.call(rbind, out)
}

#################################################################
# UTILITIES: Methodological References
#################################################################

#' @title Retrieve Standard Methodological References
#'
#' @description Provides a quick-reference table outlining standard evaluation 
#' heuristics and literature citations for PLS-SEM. This tool is intended to 
#' support researcher-led interpretation rather than enforce mechanical rule application.
#'
#' @return A data frame mapping analytical components to standard thresholds and academic references.
#' @export
get_references <- function() {
  
  data.frame(
    Component = c(
      "Indicator loadings",
      "Composite Reliability (CR)",
      "Average Variance Extracted (AVE)",
      "HTMT",
      "Structural paths",
      "Indirect effects",
      "PLSpredict",
      "Common Method Bias"
    ),
    Threshold = c(
      ">= 0.70",
      ">= 0.70",
      ">= 0.50",
      "< 0.85",
      "Bootstrap CI excludes 0",
      "Product of path coefficients",
      "RMSE vs LM benchmark",
      "VIF < 3.3"
    ),
    Reference = c(
      "Hair et al. (2017, 2022)",
      "Hair et al. (2017)",
      "Fornell & Larcker (1981)",
      "Henseler et al. (2015)",
      "Hair et al. (2017)",
      "Hair et al. (2017)",
      "Shmueli et al. (2019)",
      "Kock (2015)"
    ),
    stringsAsFactors = FALSE
  )
}

#################################################################
# VISUALISATION: Helper for Boundary Intersection
#################################################################

#' @title Calculate Exact Boundary Intersections for Rectangles
#'
#' @description A geometric helper function that calculates the precise offset required 
#' for arrow positioning, ensuring that structural paths touch the boundaries of the 
#' construct rectangles rather than overlapping them.
#'
#' @param angle Numeric. The angle of the path in radians.
#' @param bw Numeric. The half-width of the rectangle.
#' @param bh Numeric. The half-height of the rectangle.
#'
#' @return Numeric. The geometric offset distance.
#' @export
get_boundary_offset <- function(angle, bw, bh) {
  # Apply trigonometric logic to determine if the intersection occurs 
  # on the vertical or horizontal boundaries of the node
  if (abs(tan(angle)) < bh/bw) {
    return(bw / abs(cos(angle)))
  } else {
    return(bh / abs(sin(angle)))
  }
}

#################################################################
# VISUALISATION: Plot Structural Model (Base R)
#################################################################

#' @title Plot Structural Model Base Geometry
#'
#' @description Generates a pure base R plot of the structural model layout. 
#' It relies exclusively on native graphical parameters to avoid third-party 
#' dependencies, ensuring long-term reproducibility and stability.
#'
#' @param structural_model A list of formulas defining the paths.
#' @param layout Optional data frame with 'name', 'x', and 'y' coordinates for nodes.
#' @param box_width Numeric. Half-width of the construct nodes. Default is 1.1.
#' @param box_height Numeric. Half-height of the construct nodes. Default is 0.3.
#' @param cex_node Numeric. Text expansion factor for node labels. Default is 0.9.
#' @param arr_lwd Numeric. Line width for structural arrows. Default is 1.2.
#' @param box_col Character. Background colour for nodes. Default is "white".
#' @param save_plot Logical. If TRUE, exports the plot as a high-resolution PNG.
#' @param file_name Character. Target file name for export.
#' @param width Integer. Width of the exported image in pixels.
#' @param height Integer. Height of the exported image in pixels.
#' @param res Integer. Resolution of the exported image in PPI.
#'
#' @return Invisibly generates a plot on the active graphic device.
#' @export
plot_structural_model <- function(structural_model, layout = NULL, 
                                  box_width = 1.1, box_height = 0.3, 
                                  cex_node = 0.9, arr_lwd = 1.2, 
                                  box_col = "white",
                                  save_plot = FALSE, file_name = "structural_model.png",
                                  width = 2500, height = 1500, res = 300) {
  
  if (save_plot) {
    if (capabilities("cairo")) {
      png(filename = file_name, width = width, height = height, res = res, type = "cairo")
    } else {
      png(filename = file_name, width = width, height = height, res = res)
    }
  }
  
  # Extract node and edge lists from the structural model specification
  paths <- lapply(structural_model, function(eq) {
    data.frame(from = all.vars(eq)[-1], to = all.vars(eq)[1], stringsAsFactors = FALSE)
  })
  edges <- do.call(rbind, paths)
  nodes <- unique(c(edges$from, edges$to))
  
  # Compute default circular layout if no explicit coordinates are provided
  if (is.null(layout)) {
    n_nodes <- length(nodes)
    angles <- seq(pi/2, 2*pi + pi/2, length.out = n_nodes + 1)[-(n_nodes + 1)]
    radius <- max(3.0, n_nodes * 0.7) 
    node_pos <- data.frame(name = nodes, x = radius * cos(angles), y = radius * sin(angles), stringsAsFactors = FALSE)
  } else { node_pos <- layout }
  
  # Initialise the blank canvas
  plot(0, 0, type = "n", 
       xlim = c(min(node_pos$x) - box_width*1.2, max(node_pos$x) + box_width*1.2), 
       ylim = c(min(node_pos$y) - box_height*1.5, max(node_pos$y) + box_height*1.5), 
       axes = FALSE, xlab = "", ylab = "", main = "")
  
  # Draw structural paths with exact boundary intersections
  for(i in 1:nrow(edges)) {
    pos_from <- node_pos[node_pos$name == edges$from[i], ]
    pos_to <- node_pos[node_pos$name == edges$to[i], ]
    dx <- pos_to$x - pos_from$x; dy <- pos_to$y - pos_from$y
    angle <- atan2(dy, dx)
    
    d_from <- get_boundary_offset(angle, box_width, box_height) + 0.05
    d_to <- get_boundary_offset(angle + pi, box_width, box_height) + 0.05
    
    arrows(x0 = pos_from$x + d_from * cos(angle), y0 = pos_from$y + d_from * sin(angle), 
           x1 = pos_to$x - d_to * cos(angle), y1 = pos_to$y - d_to * sin(angle), 
           length = 0.12, lwd = arr_lwd, col = "gray40")
  }
  
  # Draw construct nodes
  for(i in 1:nrow(node_pos)) {
    rect(node_pos$x[i]-box_width, node_pos$y[i]-box_height, 
         node_pos$x[i]+box_width, node_pos$y[i]+box_height, 
         col=box_col, border="black", lwd=1.5)
    text(node_pos$x[i], node_pos$y[i], labels = gsub("_", " ", node_pos$name[i]), cex = cex_node, font = 2, col="black")
  }
  
  if (save_plot) {
    dev.off()
    message(paste("High-quality plot saved to:", file.path(getwd(), file_name)))
  }
}

#################################################################
# VISUALISATION: Plot Structural Model with Results
#################################################################

#' @title Plot Model Results (Beta & R-squared)
#'
#' @description Overlays the estimated path coefficients and R-squared values 
#' onto the structural model plot. In response to usability requests, R-squared 
#' values can be toggled on or off to provide maximum formatting flexibility.
#'
#' @param model An object of class 'pls_model'.
#' @param layout Optional data frame with node coordinates.
#' @param show_r2 Logical. If TRUE, displays R-squared values for endogenous constructs. Default is TRUE.
#' @param box_width Numeric. Half-width of the nodes.
#' @param box_height Numeric. Half-height of the nodes.
#' @param cex_node Numeric. Text expansion for node labels.
#' @param cex_beta Numeric. Text expansion for path coefficients.
#' @param cex_r2 Numeric. Text expansion for R-squared values.
#' @param arr_lwd Numeric. Line width for arrows.
#' @param box_col Character. Background colour for nodes.
#' @param save_plot Logical. Exports plot if TRUE.
#' @param file_name Character. File name for export.
#' @param width Integer. Image width.
#' @param height Integer. Image height.
#' @param res Integer. Image resolution.
#'
#' @return Invisibly generates a plot on the active graphic device.
#' @export
plot_model_results <- function(model, layout = NULL, 
                               show_r2 = TRUE,
                               box_width = 1.1, box_height = 0.3, 
                               cex_node = 0.9, cex_beta = 0.8, cex_r2 = 0.75,
                               arr_lwd = 1.2, box_col = "white",
                               save_plot = FALSE, file_name = "model_results.png",
                               width = 2500, height = 1500, res = 300) {
  
  if (save_plot) {
    if (capabilities("cairo")) {
      png(filename = file_name, width = width, height = height, res = res, type = "cairo")
    } else {
      png(filename = file_name, width = width, height = height, res = res)
    }
  }
  
  structural_model <- model$specification$structural
  t1 <- model$measurement_model 
  t4 <- model$structural_model  
  
  paths <- lapply(structural_model, function(eq) {
    data.frame(from = all.vars(eq)[-1], to = all.vars(eq)[1], stringsAsFactors = FALSE)
  })
  edges <- do.call(rbind, paths)
  nodes <- unique(c(edges$from, edges$to))
  
  if (is.null(layout)) {
    n_nodes <- length(nodes)
    angles <- seq(pi/2, 2*pi + pi/2, length.out = n_nodes + 1)[-(n_nodes + 1)]
    radius <- max(3.0, n_nodes * 0.7) 
    node_pos <- data.frame(name = nodes, x = radius * cos(angles), y = radius * sin(angles), stringsAsFactors = FALSE)
  } else { node_pos <- layout }
  
  plot(0, 0, type = "n", 
       xlim = c(min(node_pos$x) - box_width*1.2, max(node_pos$x) + box_width*1.2), 
       ylim = c(min(node_pos$y) - box_height*1.5, max(node_pos$y) + box_height*1.5), 
       axes = FALSE, xlab = "", ylab = "", main = "")
  
  # Drawing Arrows and Beta Coefficients
  for(i in 1:nrow(edges)) {
    pos_from <- node_pos[node_pos$name == edges$from[i], ]
    pos_to <- node_pos[node_pos$name == edges$to[i], ]
    dx <- pos_to$x - pos_from$x; dy <- pos_to$y - pos_from$y
    angle <- atan2(dy, dx)
    
    d_from <- get_boundary_offset(angle, box_width, box_height) + 0.05
    d_to <- get_boundary_offset(angle + pi, box_width, box_height) + 0.05
    
    x0_arr <- pos_from$x + d_from * cos(angle); y0_arr <- pos_from$y + d_from * sin(angle)
    x1_arr <- pos_to$x - d_to * cos(angle); y1_arr <- pos_to$y - d_to * sin(angle)
    
    arrows(x0 = x0_arr, y0 = y0_arr, x1 = x1_arr, y1 = y1_arr, length = 0.12, lwd = arr_lwd, col = "gray40")
    
    # Robust identification of the coefficient column across potential language/version changes
    beta_col <- grep("Path", colnames(t4), value = TRUE)[1]
    beta_val <- t4[[beta_col]][t4$From == edges$from[i] & t4$To == edges$to[i]]
    
    mid_x <- (x0_arr + x1_arr)/2; mid_y <- (y0_arr + y1_arr)/2
    lbl <- paste0("beta=", beta_val)
    
    # Plot coefficient label with a background mask for readability
    w <- strwidth(lbl, cex=cex_beta) * 0.75 
    h <- strheight(lbl, cex=cex_beta) * 0.9  
    rect(mid_x - w, mid_y - h, mid_x + w, mid_y + h, col="white", border=NA)
    text(mid_x, mid_y, labels = lbl, cex = cex_beta, font = 3, col = "black")
  }
  
  # Drawing Rectangles and optional R2
  for(i in 1:nrow(node_pos)) {
    rect(node_pos$x[i]-box_width, node_pos$y[i]-box_height, 
         node_pos$x[i]+box_width, node_pos$y[i]+box_height, 
         col=box_col, border="black", lwd=1.5)
    
    clean_name <- gsub("_", " ", node_pos$name[i])
    r2_val <- t1$R2[t1$Construct == node_pos$name[i]][1]
    
    if(show_r2 && !is.na(r2_val)) {
      text(node_pos$x[i], node_pos$y[i] + (box_height * 0.25), labels = clean_name, cex=cex_node, font=2, col="black")
      text(node_pos$x[i], node_pos$y[i] - (box_height * 0.35), labels = paste0("(R2=", r2_val, ")"), cex=cex_r2, font=3, col="gray30")
    } else {
      text(node_pos$x[i], node_pos$y[i], labels = clean_name, cex=cex_node, font=2, col="black")
    }
  }
  
  if (save_plot) {
    dev.off()
    message(paste("High-quality plot saved to:", file.path(getwd(), file_name)))
  }
}

#################################################################
# METHODOLOGICAL ASSESSMENT: Optional Interpretive Layer
#################################################################

#' @title Optional Methodological Assessment
#'
#' @description Provides a narrative diagnostic layer based on established PLS-SEM 
#' literature (e.g., Hair et al., 2017; Henseler et al., 2015). Crucially, this 
#' function contextualises the raw metrics without imposing mechanical or automated 
#' decision-making, keeping the final analytical judgement with the researcher.
#'
#' @param model An object of class 'pls_model'.
#'
#' @return Invisibly prints a structured console report flagging potential issues.
#' @export
interpret_model <- function(model) {
  
  if (!inherits(model, "pls_model")) {
    stop("Object must be of class 'pls_model'")
  }
  
  cat("\n=================================================================\n")
  cat(" OPTIONAL METHODOLOGICAL ASSESSMENT (HEURISTIC AIDS)\n")
  cat(" Note: The following classifications are based on standard \n")
  cat(" literature thresholds. They are provided to contextualise \n")
  cat(" results, not to replace researcher judgement.\n")
  cat("=================================================================\n\n")
  
  # --- 1. Measurement Model (Reliability & AVE)
  cat("--- 1. Reflective Measurement Model (Hair et al., 2017) ---\n")
  
  t1 <- model$measurement_model 
  
  cr_issues <- t1$Construct[!is.na(t1$`Composite Reliability (CR)`) & t1$`Composite Reliability (CR)` < 0.70]
  ave_issues <- t1$Construct[!is.na(t1$AVE) & t1$AVE < 0.50]
  
  if(length(cr_issues) == 0) {
    cat(" [*] Composite Reliability: All constructs meet the >= 0.70 guideline.\n")
  } else {
    cat(" [!] Composite Reliability: Constructs below 0.70 ->", paste(unique(cr_issues), collapse=", "), "\n")
  }
  
  if(length(ave_issues) == 0) {
    cat(" [*] Convergent Validity (AVE): All constructs meet the >= 0.50 guideline.\n")
  } else {
    cat(" [!] Convergent Validity (AVE): Constructs below 0.50 ->", paste(unique(ave_issues), collapse=", "), "\n")
  }
  
  # --- 2. Discriminant Validity (HTMT & HTMT2)
  cat("\n--- 2. Discriminant Validity (Henseler et al., 2015; Roemer et al., 2021) ---\n")
  
  htmt <- model$discriminant_validity$HTMT
  htmt2 <- model$discriminant_validity$HTMT2
  
  check_htmt_vals <- function(mat, threshold = 0.85) {
    if (is.null(mat)) return("Not calculated")
    issues <- which(mat > threshold, arr.ind = TRUE)
    if(nrow(issues) == 0) return("None")
    unique_pairs <- unique(apply(issues, 1, function(x) {
      names <- sort(c(rownames(mat)[x[1]], colnames(mat)[x[2]]))
      paste(names[1], "-", names[2])
    }))
    paste(unique(unique_pairs), collapse = ", ")
  }
  
  cat(" [*] HTMT pairs > 0.85 (Conservative):", check_htmt_vals(htmt, 0.85), "\n")
  cat(" [*] HTMT2 pairs > 0.85 (Congeneric):", check_htmt_vals(htmt2, 0.85), "\n")
  
  # --- 3. Collinearity (Common Method Bias)
  cat("\n--- 3. Full Collinearity VIF (Kock, 2015) ---\n")
  
  cmb <- model$diagnostics$common_method_bias
  
  cmb_issues <- cmb$Construct[!is.na(cmb$VIF) & cmb$VIF > 3.3]
  
  if(length(cmb_issues) == 0) {
    cat(" [*] CMB VIF: All constructs meet the <= 3.3 guideline.\n") 
  } else {
    cat(" [!] CMB VIF: Constructs > 3.3 guideline ->", paste(cmb_issues, collapse=", "), "\n")
  }
  
  cat("\n=================================================================\n")
}

#################################################################
# METHODOLOGICAL BRIDGE: Lavaan Syntax Integration
#################################################################

#' @title Export to Lavaan Syntax (CB-SEM / CFA Integration)
#'
#' @description Actively encourages cross-methodological validation by natively 
#' translating the PLSsemEngine model specification into lavaan-compatible syntax. 
#' This allows researchers to seamlessly evaluate their composite models against 
#' common factor models (CFA).
#'
#' @param measurement_model A named list defining the reflective blocks.
#' @param structural_model An optional list of formulas defining the paths.
#'
#' @return Invisibly returns the generated lavaan syntax string while printing it to the console.
#' @export
export_lavaan_syntax <- function(measurement_model, structural_model = NULL) {
  
  syntax <- c("# --- Measurement Model (CFA) ---")
  
  # Construct measurement equations
  for (construct in names(measurement_model)) {
    items <- measurement_model[[construct]]
    syntax <- c(syntax, paste(construct, "=~", paste(items, collapse = " + ")))
  }
  
  # Construct structural equations if provided
  if (!is.null(structural_model) && length(structural_model) > 0) {
    syntax <- c(syntax, "", "# --- Structural Model ---")
    for (eq in structural_model) {
      lhs <- as.character(eq[[2]])
      rhs <- all.vars(eq[[3]])
      if (length(rhs) > 0) {
        syntax <- c(syntax, paste(lhs, "~", paste(rhs, collapse = " + ")))
      }
    }
  }
  
  final_syntax <- paste(syntax, collapse = "\n")
  
  cat("\n=================================================================\n")
  cat(" LAVAAN SYNTAX GENERATOR (CB-SEM / CFA Integration)\n")
  cat(" Copy and paste this syntax to run models using the 'lavaan' package.\n")
  cat("=================================================================\n\n")
  cat(final_syntax, "\n\n")
  cat("=================================================================\n")
  
  invisible(final_syntax)
}