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

#################################################
# PLSsemEngine (Reflective Measurement, Mode A)
#################################################
pls_engine <- function(data,
                       measurement_model,
                       structural_model,
                       max_iter = 300,
                       tol = 1e-6,
                       sign_correction = FALSE,
                       inner_scheme = "factorial") { 
  
  X_raw <- as.matrix(data)
  X <- scale(X_raw)
  
  if (any(is.na(X))) {
    stop("Scaling produced NA values. Check for zero-variance indicators.")
  }
  
  blocks <- measurement_model
  constructs <- names(blocks)
  
  outer_weights <- vector("list", length(blocks))
  names(outer_weights) <- constructs
  
  # --- 1. Adjacency Matrix for Inner Model
  # Construct the connection network (S) based on the structural model
  C <- matrix(0, nrow = length(constructs), ncol = length(constructs), 
              dimnames = list(constructs, constructs))
  for (f in structural_model) {
    lhs <- as.character(f[[2]])
    rhs <- all.vars(f[[3]])
    for (r in rhs) {
      C[r, lhs] <- 1
      C[lhs, r] <- 1 # Inner relationships in PLS are bidirectional for weight calculation
    }
  }
  
  # --- Latent variable score initialization
  scores <- sapply(blocks, function(items) {
    rowMeans(X[, items, drop = FALSE], na.rm = TRUE)
  })
  scores <- scale(scores)
  colnames(scores) <- constructs
  
  # --- PLS iterative algorithm (Mode A)
  for (iter in seq_len(max_iter)) {
    
    scores_old <- scores
    
    # --- 2. Inner Approximation
    # Constructs communicate based on the chosen inner weighting scheme
    R <- cor(scores, use = "pairwise.complete.obs")
    if (inner_scheme == "centroid") {
      E <- sign(R) * C
    } else {
      E <- R * C # Factorial scheme by default
    }
    
    # Inner approximation (Z)
    Z <- scores %*% E
    
    # --- 3. Outer Approximation (Weight Update Mode A)
    for (j in seq_along(blocks)) {
      items <- blocks[[j]]
      
      # Correlation between indicators and the inner approximation
      w <- apply(X[, items, drop = FALSE], 2, function(x) cor(x, Z[, j], use = "pairwise.complete.obs"))
      
      if (all(is.na(w))) next
      
      # --- 4. Fix Composite Variance to 1 
      # Normalize by dividing by the standard deviation to ensure unit variance
      unscaled_score <- as.numeric(X[, items, drop = FALSE] %*% w)
      w <- w / sd(unscaled_score, na.rm = TRUE)
      
      outer_weights[[j]] <- w
      scores[, j] <- X[, items, drop = FALSE] %*% w
    }
    
    scores <- scale(scores)
    
    # Convergence criterion
    if (max(abs(scores - scores_old), na.rm = TRUE) < tol) break
  }
  
  # --- Sign Indeterminacy Correction
  # Disabled by default to prevent bootstrap truncation (Rönkkö et al., 2015)
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
  
  # --- Indicator loadings
  loadings <- matrix(NA, nrow = sum(lengths(blocks)), ncol = length(blocks))
  rownames(loadings) <- unlist(blocks)
  colnames(loadings) <- constructs
  
  for (j in seq_along(blocks)) {
    for (item in blocks[[j]]) {
      loadings[item, j] <- cor(X[, item], scores[, j], use = "pairwise.complete.obs")
    }
  }
  
  # --- Structural model estimation
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
    inner_scheme = inner_scheme # Save metadata
  )
}

#################################################
# BOOTSTRAP STRUCTURAL MODEL (OLS on scores)
#################################################
bootstrap_paths <- function(data,
                            measurement_model,
                            structural_model,
                            nboot = 500,
                            seed = 123,
                            sign_correction = FALSE,
                            inner_scheme = "factorial") { 
  
  set.seed(seed)
  n <- nrow(data)
  results <- list()
  
  for (b in seq_len(nboot)) {
    
    idx <- sample(seq_len(n), replace = TRUE)
    data_b <- data[idx, , drop = FALSE]
    
    # Pass parameters to the core engine
    engine_b <- pls_engine(
      data = data_b,
      measurement_model = measurement_model,
      structural_model = structural_model,
      sign_correction = sign_correction,
      inner_scheme = inner_scheme 
    )
    
    for (f in structural_model) {
      lhs <- as.character(f[[2]])
      rhs <- all.vars(f[[3]])
      betas <- engine_b$paths[[lhs]]
      
      for (i in seq_along(rhs)) {
        results[[length(results) + 1]] <- data.frame(
          Path = paste(rhs[i], "→", lhs),
          Beta = betas[i]
        )
      }
    }
  }
  
  do.call(rbind, results)
}

#################################################
# COMMON METHOD BIAS (Full collinearity VIF)
#################################################
# Full collinearity VIF following Kock (2015).
# Values above 3.3 may indicate potential common method bias.
# No automatic classification is applied.

compute_cmb_vif <- function(scores, digits = 2) {
  
  constructs <- colnames(scores)
  vif <- numeric(length(constructs))
  names(vif) <- constructs
  
  for (i in seq_along(constructs)) {
    
    y <- scores[, i]
    x <- scores[, -i, drop = FALSE]
    
    fit <- lm(y ~ x)
    r2 <- summary(fit)$r.squared
    vif[i] <- 1 / (1 - r2)
  }
  
  data.frame(
    Construct = constructs,
    VIF = round(vif, digits),
    row.names = NULL
  )
}

#################################################
# HTMT & HTMT2 (Heterotrait-Monotrait Ratios)
#################################################
# Computes HTMT (Henseler et al., 2015) and HTMT2 (Roemer et al., 2021).
# HTMT2 uses the geometric mean, relaxing the tau-equivalence assumption, 
# making it a consistent estimator for congeneric models.
# No threshold-based classification is applied.

compute_htmt_metrics <- function(data, measurement_model, digits = 2) {
  
  X <- scale(as.matrix(data))
  constructs <- names(measurement_model)
  
  htmt_matrix <- matrix(NA, nrow = length(constructs), ncol = length(constructs), dimnames = list(constructs, constructs))
  htmt2_matrix <- matrix(NA, nrow = length(constructs), ncol = length(constructs), dimnames = list(constructs, constructs))
  
  # Helper for geometric mean (strictly positive values)
  geom_mean <- function(x) {
    x_pos <- x[x > 0 & !is.na(x)]
    if (length(x_pos) == 0) return(NA)
    exp(mean(log(x_pos)))
  }
  
  for (i in seq_along(constructs)) {
    for (j in seq_along(constructs)) {
      
      if (i >= j) next
      
      items_i <- measurement_model[[constructs[i]]]
      items_j <- measurement_model[[constructs[j]]]
      
      # Heterotrait correlations (absolute)
      cor_ij <- abs(cor(X[, items_i, drop = FALSE], X[, items_j, drop = FALSE], use = "pairwise.complete.obs"))
      
      # Monotrait correlations (within construct, absolute)
      cor_ii <- abs(cor(X[, items_i, drop = FALSE], use = "pairwise.complete.obs"))
      cor_jj <- abs(cor(X[, items_j, drop = FALSE], use = "pairwise.complete.obs"))
      
      cor_ii_vec <- cor_ii[lower.tri(cor_ii)]
      cor_jj_vec <- cor_jj[lower.tri(cor_jj)]
      
      # HTMT (Arithmetic mean)
      mean_hetero_a <- mean(cor_ij, na.rm = TRUE)
      mean_mono_i_a <- mean(cor_ii_vec, na.rm = TRUE)
      mean_mono_j_a <- mean(cor_jj_vec, na.rm = TRUE)
      
      htmt_val <- mean_hetero_a / sqrt(mean_mono_i_a * mean_mono_j_a)
      htmt_matrix[i, j] <- htmt_val
      htmt_matrix[j, i] <- htmt_val
      
      # HTMT2 (Geometric mean)
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

#################################################
# MODEL FIT (SRMR, d_ULS, d_G)
#################################################
# Implements model fit metrics following Henseler et al. (2014).
# Calculates SRMR, exact d_ULS, and exact d_G based on the 
# saturated model-implied correlation matrix.
compute_model_fit <- function(engine, data, measurement_model, digits = 3) {
  
  items <- unlist(measurement_model)
  S <- cor(data[, items, drop = FALSE], use = "pairwise.complete.obs")
  
  # Build strict loading matrix (0 for cross-loadings)
  constructs <- names(measurement_model)
  Lambda <- matrix(0, nrow = length(items), ncol = length(constructs),
                   dimnames = list(items, constructs))
  
  for (cn in constructs) {
    inds <- measurement_model[[cn]]
    Lambda[inds, cn] <- engine$loadings[inds, cn]
  }
  
  # Construct correlation matrix (saturated inner model)
  R <- cor(engine$scores, use = "pairwise.complete.obs")
  
  # Model-implied correlation matrix
  Sigma_hat <- Lambda %*% R %*% t(Lambda)
  diag(Sigma_hat) <- 1
  
  # Number of off-diagonal elements
  p <- nrow(S)
  lower_tri_idx <- lower.tri(S)
  e <- S[lower_tri_idx] - Sigma_hat[lower_tri_idx]
  
  # 1. SRMR (Standardized Root Mean Square Residual)
  srmr <- sqrt(mean(e^2))
  
  # 2. d_ULS (Squared Euclidean distance)
  d_uls <- sum(e^2)
  
  # 3. d_G (Geodesic distance)
  # Calculated as half the sum of squared logs of the eigenvalues of S * Sigma_hat^-1
  d_g <- NA
  tryCatch({
    eig_vals <- eigen(S %*% solve(Sigma_hat), only.values = TRUE)$values
    eig_vals <- Re(eig_vals)
    eig_vals <- eig_vals[eig_vals > 0] # Keep only positive real parts
    d_g <- 0.5 * sum((log(eig_vals))^2)
  }, error = function(e) {
    d_g <<- NA
  })
  
  data.frame(
    Metric = c("SRMR", "d_ULS", "d_G"),
    Value = c(round(srmr, digits), round(d_uls, digits), round(d_g, digits)),
    stringsAsFactors = FALSE
  )
}

#################################################
# PLSpredict
#################################################
# PLSpredict implementation following Shmueli et al. (2019).
#
# Predictive performance is assessed using k-fold cross-validation,
# comparing PLS-based predictions against a linear model benchmark.
#
# Predictive power is reported without automatic decision rules.
#################################################
plspredict <- function(data,
                       measurement_model,
                       structural_model,
                       k = 5,
                       seed = 123,
                       sign_correction = FALSE,
                       inner_scheme = "factorial") { 
  
  set.seed(seed)
  n <- nrow(data)
  folds <- sample(rep(1:k, length.out = n))
  
  results <- list()
  
  for (fold in seq_len(k)) {
    
    train <- data[folds != fold, , drop = FALSE]
    test  <- data[folds == fold,  , drop = FALSE]
    
    # Pass parameters to the core engine
    model <- pls_engine(train, measurement_model, structural_model, 
                        sign_correction = sign_correction,
                        inner_scheme = inner_scheme) 
    
    # Standardize test data using training mean and sd
    X_train <- scale(as.matrix(train))
    X_test  <- scale(as.matrix(test),
                     center = attr(X_train, "scaled:center"),
                     scale  = attr(X_train, "scaled:scale"))
    
    # Compute latent scores for test data using training outer weights
    test_scores <- matrix(NA,
                          nrow = nrow(X_test),
                          ncol = length(measurement_model))
    
    colnames(test_scores) <- names(measurement_model)
    
    for (j in seq_along(measurement_model)) {
      
      items <- measurement_model[[j]]
      w <- model$outer_weights[[j]]
      
      test_scores[, j] <- X_test[, items, drop = FALSE] %*% w
    }
    
    test_scores <- scale(test_scores)
    
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
        
        # PLS prediction (standardized scale)
        y_pls_std <- as.numeric(y_hat * lambda)
        
        # Re-scale to original indicator scale using TRAIN statistics
        mean_train <- mean(train[[ind]], na.rm = TRUE)
        sd_train   <- sd(train[[ind]], na.rm = TRUE)
        
        y_pls <- y_pls_std * sd_train + mean_train
        
        y_obs <- test[[ind]]
        
        rmse_pls <- sqrt(mean((y_obs - y_pls)^2, na.rm = TRUE))
        
        # LM benchmark (latent-score based)
        
        # Training data for LM
        df_train <- data.frame(
          y = train[[ind]],
          model$scores[, rhs, drop = FALSE]
        )
        
        lm_fit <- lm(y ~ ., data = df_train)
        
        # Test data for LM (use test_scores)
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
  
  # --- Table 5
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
  
  table5$Q2_predict <- ifelse(
    table5$RMSE_LM == 0,
    NA,
    1 - (table5$RMSE_PLS^2 / table5$RMSE_LM^2)
  )
  
  # Warning is issued only when Q²_predict is negative,
  # indicating that PLS predictions are outperformed
  # by the linear benchmark.
  warnings <- character(0)
  if (any(table5$Q2_predict < 0, na.rm = TRUE)) {
    warnings <- paste(
      "Negative Q²_predict detected:",
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
# WRAPPER – High-Level Analysis Workflow (Version 1.2.0)
#################################################################
# Main function of the PLSsemEngine. It coordinates the estimation,
# assessment, and predictive evaluation stages, returning results
# in a structured, publication-ready S3 object.

#' Estimate a PLS-SEM Model
#'
#' @description This is the main function of the PLSsemEngine. It estimates 
#' latent variable scores, measurement models, and structural paths.
#' @param data A data frame containing the indicators.
#' @param measurement_model A named list defining the reflective blocks.
#' @param structural_model A list of formulas defining the paths.
#' @param k Number of folds for cross-validation (default is 5).
#' @param nboot Number of bootstrap resamples (default is 500).
#' @param digits Number of decimal places for output tables.
#' @param sign_correction Logical; if TRUE, applies deterministic sign alignment.
#' @param inner_scheme The weighting scheme for the inner model ("factorial" or "centroid").
#' @return A list of class 'pls_model' containing structured tables and metadata.
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
  # Executes the iterative PLS algorithm (Mode A) to obtain latent scores
  engine <- pls_engine(
    data = data,
    measurement_model = measurement_model,
    structural_model = structural_model,
    sign_correction = sign_correction,
    inner_scheme = inner_scheme 
  )
  
  # =============================================================
  # 2. MEASUREMENT MODEL ASSESSMENT (Manuscript Table 1)
  # =============================================================
  # Extract item loadings
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
  
  # Calculate Construct Reliability (CR), AVE, and R-squared
  t1_constructs <- do.call(
    rbind,
    lapply(names(measurement_model), function(cn) {
      items <- measurement_model[[cn]]
      lambda <- engine$loadings[items, cn]
      lambda2 <- lambda^2
      
      # Composite Reliability (CR) calculation
      den <- (sum(lambda))^2 + sum(1 - lambda2)
      CR <- ifelse(den == 0, NA, (sum(lambda))^2 / den)
      
      # Average Variance Extracted (AVE) calculation
      AVE <- mean(lambda2)
      
      # Retrieve R2 if construct is endogenous
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
  
  # Merge item and construct metrics for Table 1
  measurement_results <- merge(t1_items, t1_constructs, by = "Construct", all.x = TRUE)
  measurement_results <- measurement_results[order(measurement_results$Construct, measurement_results$Item), ]
  rownames(measurement_results) <- NULL
  
  # =============================================================
  # 3. DISCRIMINANT VALIDITY (Manuscript Table 3)
  # =============================================================
  # Computes HTMT and HTMT2 (congeneric models) as requested by Reviewer 4
  htmt_results <- compute_htmt_metrics(data, measurement_model, digits)
  
  # =============================================================
  # 4. STRUCTURAL MODEL & INFERENCE (Manuscript Table 4)
  # =============================================================
  # Execute bootstrap resampling
  boot <- bootstrap_paths(
    data = data,
    measurement_model = measurement_model,
    structural_model = structural_model,
    nboot = nboot,
    sign_correction = sign_correction,
    inner_scheme = inner_scheme 
  )
  
  # Compile structural path coefficients, confidence intervals, and f2 effect sizes
  structural_results <- do.call(
    rbind,
    lapply(structural_model, function(f) {
      
      lhs <- as.character(f[[2]])
      rhs <- all.vars(f[[3]])
      
      betas <- engine$paths[[lhs]]
      r2_incl <- engine$r2[lhs]
      
      do.call(
        rbind,
        lapply(seq_along(rhs), function(i) {
          
          path_name <- paste(rhs[i], "→", lhs)
          bdist <- boot$Beta[boot$Path == path_name]
          
          # f2 Effect Size calculation (R2 increase)
          rhs_excl <- setdiff(rhs, rhs[i])
          if (length(rhs_excl) == 0) {
            r2_excl <- 0
          } else {
            fit_excl <- lm(engine$scores[, lhs] ~ engine$scores[, rhs_excl, drop = FALSE])
            r2_excl <- summary(fit_excl)$r.squared
          }
          
          # Mathematical protection against R2 values >= 1
          f2_val <- ifelse(r2_incl >= 1, NA, (r2_incl - r2_excl) / (1 - r2_incl))
          
          data.frame(
            From = rhs[i],
            To = lhs,
            `Path Coefficient (β)` = round(betas[i], digits),
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
  # 5. PREDICTIVE EVALUATION (Manuscript Table 5)
  # =============================================================
  # Executes k-fold cross-validation (PLSpredict)
  pred <- plspredict(
    data = data,
    measurement_model = measurement_model,
    structural_model = structural_model,
    k = k,
    sign_correction = sign_correction,
    inner_scheme = inner_scheme 
  )
  
  # =============================================================
  # 6. COMPLEMENTARY DIAGNOSTICS (Reviewer-requested features)
  # =============================================================
  # Full Collinearity VIF for Common Method Bias
  cmb_table <- compute_cmb_vif(engine$scores)
  
  # Global Model Fit indices (SRMR, d_ULS, d_G)
  fit_table <- compute_model_fit(engine, data, measurement_model, digits)
  
  # =============================================================
  # 7. FINAL OBJECT ASSEMBLY
  # =============================================================
  # Organized into descriptive elements to satisfy 'usability' concerns
  model <- list(
    measurement_model = measurement_results,     # Corresponds to Table 1
    discriminant_validity = htmt_results,        # Corresponds to Table 3 (HTMT & HTMT2)
    structural_model = structural_results,       # Corresponds to Table 4
    predictive_relevance = pred$table,           # Corresponds to Table 5
    diagnostics = list(
      common_method_bias = cmb_table,            # Full Collinearity VIF
      global_fit = fit_table                     # SRMR, d_ULS, d_G
    ),
    latent_scores = as.data.frame(round(engine$scores, digits)),
    raw_engine = engine,                         # Stored for plotting/internal diagnostics
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
# =====================
# Export Scores (helper)
# =====================
#' @export
export_scores <- function(model, file = "latent_scores.csv") {
  
  # CAMBIO: Ahora accede directamente a model$latent_scores
  scores <- as.data.frame(model$latent_scores)
  
  # Se añade la columna de casos para mejorar la trazabilidad
  scores$Case <- seq_len(nrow(scores))
  scores <- scores[, c("Case", setdiff(colnames(scores), "Case"))]
  
  write.csv(scores, file, row.names = FALSE)
  invisible(scores)
}
# =====================
# htmt_item_diagnostics
# =====================
# HTMT-guided item diagnostics.
# Updated for v1.2.0 compatibility.
#' @export
htmt_item_diagnostics <- function(model, threshold = 0.85, digits = 2) {
  
  # CAMBIOS: Actualización de rutas para coincidir con el nuevo Wrapper
  # Acceso a HTMT dentro de discriminant_validity
  htmt <- model$discriminant_validity$HTMT 
  # Acceso directo a latent_scores y measurement_model
  scores <- model$latent_scores
  loadings <- model$measurement_model
  
  constructs <- colnames(htmt)
  problems <- which(htmt > threshold, arr.ind = TRUE)
  
  if (nrow(problems) == 0) return(NULL)
  
  out <- list()
  
  for (k in seq_len(nrow(problems))) {
    
    c1 <- constructs[problems[k, 1]]
    c2 <- constructs[problems[k, 2]]
    
    items_c1 <- loadings$Item[loadings$Construct == c1]
    items_c2 <- loadings$Item[loadings$Construct == c2]
    
    # Diagnóstico de correlación cruzada ítem-constructo
    for (it in items_c1) {
      # Usamos pairwise.complete.obs por seguridad estadística
      r <- cor(scores[, c2], scores[, c1], use = "pairwise.complete.obs")
      out[[length(out) + 1]] <- data.frame(
        Construct_A = c1,
        Construct_B = c2,
        Item = it,
        `Construct-level correlation` = round(abs(r), digits)
      )
    }
    
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
# =====================
# Indirect effects
# =====================
# Indirect effects computed as the product of path coefficients[cite: 85, 228].
# Consistent with the engine's philosophy, no automated mediation labels 
# are provided, leaving interpretation to the researcher[cite: 86, 87, 230].
#' @export
get_indirect_effects <- function(model, digits = 3) {
  
  # CAMBIO: Ahora accede directamente a model$structural_model (antes table4)
  t4 <- model$structural_model
  
  if (is.null(t4)) return(NULL)
  
  # Identificación robusta de la columna de coeficientes beta
  beta_col <- grep("Path", colnames(t4), value = TRUE)[1]
  
  out <- list()
  
  # Identificación de mediadores reales (constructos que son tanto origen como destino)
  mediators <- intersect(unique(t4$To), unique(t4$From))
  
  for (m in mediators) {
    
    incoming <- t4[t4$To == m, , drop = FALSE]
    outgoing <- t4[t4$From == m, , drop = FALSE]
    
    if (nrow(incoming) == 0 || nrow(outgoing) == 0) next
    
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
# =====================
# References table
# =====================
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
      "≥ 0.70",
      "≥ 0.70",
      "≥ 0.50",
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

# =============================================================
# Helper: Exact Boundary Intersection for Rectangles
# =================--------------------------------------------
# Calculates the offset for arrow positioning to ensure they 
# touch the rectangle boundaries precisely.
#' @export
get_boundary_offset <- function(angle, bw, bh) {
  if (abs(tan(angle)) < bh/bw) {
    return(bw / abs(cos(angle)))
  } else {
    return(bh / abs(sin(angle)))
  }
}
# =============================================================
# Plot Structural Model (Base R - Pure Rectangles & Export Ready)
# =============================================================
#' @export
plot_structural_model <- function(structural_model, layout = NULL, 
                                  box_width = 1.1, box_height = 0.3, 
                                  cex_node = 0.9, arr_lwd = 1.2, 
                                  box_col = "white",
                                  save_plot = FALSE, file_name = "structural_model.png",
                                  width = 2500, height = 1500, res = 300) {
  
  # [El resto del código se mantiene idéntico]
  if (save_plot) {
    if (capabilities("cairo")) {
      png(filename = file_name, width = width, height = height, res = res, type = "cairo")
    } else {
      png(filename = file_name, width = width, height = height, res = res)
    }
  }
  
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

# =============================================================
# Plot Structural Model with Results (Base R - Configurable)
# =============================================================
# Updated for v1.2.0: Aligned with new modular output structure.
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
  
  # CAMBIOS: Actualización de rutas internas
  structural_model <- model$specification$structural
  t1 <- model$measurement_model # Anteriormente table1
  t4 <- model$structural_model  # Anteriormente table4
  
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
  
  # Dibujo de flechas y Coeficientes Beta
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
    
    # CAMBIO: Identificación robusta de la columna de coeficientes
    beta_col <- grep("Path", colnames(t4), value = TRUE)[1]
    beta_val <- t4[[beta_col]][t4$From == edges$from[i] & t4$To == edges$to[i]]
    
    mid_x <- (x0_arr + x1_arr)/2; mid_y <- (y0_arr + y1_arr)/2
    lbl <- paste0("beta=", beta_val)
    
    w <- strwidth(lbl, cex=cex_beta) * 0.75 
    h <- strheight(lbl, cex=cex_beta) * 0.9  
    rect(mid_x - w, mid_y - h, mid_x + w, mid_y + h, col="white", border=NA)
    text(mid_x, mid_y, labels = lbl, cex = cex_beta, font = 3, col = "black")
  }
  
  # Dibujo de Rectángulos y R2
  for(i in 1:nrow(node_pos)) {
    rect(node_pos$x[i]-box_width, node_pos$y[i]-box_height, 
         node_pos$x[i]+box_width, node_pos$y[i]+box_height, 
         col=box_col, border="black", lwd=1.5)
    
    clean_name <- gsub("_", " ", node_pos$name[i])
    # CAMBIO: Búsqueda del R2 en el nuevo objeto t1
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
# =============================================================
# METHODOLOGICAL ASSESSMENT (Optional Interpretation Layer)
# =============================================================
# Provides a narrative diagnostic layer based on established 
# literature (Hair et al., Henseler et al.) without automated decision-making.
# Updated for v1.2.0 to match the descriptive output structure.
#' @export
interpret_model <- function(model) {
  
  if (!inherits(model, "pls_model")) {
    stop("Object must be of class 'pls_model'")
  }
  
  cat("\n=================================================================\n")
  cat(" OPTIONAL METHODOLOGICAL ASSESSMENT (HEURISTIC AIDS)\n")
  cat(" Note: The following classifications are based on standard \n")
  cat(" literature thresholds. They are provided to contextualise \n")
  cat(" results, not to replace researcher judgment.\n")
  cat("=================================================================\n\n")
  
  # 1. Measurement Model (Reliability & AVE)
  cat("--- 1. Reflective Measurement Model (Hair et al., 2017) ---\n")
  
  # CAMBIO: Ahora accede directamente a model$measurement_model
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
  
  # 2. Discriminant Validity (HTMT & HTMT2)
  cat("\n--- 2. Discriminant Validity (Henseler et al., 2015; Roemer et al., 2021) ---\n")
  
  # CAMBIO: Acceso a la lista discriminant_validity
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
  
  # 3. Collinearity (Common Method Bias)
  cat("\n--- 3. Full Collinearity VIF (Kock, 2015) ---\n")
  
  # CAMBIO: Acceso a diagnostics$common_method_bias
  cmb <- model$diagnostics$common_method_bias
  
  cmb_issues <- cmb$Construct[!is.na(cmb$VIF) & cmb$VIF > 3.3]
  
  if(length(cmb_issues) == 0) {
    cat(" [*] CMB VIF: All constructs meet the <= 3.3 guideline.\n") 
  } else {
    cat(" [!] CMB VIF: Constructs > 3.3 guideline ->", paste(cmb_issues, collapse=", "), "\n")
  }
  
  cat("\n=================================================================\n")
}

# =====================
# METHODOLOGICAL INTEGRATION (CB-SEM / CFA)
# =====================
# Bridges the workflow between variance-based PLS-SEM and 
# covariance-based SEM (CB-SEM) or Confirmatory Factor Analysis (CFA).
# Automatically translates the native engine specification into lavaan syntax,
# fostering comparative methodological education and cross-validation.

#' Export to Lavaan Syntax
#'
#' @description Translates the PLSsemEngine model into lavaan-compatible syntax
#' for CFA or CB-SEM validation.
#' @param measurement_model The measurement model list.
#' @param structural_model The structural model list.
#' @export
export_lavaan_syntax <- function(measurement_model, structural_model = NULL) {
  
  syntax <- c("# --- Measurement Model (CFA) ---")
  
  for (construct in names(measurement_model)) {
    items <- measurement_model[[construct]]
    syntax <- c(syntax, paste(construct, "=~", paste(items, collapse = " + ")))
  }
  
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