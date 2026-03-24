# Version: 1.1.0
# Date: 2026-03-24
# Final symmetric release for SoftwareX submission

#################################################
# PLS-SEM ENGINE (Reflective Measurement, Mode A)
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
# PLS-SEM ENGINE (Reflective Measurement, Mode A)
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
# HTMT (Heterotrait-Monotrait Ratio)
#################################################
# HTMT is computed following Henseler et al. (2015).
# It is defined as the ratio of:
# - the mean of heterotrait-heteromethod correlations
# - to the geometric mean of monotrait-heteromethod correlations.
#
# No threshold-based classification is applied.

compute_htmt <- function(data, measurement_model, digits = 2) {
  
  X <- scale(as.matrix(data))
  constructs <- names(measurement_model)
  
  htmt_matrix <- matrix(NA,
                        nrow = length(constructs),
                        ncol = length(constructs))
  
  rownames(htmt_matrix) <- constructs
  colnames(htmt_matrix) <- constructs
  
  for (i in seq_along(constructs)) {
    for (j in seq_along(constructs)) {
      
      if (i >= j) next
      
      items_i <- measurement_model[[constructs[i]]]
      items_j <- measurement_model[[constructs[j]]]
      
      # Heterotrait correlations
      cor_ij <- abs(cor(X[, items_i],
                        X[, items_j],
                        use = "pairwise.complete.obs"))
      
      mean_hetero <- mean(cor_ij, na.rm = TRUE)
      
      # Monotrait correlations (within construct)
      cor_ii <- abs(cor(X[, items_i],
                        use = "pairwise.complete.obs"))
      cor_jj <- abs(cor(X[, items_j],
                        use = "pairwise.complete.obs"))
      
      mean_mono_i <- mean(cor_ii[lower.tri(cor_ii)], na.rm = TRUE)
      mean_mono_j <- mean(cor_jj[lower.tri(cor_jj)], na.rm = TRUE)
      
      htmt_value <- mean_hetero / sqrt(mean_mono_i * mean_mono_j)
      
      htmt_matrix[i, j] <- htmt_value
      htmt_matrix[j, i] <- htmt_value
    }
  }
  
  diag(htmt_matrix) <- NA
  
  round(htmt_matrix, digits)
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

#################################################
# WRAPPER – Tables 1 to 5
#################################################
# High-level wrapper function to execute the full PLS-SEM workflow.
# It delegates estimation to the core engine and organizes the outputs
# into structured, publication-ready tables.
pls_sem <- function(data,
                    measurement_model,
                    structural_model,
                    k = 5,
                    nboot = 500,
                    digits = 2,
                    sign_correction = FALSE,
                    inner_scheme = "factorial") { 
  
  # =====================
  # Base estimate
  # =====================
  engine <- pls_engine(
    data = data,
    measurement_model = measurement_model,
    structural_model = structural_model,
    sign_correction = sign_correction,
    inner_scheme = inner_scheme 
  )
  
  # =====================
  # CMB – Full collinearity VIF
  # =====================
  cmb_table <- compute_cmb_vif(engine$scores)
  
  # =====================
  # Latent scores (table)
  # =====================
  scores_table <- as.data.frame(round(engine$scores, digits))
  scores_table$Case <- seq_len(nrow(scores_table))
  scores_table <- scores_table[, c("Case", setdiff(colnames(scores_table), "Case"))]
  
  # =====================
  # Table 1 – Loadings
  # =====================
  table1 <- data.frame(
    Construct = apply(engine$loadings, 1, function(x)
      colnames(engine$loadings)[which.max(abs(x))]),
    Item = rownames(engine$loadings),
    Loading = round(
      apply(engine$loadings, 1, function(x) max(x, na.rm = TRUE)),
      digits
    ),
    row.names = NULL
  )
  
  # =====================
  # Table 2 – CR, AVE and R2
  # =====================
  table2 <- do.call(
    rbind,
    lapply(names(measurement_model), function(cn) {
      
      items <- measurement_model[[cn]]
      lambda <- engine$loadings[items, cn]
      lambda2 <- lambda^2
      
      den <- (sum(lambda))^2 + sum(1 - lambda2)
      CR <- ifelse(den == 0, NA, (sum(lambda))^2 / den)
      
      AVE <- mean(lambda2)
      
      # Assign R2 only if the construct is endogenous
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
  
  # =====================
  # Table 3 – HTMT
  # =====================
  table3 <- compute_htmt(data, measurement_model, digits)
  
  # =====================
  # Table 4 - Paths and f2 Effect Size
  # =====================
  boot <- bootstrap_paths(
    data = data,
    measurement_model = measurement_model,
    structural_model = structural_model,
    nboot = nboot,
    sign_correction = sign_correction,
    inner_scheme = inner_scheme 
  )
  
  table4 <- do.call(
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
          
          # f2 Effect Size calculation
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
  
  # =====================
  # Table 5 – PLSpredict
  # =====================
  pred <- plspredict(
    data = data,
    measurement_model = measurement_model,
    structural_model = structural_model,
    k = k,
    sign_correction = sign_correction,
    inner_scheme = inner_scheme 
  )
  
  # =====================
  # Final Object Assembly
  # =====================
  model <- list(
    tables = list(
      table1 = table1,
      table2 = table2,
      table3 = table3,
      table4 = table4,
      table5 = pred$table,
      cmb = cmb_table,
      scores = scores_table
    ),
    warnings = list(
      measurement = NULL,
      structural = NULL,
      predictive = pred$warnings
    ),
    meta = list(
      n_obs = nrow(data),
      constructs = names(measurement_model)
    ),
    # Core components saved for plotting and diagnostic functions
    engine = engine, 
    structural_model = structural_model 
  )
  
  class(model) <- "pls_model"
  model
}

# =====================
# Export Scores (helper)
# =====================

export_scores <- function(model, file = "latent_scores.csv") {
  
  scores <- as.data.frame(model$tables$scores)
  scores$Case <- seq_len(nrow(scores))
  scores <- scores[, c("Case", colnames(scores)[colnames(scores) != "Case"])]
  
  write.csv(scores, file, row.names = FALSE)
  invisible(scores)
}

# =====================
# htmt_item_diagnostics
# =====================
# HTMT-guided item diagnostics.
#
# This function is only intended for exploratory diagnostics
# when HTMT thresholds are exceeded.
# It does not suggest item removal or model modification.
# Any model refinement decisions must be theoretically justified.
htmt_item_diagnostics <- function(model, threshold = 0.85, digits = 2) {
  
  htmt <- model$tables$table3
  scores <- model$tables$scores
  loadings <- model$tables$table1
  
  constructs <- colnames(htmt)
  problems <- which(htmt > threshold, arr.ind = TRUE)
  
  if (nrow(problems) == 0) return(NULL)
  
  out <- list()
  
  for (k in seq_len(nrow(problems))) {
    
    c1 <- constructs[problems[k, 1]]
    c2 <- constructs[problems[k, 2]]
    
    items_c1 <- loadings$Item[loadings$Construct == c1]
    items_c2 <- loadings$Item[loadings$Construct == c2]
    
    # cross-item–construct correlation (diagnostic)
    for (it in items_c1) {
      r <- cor(scores[, c2], scores[, c1])
      out[[length(out) + 1]] <- data.frame(
        Construct_A = c1,
        Construct_B = c2,
        Item = it,
        `Construct-level correlation` = round(abs(r), digits)
      )
    }
    
    for (it in items_c2) {
      r <- cor(scores[, c1], scores[, c2])
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
# Indirect effects computed as the product of path coefficients.
#
# No bootstrap inference or mediation classification
# is performed by design.
get_indirect_effects <- function(model, digits = 3) {
  
  t4 <- model$tables$table4
  
  # identify beta column robustly
  beta_col <- grep("Path", colnames(t4), value = TRUE)[1]
  
  out <- list()
  
  # real mediators
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

# =====================
# Helper: Exact Boundary Intersection for Rectangles
# =====================
get_boundary_offset <- function(angle, bw, bh) {
  if (abs(tan(angle)) < bh/bw) {
    return(bw / abs(cos(angle)))
  } else {
    return(bh / abs(sin(angle)))
  }
}

# =====================
# Plot Structural Model (Base R - Pure Rectangles & Export Ready)
# =====================
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

# =====================
# Plot Structural Model with Results (Base R - Configurable & Export Ready)
# =====================
plot_model_results <- function(model, layout = NULL, 
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
  
  structural_model <- model$structural_model
  t2 <- model$tables$table2; t4 <- model$tables$table4 
  
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
    
    x0_arr <- pos_from$x + d_from * cos(angle); y0_arr <- pos_from$y + d_from * sin(angle)
    x1_arr <- pos_to$x - d_to * cos(angle); y1_arr <- pos_to$y - d_to * sin(angle)
    
    arrows(x0 = x0_arr, y0 = y0_arr, x1 = x1_arr, y1 = y1_arr, length = 0.12, lwd = arr_lwd, col = "gray40")
    
    beta_val <- t4$`Path Coefficient (β)`[t4$From == edges$from[i] & t4$To == edges$to[i]]
    mid_x <- (x0_arr + x1_arr)/2; mid_y <- (y0_arr + y1_arr)/2
    lbl <- paste0("β=", beta_val)
    
    # Text Masking - Expanded significantly to completely cut the arrow line
    w <- strwidth(lbl, cex=cex_beta) * 0.75 
    h <- strheight(lbl, cex=cex_beta) * 0.9  
    rect(mid_x - w, mid_y - h, mid_x + w, mid_y + h, col="white", border=NA)
    text(mid_x, mid_y, labels = lbl, cex = cex_beta, font = 3, col = "black")
  }
  
  for(i in 1:nrow(node_pos)) {
    rect(node_pos$x[i]-box_width, node_pos$y[i]-box_height, 
         node_pos$x[i]+box_width, node_pos$y[i]+box_height, 
         col=box_col, border="black", lwd=1.5)
    
    clean_name <- gsub("_", " ", node_pos$name[i])
    r2_val <- t2$R2[t2$Construct == node_pos$name[i]]
    
    # Dynamic spacing based on box height to prevent overlaps regardless of resolution
    if(!is.na(r2_val)) {
      text(node_pos$x[i], node_pos$y[i] + (box_height * 0.25), labels = clean_name, cex=cex_node, font=2, col="black")
      text(node_pos$x[i], node_pos$y[i] - (box_height * 0.35), labels = paste0("(R²=", r2_val, ")"), cex=cex_r2, font=3, col="gray30")
    } else {
      text(node_pos$x[i], node_pos$y[i], labels = clean_name, cex=cex_node, font=2, col="black")
    }
  }
  
  if (save_plot) {
    dev.off()
    message(paste("High-quality plot saved to:", file.path(getwd(), file_name)))
  }
}