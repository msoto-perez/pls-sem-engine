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
pls_engine <- function(data,
                       measurement_model,
                       structural_model,
                       max_iter = 300,
                       tol = 1e-6) {
  
  X <- scale(as.matrix(data))
  blocks <- measurement_model
  constructs <- names(blocks)
  
  # --- Latent variable score initialization
  scores <- sapply(blocks, function(items) {
    rowMeans(X[, items, drop = FALSE], na.rm = TRUE) # Agregado na.rm por seguridad
  })
  scores <- scale(scores)
  colnames(scores) <- constructs
  
  # --- PLS iterative algorithm (Mode A)
  for (iter in seq_len(max_iter)) {
    
    scores_old <- scores
    
    for (j in seq_along(blocks)) {
      items <- blocks[[j]]
      y <- scores[, j]
      
      # Update outer weights using indicator-score correlations.
      # Missing correlations are handled using pairwise deletion to preserve data.
      w <- apply(X[, items, drop = FALSE], 2, function(x) cor(x, y, use = "pairwise.complete.obs"))
      
      if (all(is.na(w))) next
      
      w <- w / sqrt(sum(w^2))
      scores[, j] <- X[, items, drop = FALSE] %*% w
    }
    
    scores <- scale(scores)
    
    # Convergence criterion
    if (max(abs(scores - scores_old), na.rm = TRUE) < tol) break
  }
  
  # --- Sign Indeterminacy Correction
  # Ensures that the dominant indicator (highest absolute loading) within each block
  # has a positive sign. If negative, the latent score sign is flipped.
  for (j in seq_along(blocks)) {
    items <- blocks[[j]]
    # Compute temporary loadings to check signs (ignoring NAs)
    temp_loadings <- cor(X[, items, drop=FALSE], scores[, j], use = "pairwise.complete.obs")
    
    # Identify the indicator with the maximum absolute loading
    max_ind <- which.max(abs(temp_loadings))
    
    # If the dominant indicator's loading is negative, flip the score sign
    if (!is.na(temp_loadings[max_ind]) && temp_loadings[max_ind] < 0) {
      scores[, j] <- scores[, j] * -1
    }
  }
  
  # --- Indicator loadings
  # Loadings are computed as correlations between standardized
  # indicators and their corresponding latent variable scores.
  loadings <- matrix(NA,
                     nrow = sum(lengths(blocks)),
                     ncol = length(blocks))
  rownames(loadings) <- unlist(blocks)
  colnames(loadings) <- constructs
  
  for (j in seq_along(blocks)) {
    for (item in blocks[[j]]) {
      # Loadings computed with pairwise deletion for missing data stability
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
    r2 = r2
  )
}
#################################################
# BOOTSTRAP STRUCTURAL MODEL (OLS on scores)
#################################################
# Bootstrap procedure for structural path coefficients.
#
# Bootstrap samples are generated at the observation level.
# For each resample, the full PLS-SEM model is re-estimated,
# and path coefficients are stored for inference.
#
# No bootstrap-based classification or hypothesis testing
# is performed automatically.
bootstrap_paths <- function(data,
                            measurement_model,
                            structural_model,
                            nboot = 500,
                            seed = 123) {
  
  set.seed(seed)
  n <- nrow(data)
  results <- list()
  
  for (b in seq_len(nboot)) {
    
    idx <- sample(seq_len(n), replace = TRUE)
    data_b <- data[idx, , drop = FALSE]
    
    engine_b <- pls_engine(
      data = data_b,
      measurement_model = measurement_model,
      structural_model = structural_model
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
# Common Method Bias assessment using full collinearity VIF
# as proposed by Kock (2015).
#
# Each latent variable is regressed on all others to compute
# variance inflation factors (VIF).

compute_cmb_vif <- function(scores) {
  
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
    VIF = round(vif, 2),
    CMB = ifelse(vif < 3.3, "Below recommended threshold", "Above recommended threshold"),
    row.names = NULL
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

plspredict <- function(data,
                       measurement_model,
                       structural_model,
                       k = 5,
                       seed = 123) {
  
  set.seed(seed)
  n <- nrow(data)
  folds <- sample(rep(1:k, length.out = n))
  
  results <- list()
  
  for (fold in seq_len(k)) {
    
    train <- data[folds != fold, , drop = FALSE]
    test  <- data[folds == fold,  , drop = FALSE]
    
    model <- pls_engine(train, measurement_model, structural_model)
    
    for (f in structural_model) {
      
      lhs <- as.character(f[[2]])
      rhs <- all.vars(f[[3]])
      if (length(rhs) == 0) next
      
      y_hat <- as.matrix(model$scores[, rhs, drop = FALSE]) %*%
        model$paths[[lhs]]
      
      inds <- measurement_model[[lhs]]
      
      for (ind in inds) {
        
        lambda <- model$loadings[ind, lhs]
        if (is.na(lambda)) next
        
        # PLS prediction
        y_pls <- as.numeric(y_hat * lambda)
        y_obs <- test[[ind]]
        rmse_pls <- sqrt(mean((y_obs - y_pls)^2, na.rm = TRUE))
        
        # LM benchmark (latent-score based)
        df_train <- data.frame(
          y = train[[ind]],
          model$scores[, rhs, drop = FALSE]
        )
        lm_fit <- lm(y ~ ., data = df_train)
        
        df_test <- data.frame(model$scores[, rhs, drop = FALSE])
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
  
  table5$Q2_predict <- 1 - (table5$RMSE_PLS^2 / table5$RMSE_LM^2)
  
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
# High-level wrapper that produces paper-ready tables
# commonly reported in PLS-SEM studies.
#
# The function returns numerical results only.
# Interpretation and model evaluation remain the
# responsibility of the researcher.
pls_sem <- function(data,
                    measurement_model,
                    structural_model,
                    k = 5,
                    nboot = 500,
                    digits = 2) {
  
  # =====================
  # Base estimate
  # =====================
  engine <- pls_engine(
    data = data,
    measurement_model = measurement_model,
    structural_model = structural_model
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
  # Table 2 – CR y AVE
  # =====================
  table2 <- do.call(
    rbind,
    lapply(names(measurement_model), function(cn) {
      
      items <- measurement_model[[cn]]
      lambda <- engine$loadings[items, cn]
      lambda2 <- lambda^2
      
      CR <- (sum(lambda))^2 / ((sum(lambda))^2 + sum(1 - lambda2))
      AVE <- mean(lambda2)
      
      data.frame(
        Construct = cn,
        `Composite Reliability (CR)` = round(CR, digits),
        AVE = round(AVE, digits),
        row.names = NULL
      )
    })
  )
  
  # =====================
  # Table 3 – HTMT (proxy)
  # =====================
  table3 <- round(cor(engine$scores), digits)
  diag(table3) <- NA
  
  # =====================
  # Table 4
  # =====================
  boot <- bootstrap_paths(
    data = data,
    measurement_model = measurement_model,
    structural_model = structural_model,
    nboot = nboot
  )
  
  table4 <- do.call(
    rbind,
    lapply(structural_model, function(f) {
      
      lhs <- as.character(f[[2]])
      rhs <- all.vars(f[[3]])
      
      betas <- engine$paths[[lhs]]
      
      do.call(
        rbind,
        lapply(seq_along(rhs), function(i) {
          
          path_name <- paste(rhs[i], "→", lhs)
          bdist <- boot$Beta[boot$Path == path_name]
          
          data.frame(
            From = rhs[i],
            To = lhs,
            `Path Coefficient (β)` = round(betas[i], digits),
            CI_low = round(quantile(bdist, 0.025), digits),
            CI_high = round(quantile(bdist, 0.975), digits),
            R2 = round(engine$r2[lhs], digits),
            row.names = NULL
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
    k = k
  )
  
  # =====================
  # Final object
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
    )
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
