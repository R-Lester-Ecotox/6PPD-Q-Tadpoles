median_aov <- function(data, cols, class_col, nboot = 1000) {
  # ---- Error checking ----
  if (!is.data.frame(data)) {
    stop("`data` must be a data.frame.", call. = FALSE)
  }
  
  if (!requireNamespace("WRS2", quietly = TRUE)) {
    stop("Package 'WRS2' is required but not installed.", call. = FALSE)
  }
  
  # Check that all specified columns exist
  missing_cols <- setdiff(cols, colnames(data))
  if (length(missing_cols) > 0) {
    stop(
      "The following `cols` are not found in `data`: ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }
  
  for (arg_name in c("class_col")) {
    arg_val <- get(arg_name)
    if (length(arg_val) != 1 || !arg_val %in% colnames(data)) {
      stop(
        sprintf("`%s` ('%s') not found in `data`.", arg_name, arg_val),
        call. = FALSE
      )
    }
  }
  
  # Ensure grouping columns are factors
  if (!is.factor(data[[class_col]])) {
    message(
      sprintf("Coercing '%s' to factor.", class_col)
    )
    data[[class_col]] <- as.factor(data[[class_col]])
  }
  
  # Check that metabolite columns are numeric
  non_numeric <- cols[!vapply(data[cols], is.numeric, logical(1))]
  if (length(non_numeric) > 0) {
    stop(
      "The following `cols` are not numeric: ",
      paste(non_numeric, collapse = ", "),
      call. = FALSE
    )
  }
  
  # Check minimum group sizes
  n_class <- nlevels(data[[class_col]])
  if (n_class < 2) {
    stop("`class_col` must have at least 2 levels.", call. = FALSE)
  }
  
  # ---- Record the call and formula ----
  
  matched_call <- match.call()
  model_formula <- stats::as.formula(
    paste("y ~", class_col)
  )
  message("Model formula: ", deparse(model_formula))
  
  # ---- Fit models ----
  
  mat   <- as.matrix(data[, cols])
  class <- data[[class_col]]
  
  res <- lapply(seq_along(cols), function(i) {
    
    y <- mat[, i]
    
    fit <- tryCatch(
      WRS2::med1way(
        y ~ class,
        iter = nboot
      ),
      error = function(e) {
        warning(
          sprintf(
            "med1way failed for '%s': %s", cols[i], conditionMessage(e)
          ),
          call. = FALSE
        )
        return(NULL)
      }
    )
    
    if (is.null(fit)) {
      return(
        data.frame(
          variable             = cols[i],
          test                 = NA_real_,
          crit_val             = NA_real_,
          p_val                = NA_real_
        )
      )
    }
    
    data.frame(
      variable             = cols[i],
      test                 = fit$test,
      crit_val             = fit$crit.val,
      p_val                = fit$p.value
    )
  })
  
  binded_res <- do.call(rbind, res)
  
  binded_res$p_val_adj <- 
    p.adjust(binded_res$p_val, method = "BH")
  
  attr(binded_res, "call") <- matched_call
  attr(binded_res, "formula") <- model_formula
  
  binded_res
}

median_aov_parallel <- function(
    data,
    cols,
    class_col,
    nboot   = 1000,
    n_cores = parallel::detectCores() - 1
) {
  # ---- Record the call and formula ----
  
  matched_call  <- match.call()
  model_formula <- stats::as.formula(paste("y ~", class_col))
  message("Model formula: ", deparse(model_formula))
  
  # ---- Register parallel backend ----
  
  n_cores <- max(1L, as.integer(n_cores))
  cl      <- parallel::makeCluster(n_cores)
  doParallel::registerDoParallel(cl)
  
  # Ensure the cluster is always torn down on function exit,
  # even if an error is thrown, and reset to sequential mode.
  on.exit(
    {
      parallel::stopCluster(cl)
      foreach::registerDoSEQ()
    },
    add = TRUE
  )
  
  message(sprintf("Running on %d core(s).", n_cores))
  
  # ---- Fit models ----
  
  mat   <- as.matrix(data[, cols])
  class <- data[[class_col]]
  
  # Suppress R CMD check note: "no visible binding for variable 'i'"
  i <- NULL
  
  res <- foreach::foreach(
    i         = seq_along(cols),
    .packages = "WRS2",  # loaded on every worker
    .combine  = rbind    # replaces do.call(rbind, res)
  ) %dopar% {
    
    y <- mat[, i]
    
    fit <- tryCatch(
      WRS2::med1way(y ~ class, iter = nboot),
      error = function(e) {
        warning(
          sprintf(
            "med1way failed for '%s': %s",
            cols[i], conditionMessage(e)
          ),
          call. = FALSE
        )
        NULL
      }
    )
    
    if (is.null(fit)) {
      data.frame(
        variable = cols[i],
        test     = NA_real_,
        crit_val = NA_real_,
        p_val    = NA_real_
      )
    } else {
      data.frame(
        variable = cols[i],
        test     = fit$test,
        crit_val = fit$crit.val,
        p_val    = fit$p.value
      )
    }
  }
  
  rownames(res) <- NULL
  
  res$p_val_adj <- p.adjust(res$p_val, method = "BH")
  
  attr(res, "call")    <- matched_call
  attr(res, "formula") <- model_formula
  
  res
}
