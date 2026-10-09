# Calculates the relative standard deviation of a given vector
calculate_RSD <- function(data) {
  n <- length(data)
  pop_sd <- sqrt((n-1)/n) * sd(data)
  rsd <- pop_sd/mean(data)*100
  rsd
}

# Helper function, changes the cols variable to only include those present in
# the names of the dataframe.
cleanup_cols <- function(data, cols) {
  
  # Check cols is a character vector
  if (!is.character(cols)) {
    stop("ERROR: `cols` should be a character vector with the names of the target columns")
  }
  
  # Get the names of the dataframe, to save evaluating it twice
  nm <- names(data)
  
  # Subset nm for anything in cols
  result <- nm[nm %in% cols]
  result
}

# Takes the names from each of two data columns that match the given pattern. 
# Then checks that the names in both match one another (although not the order)
# If names do match, it outputs a vector of the names, otherwise it throws an error
get_metabolite_columns <- function(neg_mode, pos_mode, pattern) {
  # Get each set of names
  data_cols_neg <- names(neg_mode)[grep(pattern, names(neg_mode))]
  data_cols_pos <- names(pos_mode)[grep(pattern, names(pos_mode))]
  
  # Logic to check that the vectors contain the same names
  if(sum(!data_cols_neg %in% data_cols_pos) > 0 ||
     sum(!data_cols_pos %in% data_cols_neg) > 0) {
    stop("names for positive and negative modes do not match")
  }
  # Return a vector of names
  data_cols_neg
}

# Replaces NA values in the "Names" column with "Unknown_1", "Unknown_2", 
# "Unknown_3"...
NA_2_Unknown <- function(.data, name_col = Name, new_names = "Unknown_") {
  na_free_data <- .data %>% 
    mutate({{name_col}} := if_else(
      is.na({{name_col}}),
      paste0(new_names, cumsum(is.na({{name_col}}))),
      as.character({{name_col}})
    ))
  
  na_free_data
}

# Removes metabolites with a low ID score
remove_low_scores <- function(data, scores, threshold = 80) {
  # Work out which rows have scores above the threshold
  logic <- scores >= threshold
  
  if (sum(!logic) > 0) {
    # If there is at least one score below the threshold, subset the dataset
    result <- data[logic, ]
  } else {
    # If there are no scores below the threshold, simply return the dataset
    result <- data
  }
  result
}

# Removes metabolites with a mass error above the threshold
remove_high_mass_error <- function(data, masserror, threshold = 5) {
  # Work out which rows have ME below the threshold
  logic <- masserror <= threshold
  
  # If there is at least one score below the threshold, subset the dataset
  if (sum(!logic) > 0) {
    result <- data[logic, ]
  # If there are no scores below the threshold, simply return the dataset
  } else {
    result <- data
  }
  result
}

# Calculates the baseline for a dataset by calculating row means for any columns
# with `pattern` in the name
calculate_baseline <- function(data, pattern = "Recon") {
  cols <- grep(pattern, names(data))   # indices of matching columns
  if(length(cols) == 0) {
    warning("Pattern does not match any columns. Returning a baseline vector of
            only zeroes")
    baseline <- c(rep(0, nrow(data)))
  } else {
    baseline <- rowMeans(data[, cols, drop = FALSE], na.rm = TRUE) 
  }
  baseline
}

# Removes metabolites from analysis that have excessive detection in blanks
# Uses a pattern in dataframe names to detect blanks ("Proc" by default).
# The condition for removing a metabolite is if the blanks' mean is greater than
# or equal to a fraction of the mean of the actual samples. 1/5 mean is default.

remove_noisy_baseline <- function(
    data, cols,
    pattern = "Proc",
    threshold = 0.2) {
  # Removes any column names that are not present in names(data)
  clean <- cleanup_cols(data = data, cols = cols)
  # Find names of baseline columns
  baseline_cols <- names(data)[str_detect(names(data), pattern)]
  # Error checking
  if(length(baseline_cols) < 1) {
    warning("`pattern` does not match any columns. Returning `data` unchanged.")
    return(data)
  }
  # Isolate the data columns as their own dataframe
  metabolites <- data %>% 
    select(any_of(clean)) %>% 
    # Calculate the mean of each row multiplied by the threshold
    rowMeans(na.rm = TRUE) * threshold
  # Isolate the baseline columns as their own dataframe
  baseline <- data %>% 
    select(all_of(baseline_cols)) %>% 
    # Calculate the mean of each row
    replace(is.na(.), 0) %>% 
    rowMeans()
  # Work out which metabolites to keep
  logic <- !baseline > metabolites
  # Filter the data based on `logic`
  result <- data[logic, ]
} 

# Subtracts the baseline from each column
# Baseline should be a named vector with names equivalent to the columns you
# want to correct
subtract_baseline <- function(data, cols, baseline_col = "baseline") {
  # Removes any column names that are not present in names(data)
  clean <- cleanup_cols(data = data, cols = cols)
  
  # Get the baseline column as a vector to prevent having to subset over and over
  base_vec <- data[[baseline_col]]
  
  # Substract the baseline from each column
  for (i in clean) {
    data[[i]] <- data[[i]] - base_vec
  }
  
  # Return the modified data without the baseline column
  data[[baseline_col]] <- NULL
  data
}

# Normalises each column to the weight
# `weights` should be a named vector with names equivalent to the columns you
# want to normalise
weight_normalise <- function(data, cols, weights) {
  # Removes any column names that are not present in names(data)
  clean <- cleanup_cols(data = data, cols = cols)
  
  # Multiply each column by the weight
  for (i in clean) {
    data[[i]] <- data[[i]] * weights[[i]]
  }
  
  # Return the modified data without the baseline column
  data
}

# QToF analyses detect the same metabolites multiple times, sometimes in both methods.
# This funtion removes duplicates by summing. It also updates the "Method" column
# to indicate whether a metabolite was detected only in POS mode, only in NEG,
# or in both. 

sum_same_metabolites <- 
  function(data, cols, metabolite_names = "Name", mode_col = "Method") {
    
    # Split numeric and non-numeric aggregation targets
    agg_numeric <- aggregate(
      data[cols],
      by = list(Name = data[[metabolite_names]]),
      FUN = sum,
      na.rm = TRUE
    )
    
    agg_mode <- aggregate(
      data[[mode_col]],
      by = list(Name = data[[metabolite_names]]),
      FUN = function(x) paste(sort(unique(x)), collapse = "/")
    )
    
    # Merge results
    out <- merge(agg_mode, agg_numeric, by = "Name")
    
    names(out)[names(out) == "x"] <- mode_col
    out
  }

# Iterates through specified columns and replaces any values below the threshold
# with `NA`
low_count_filter <- function(data, cols = NULL, threshold = 1500) {
  # If no columns are specified, cols == all columns
  if (is.null(cols)) {
    cols <- names(data)
  }
  
  # Removes any column names that are not present in names(data)
  clean <- cleanup_cols(data = data, cols = cols)
  
  # Sets all values that are below the threshold to `NA`
  for(i in clean) {
    data[,i][data[,i] < threshold] <- NA_real_
  }
  data
}

# Transposes the QToF data into wide format
transpose_selected_cols <- function(data, cols, names = "Name") {
  # Get the column specified by "Name" as a vector to supply as a column names later
  nm <- data[[names]]
  
  # Check values in nm are unique. Throw an error if not
  if(! length(nm) == length(unique(nm))) {
    stop("ERROR: names are not unique!")
  }
  
  # Removes any column names that are not present in names(data)
  clean <- cleanup_cols(data = data, cols = cols)
  
  # Transpose dataframe. This converts it to a matrix, so convert it back to
  # a dataframe, then use nm to name the columns.
  result <- t(data[clean])
  result <- as.data.frame(result)
  names(result) <- nm
  
  result
}

# Merges the metadata and data frames provided (in wide format).
# Uses the rownames of both to ensure that they match up.
# If rownames are not the same in both dataframes (ignoring order), then the 
# function throws an error
merge_metadata <- function(data, metadata) {
  # Get rownames into vectors
  data_nm <- rownames(data)
  metadata_nm <- rownames(metadata)
  
  # Logic to check that the vectors contain the same names
  if(sum(!data_nm %in% metadata_nm) > 0 ||
     sum(!metadata_nm %in% data_nm) > 0) {
    stop("Row names for data and metadata do not match")
  }
  
  # Merge the dataframes by rownames, using the order specified in the metadata
  result <- cbind(metadata, data[metadata_nm, ])
  
  # Remove row names
  rownames(result) <- NULL
  
  # Return result
  result
}

# Groups must be a named vector of the same length as cols
eighty_per_cent_rule <- function(data, cols, groups = NULL, threshold = 0.8) {
  # Create a subsetted dataframe with only the data columns
  no_metadata <- data[cols]
  
  if (is.null(groups)) {
    # Ungrouped logic
    # Create a logical vector indicating which rows have more than 80% of values NOT NA
    keep <- rowMeans(!is.na(no_metadata)) >= threshold
    # Subest the original dataframe and save the result
    result <- data[keep, ]
    
  } else if (!is.null(groups) && 
             length(groups) == length(cols) &&
             length(names(groups)) != 0 &&
             is.character(groups)) {
    # Grouped logic
    
    if (!all(names(groups) %in% names(data))) {
      stop("All class names must correspond to data columns")
    }
    
    group_levels <- unique(groups)
    
    keep <- vapply(
      seq_len(nrow(no_metadata)),
      function(i) {
        row <- no_metadata[i, ]
        any(vapply(
          group_levels,
          function(gr) {
            subset_one_group <- names(groups)[groups == gr]
            mean(!is.na(row[subset_one_group])) >= threshold
          },
          logical(1)
        ))
      },
      logical(1)
    )
    
    result <- data[keep, ]
    
  } else {
    stop("`groups` must be either `NULL` or a named vector with the same length 
        as cols")
  }
  result
}

# Where there are identical metabolites in either/both positive and negative
# mode, this function chooses the one with the highest score

POS_NEG_Merge <- function(pos_mode_data, neg_mode_data) {
  merged <- bind_rows(pos_mode_data, neg_mode_data, .id = "Method") %>% 
    mutate(Method = as.character(Method),
           Method = str_replace_all(Method, "1", "QToF_POS"),
           Method = str_replace_all(Method, "2", "QToF_NEG"))
  merged
}

# Create function to filter low values out of treatment columns. Low values are 
# specified by the threshold. These are appropriate to the device used for 
# metabolomics. Basically exists to remove noise.
low_filter <- function(data, threshold = 1500) {
  filtered <- ifelse(data < threshold,
                     NA,
                     data)
  
  filtered
}
# Create a function to impute low values to 20% (or any number you choose) of
# the median.
impute_lows <- function(data, impute_multiplier = 0.2) {
  imputed_vals <- median(data, na.rm = TRUE) * impute_multiplier
  
  
  imputed <- ifelse(is.na(data),
                    imputed_vals,
                    data)
  imputed
}

remove_duplicate_metabolites <- function(data, names = "Name") {
  # Display a warning to ensure users have already arranged merged data in
  # the hierarchy of most reliable
  warning("Metabolomics data must be in a single dataframe with data arranged in 
          a hierarchy of most to least confident")
  # Create a logical vector of duplicates in the names vector.
  dups <- duplicated(data[names])
  # Remove rows from `data` that are TRUE in logical test
  result <- data[!dups,]
}

# Accessory function. Applies a scaler specified using `function(x) {...}` notation
apply_scaler <- function(data, cols = NULL, scaler) {
  if(is.null(cols)) {
    # If no columns are specified, use the names of `data
    clean <- names(data[ , unlist(lapply(data, is.numeric))])
  } else {
    # Remove elements of the cols vector that aren't in names(data)
    clean <- cleanup_cols(data = data, cols = cols)
  }
  
  # Find the metadata cols (all those not in `cols`)
  metadata <- names(data)[! names(data) %in% clean]
  
  # Apply a scaler to the data columns
  # scaled <- sapply(data[, clean, drop = FALSE],
  #                  function (x) scaler(x))
  
  scaled <- t(apply(
    data[, clean, drop = FALSE],
    MARGIN = 1,
    FUN = scaler
  ))
  # Bind metadata and modified data
  result <- cbind(data[metadata], as.data.frame(scaled))
  
}

# Accessory function. Returns NA instead of NaN for log10. Also suppresses 
# warnings to avoid now-redundant `NaN` warnings.
log10_na <- function(x) {
  y <- suppressWarnings(log10(x))
  y[is.nan(y)] <- NA
  y
}

# Log10 scales numeric columns, or the columns provided to cols
log10_scale <- function(data, cols = NULL) {
  apply_scaler(
    data,
    cols,
    scaler = function(x) {log10_na(x)}
  )
}

# Accessory function. Returns NA instead of NaN for natural log. Also suppresses 
# warnings to avoid now-redundant `NaN` warnings.
ln_na <- function(x) {
  y <- suppressWarnings(log(x))
  y[is.nan(y)] <- NA
  y
}

# Natural log scales numeric columns, or the columns provided to cols
ln_scale <- function(data, cols = NULL) {
  apply_scaler(
    data,
    cols,
    scaler = function(x) {ln_na(x)}
  )
}

# Median normalises numeric columns, or the columns provided to cols
median_normalise <- function(data, cols = NULL) {
  apply_scaler(
    data,
    cols,
    scaler = function(x) {x / median(x, na.rm = TRUE)}
  )
}

# Pareto scales numeric columns, or the columns provided to cols
pareto_scale <- function(data, cols = NULL) {
  apply_scaler(
    data,
    cols,
    scaler = function(x) {
      (x - mean(x, na.rm = TRUE)) / sqrt(sd(x, na.rm = TRUE))
    }
  )
}

# Performs a Quantile Regression Imputation of Left-Consored data while retaining
# metadata. Requires pkg `imputeLCMD`
# • `data` is a dataframe with some columns to be imputed
# • `cols` are the columns to be imputed. All otehr cols will be returned 
#    unchanged
QRILC_impute_with_metadata <- function(data, cols = NULL) {
  if (!requireNamespace("imputeLCMD", quietly = TRUE)) {
    stop("Package 'imputeLCMD' is required but not installed.")
  }
  
  if(is.null(cols)) {
    # If no columns are specified, use the names of `data
    clean <- names(data[ , unlist(lapply(data, is.numeric))])
  } else {
    # Remove elements of the cols vector that aren't in names(data)
    clean <- cleanup_cols(data = data, cols = cols)
  }
  
  # Find the metadata (should be any data whose names are not in `cols`)
  metadata_nm <- names(data)[! names(data) %in% clean]
  metadata <- data[metadata_nm]
  
  # Impute missing data using Quantile Regression Imputation of Left-Censored
  # data on the data columns only, then bind them to the metadata
  imputed <- imputeLCMD::impute.QRILC(data[, clean])[[1]]
  
  # Bind the metadata and the imputed data
  result <- cbind(metadata, imputed)
  
}

# Helper function for WRS2 ANOVA-style test wrappers (to follow)


.validate_wrs2_inputs <- function(
    data,
    cols,
    single_cols = list(),
    factor_cols = character(0),
    min_levels  = 2L
) {
  # --- data.frame check ---
  if (!is.data.frame(data)) {
    stop("`data` must be a data.frame.", call. = FALSE)
  }
  
  # --- WRS2 availability ---
  if (!requireNamespace("WRS2", quietly = TRUE)) {
    stop(
      "Package 'WRS2' is required but not installed.",
      call. = FALSE
    )
  }
  
  # --- metabolite columns present in data ---
  missing_cols <- setdiff(cols, colnames(data))
  if (length(missing_cols) > 0) {
    stop(
      "The following `cols` are not found in `data`: ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }
  
  # --- single-column arguments present in data ---
  # single_cols is a named list: list(class_col = "treatment", ...)
  for (arg_name in names(single_cols)) {
    arg_val <- single_cols[[arg_name]]
    if (length(arg_val) != 1L || !arg_val %in% colnames(data)) {
      stop(
        sprintf("`%s` ('%s') not found in `data`.", arg_name, arg_val),
        call. = FALSE
      )
    }
  }
  
  # --- coerce factor columns ---
  # factor_cols is a character vector of keys from single_cols
  for (arg_name in factor_cols) {
    col_name <- single_cols[[arg_name]]
    if (!is.factor(data[[col_name]])) {
      message(sprintf("Coercing '%s' to factor.", col_name))
      data[[col_name]] <- as.factor(data[[col_name]])
    }
  }
  
  # --- metabolite columns numeric ---
  non_numeric <- cols[!vapply(data[cols], is.numeric, logical(1))]
  if (length(non_numeric) > 0) {
    stop(
      "The following `cols` are not numeric: ",
      paste(non_numeric, collapse = ", "),
      call. = FALSE
    )
  }
  
  # --- minimum levels for each factor column ---
  # min_levels can be a scalar (applied to all) or a named integer vector
  # e.g. c(class_col = 2L, time_col = 3L)
  for (arg_name in factor_cols) {
    col_name <- single_cols[[arg_name]]
    min_lvl <- if (
      !is.null(names(min_levels)) && arg_name %in% names(min_levels)
    ) {
      min_levels[[arg_name]]
    } else {
      min_levels[[1L]]
    }
    if (nlevels(data[[col_name]]) < min_lvl) {
      stop(
        sprintf(
          "`%s` must have at least %d level(s).",
          arg_name, min_lvl
        ),
        call. = FALSE
      )
    }
  }
  
  data
}

# Uses Wilcox's trimmed mean with repeated measures to find significantly varying
# variables. Outputs all the results in a simple dataframe and also calculates
# Benjamini-Hochberg adjusted p-values.
#
# data      - a data.frame containing all relevant columns
# cols      - indicates cols where metabolite data are stored
# class_col - the column with class information (i.e. Treatment groups)
# time_col  - the column with time variable information
# id_col    - indicates the experimental units. There must be a unique value for
#             each experimental unit (i.e. 3 treatments, 10 reps == 30 unique values)
trimmed_mean_1way_aov_bt_parallel <- function(
  data,
  cols,
  class_col,
  tr = 0.2,
  nboot   = 599,
  n_cores = parallel::detectCores() - 1
) {
  
  # --- Error checking ---
  
  data <- .validate_wrs2_inputs(
    data        = data,
    cols        = cols,
    single_cols = list(class_col = class_col),
    factor_cols = c("class_col"),
    min_levels  = 2L            # scalar: both must have >= 2 levels
  )
  
  # --- Record the call and formula ---
  
  matched_call  <- match.call()
  model_formula <- stats::as.formula(paste("y ~", class_col))
  message("Model formula: ", deparse(model_formula))
  
  # --- Register parallel backend ---
  
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
  
  # --- Fit models ---
  
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
      WRS2::t1waybt(y ~ class,
                    tr = tr,
                    nboot = nboot
                  ),
      error = function(e) {
        warning(
          sprintf(
            "t1waybt failed for '%s': %s",
            cols[i], conditionMessage(e)
          ),
          call. = FALSE
        )
        NULL
      }
    )
    
    if (is.null(fit)) {
      data.frame(
        variable      = cols[i],
        test          = NA_real_,
        Var.Explained = NA_real_,
        Effect.Size   = NA_real_,
        p_val         = NA_real_
      )
    } else {
      data.frame(
        variable      = cols[i],
        test          = fit$test,
        Var.Explained = fit$Var.Explained,
        Effect.Size   = fit$Effect.Size,
        p_val         = fit$p.value
      )
    }
  }
  
  rownames(res) <- NULL
  
  res$p_val_adj <- p.adjust(res$p_val, method = "BH")
  
  attr(res, "call")    <- matched_call
  attr(res, "formula") <- model_formula
  
  res
}

trimmed_mean_1way_aov <- function(data, cols, class_col, tr = 0.2) {
  
  # --- Error checking ---
  
  data <- .validate_wrs2_inputs(
    data        = data,
    cols        = cols,
    single_cols = list(class_col = class_col),
    factor_cols = c("class_col"),
    min_levels  = 2L            # scalar: both must have >= 2 levels
  )
  
  # --- Record the call and formula ---
  
  matched_call <- match.call()
  model_formula <- stats::as.formula(
    paste("y ~", class_col)
  )
  message("Model formula: ", deparse(model_formula))
  
  # --- Fit models ---
  
  mat   <- as.matrix(data[, cols])
  class <- data[[class_col]]
  
  res <- lapply(seq_along(cols), function(i) {
    
    y <- mat[, i]
    
    fit <- tryCatch(
      WRS2::t1way(
        y ~ class,
        id = id,
        tr = 0.2
      ),
      error = function(e) {
        warning(
          sprintf(
            "t1way failed for '%s': %s", cols[i], conditionMessage(e)
          ),
          call. = FALSE
        )
        return(NULL)
      }
    )
    
    if (is.null(fit)) {
      return(
        data.frame(
          variable       = cols[i],
          test           = NA_real_,
          df1            = NA_real_,
          df2            = NA_real_,
          p_val          = NA_real_,
          effsize        = NA_real_,
          effsize_ci_low = NA_real_,
          effsize_ci_hi  = NA_real_
        )
      )
    }
    
    data.frame(
      variable       = cols[i],
      test           = fit$test,
      df1            = fit$df1,
      df2            = fit$df2,
      p_val          = fit$p.value,
      effsize        = fit$effsize,
      effsize_ci_low = fit$effsize_ci[[1]],
      effsize_ci_hi  = fit$effsize_ci[[2]]
    )
  })
  
  binded_res <- do.call(rbind, res)
  
  # --- BH-adjusted p-values ---
  
  binded_res$p_val_adj <-
    p.adjust(
      binded_res$p_val,
      method = "BH")
  
  # Hard-coded reorder: place each _adj column right after its _p_val column
  binded_res <- binded_res[, c(
    "variable",
    "test",
    "df1",
    "df2",
    "p_val",
    "p_val_adj",
    "effsize",
    "effsize_ci_low",
    "effsize_ci_hi"
  )]
  
  # --- Attach metadata as attributes ---
  
  attr(binded_res, "call")    <- matched_call
  attr(binded_res, "formula") <- model_formula
  attr(binded_res, "tr")      <- 0.2
  
  binded_res
}

trimmed_mean_rep_meas_aov <- function(data, cols, class_col, time_col, id_col) {
  
  # --- Error checking ---
  
  data <- .validate_wrs2_inputs(
    data        = data,
    cols        = cols,
    single_cols = list(class_col = class_col),
    factor_cols = c("class_col"),
    min_levels  = 2L            # scalar: both must have >= 2 levels
  )
  
  # --- Record the call and formula ---
  
  matched_call <- match.call()
  model_formula <- stats::as.formula(
    paste("y ~", class_col, "*", time_col)
  )
  message("Model formula: ", deparse(model_formula))
  
  # --- Fit models ---
  
  mat   <- as.matrix(data[, cols])
  class <- data[[class_col]]
  time  <- data[[time_col]]
  id    <- data[[id_col]]
  
  res <- lapply(seq_along(cols), function(i) {
    
    y <- mat[, i]
    
    fit <- tryCatch(
      WRS2::bwtrim(
        y ~ class * time,
        id = id,
        tr = 0.2
      ),
      error = function(e) {
        warning(
          sprintf(
            "bwtrim failed for '%s': %s", cols[i], conditionMessage(e)
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
          class_Q              = NA_real_,
          class_df1            = NA_real_,
          class_df2            = NA_real_,
          class_p_val          = NA_real_,
          time_Q               = NA_real_,
          time_df1             = NA_real_,
          time_df2             = NA_real_,
          time_p_val           = NA_real_,
          interaction_Q        = NA_real_,
          interaction_df1      = NA_real_,
          interaction_df2      = NA_real_,
          interaction_p_val    = NA_real_
        )
      )
    }
    
    data.frame(
      variable            = cols[i],
      class_Q             = fit$Qa,
      class_df1           = fit$A.df[1],
      class_df2           = fit$A.df[2],
      class_p_val         = fit$A.p.value,
      time_Q              = fit$Qb,
      time_df1            = fit$B.df[1],
      time_df2            = fit$B.df[2],
      time_p_val          = fit$B.p.value,
      interaction_Q       = fit$Qab,
      interaction_df1     = fit$AB.df[1],
      interaction_df2     = fit$AB.df[2],
      interaction_p_val   = fit$AB.p.value
    )
  })
  
  binded_res <- do.call(rbind, res)
  
  # --- BH-adjusted p-values ---
  
  binded_res$class_p_val_adj       <- p.adjust(binded_res$class_p_val,       method = "BH")
  binded_res$time_p_val_adj        <- p.adjust(binded_res$time_p_val,        method = "BH")
  binded_res$interaction_p_val_adj <- p.adjust(binded_res$interaction_p_val, method = "BH")
  
  # Hard-coded reorder: place each _adj column right after its _p_val column
  binded_res <- binded_res[, c(
    "variable",
    "class_Q", "class_df1", "class_df2",
    "class_p_val", "class_p_val_adj",
    "time_Q", "time_df1", "time_df2",
    "time_p_val", "time_p_val_adj",
    "interaction_Q", "interaction_df1", "interaction_df2",
    "interaction_p_val", "interaction_p_val_adj"
  )]
  
  # --- Attach metadata as attributes ---
  
  attr(binded_res, "call")    <- matched_call
  attr(binded_res, "formula") <- model_formula
  attr(binded_res, "tr")      <- 0.2
  
  binded_res
}

median_1way_aov <- function(data, cols, class_col, nboot = 1000) {
  # --- Error checking ---
  
  data <- .validate_wrs2_inputs(
    data        = data,
    cols        = cols,
    single_cols = list(class_col = class_col),
    factor_cols = c("class_col"),
    min_levels  = 2L            # scalar: both must have >= 2 levels
  )
  
  # --- Record the call and formula ---
  
  matched_call <- match.call()
  model_formula <- stats::as.formula(
    paste("y ~", class_col)
  )
  message("Model formula: ", deparse(model_formula))
  
  # --- Fit models ---
  
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

median_1way_aov_parallel <- function(
    data,
    cols,
    class_col,
    nboot   = 1000,
    n_cores = parallel::detectCores() - 1
) {
  # --- Error checking ---
  
  data <- .validate_wrs2_inputs(
    data        = data,
    cols        = cols,
    single_cols = list(class_col = class_col),
    factor_cols = c("class_col"),
    min_levels  = 2L            # scalar: both must have >= 2 levels
  )
  
  # --- Record the call and formula ---
  
  matched_call  <- match.call()
  model_formula <- stats::as.formula(paste("y ~", class_col))
  message("Model formula: ", deparse(model_formula))
  
  # --- Register parallel backend ---
  
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
  
  # --- Fit models ---
  
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

yuen_trimmed_t_test <- function(data, cols, class_col, tr = 0.2) {
  
  # --- Error checking ---
  
  data <- .validate_wrs2_inputs(
    data        = data,
    cols        = cols,
    single_cols = list(class_col = class_col),
    factor_cols = c("class_col"),
    min_levels  = 2L            # scalar: both must have >= 2 levels
  )
  
  # Check number of levels in class class_col is two
  
  if(length(levels(data[[class_col]])) != 2) {
    stop("`class_col` doesn't have two levels.")
  }
  
  # --- Record the call and formula ---
  
  matched_call <- match.call()
  model_formula <- stats::as.formula(
    paste("y ~", class_col)
  )
  message("Model formula: ", deparse(model_formula))
  
  # --- Fit models ---
  
  mat   <- as.matrix(data[, cols])
  class <- data[[class_col]]
  
  res <- lapply(seq_along(cols), function(i) {
    
    y <- mat[, i]
    
    fit <- tryCatch(
      WRS2::yuen(
        y ~ class,
        tr = tr
      ),
      error = function(e) {
        warning(
          sprintf(
            "yuen failed for '%s': %s", cols[i], conditionMessage(e)
          ),
          call. = FALSE
        )
        return(NULL)
      }
    )
    
    if (is.null(fit)) {
      return(
        data.frame(
          variable       = cols[i],
          test           = NA_real_,
          df             = NA_real_,
          difference     = NA_real_,
          p_val          = NA_real_,
          effsize        = NA_real_,
          ci_low         = NA_real_,
          ci_hi          = NA_real_
        )
      )
    }
    
    data.frame(
      variable       = cols[i],
      test           = fit$test,
      df             = fit$df,
      difference     = fit$diff,
      p_val          = fit$p.value,
      effsize        = fit$effsize,
      ci_low         = fit$conf.int[[1]],
      ci_hi          = fit$conf.int[[2]]
    )
  })
  
  binded_res <- do.call(rbind, res)
  
  
  
  binded_res$p_val_adj <- 
    p.adjust(binded_res$p_val, method = "BH")
  
  attr(binded_res, "call") <- matched_call
  attr(binded_res, "formula") <- model_formula
  
  binded_res <- binded_res[ , c(
    "variable",
    "test",
    "df",
    "difference",
    "effsize",
    "ci_low",
    "ci_hi",
    "p_val",
    "p_val_adj"
  )]
}  

one_way_mcppb20 <- function(metabolites, data) {
  
  # ERROR CHECKING
  if (!is.character(metabolites)) {
    stop("`metabolites` must be a character vector.")
  }
  
  if (length(metabolites) == 0) {
    stop("`metabolites` contains no metabolites.")
  }
  
  if (!is.data.frame(data)) {
    stop("`data` must be a data frame.")
  }
  
  if (!"Treatment" %in% names(data)) {
    stop("`data` must contain a column called `Treatment`.")
  }
  
  missing_metabolites <- setdiff(metabolites, names(data))
  
  if (length(missing_metabolites) > 0) {
    stop("The following metabolites are not columns in `data`: ",
         paste(missing_metabolites, collapse = ", "))
  }
  
  # PRE-BUILD LISTS FOR RESULTS
  post_hocs <- list()
  clds <- list()
  
  for (metabolite in metabolites) {
    # PERFORM % BOOTSTRAPPED MULTIPLE COMPARISONS TEST
    # Create a formula to pass to our multiple comparisons test
    formula <- as.formula(paste0("`", metabolite, "`~Treatment"))
    
    # Calculate robust multiple comparisons
    model <- WRS2::mcppb20(formula, data = data)
    
    # Isolate the group names. It is computationally slower to do this 
    # repeatedly but I'm doing it so the code is more robust.
    fnames <- model$fnames
    
    # Isolate the actual test results and put them into a workable dataframe
    post_hoc_result <- as.data.frame(model$comp)
    names(post_hoc_result)[1:2] <- c("Comparison_Group_1", 
                                     "Comparison_Group_2")
    
    # Replace indices in "Group" columns with actual names
    for (column in 1:2) {
      for (row in 1:nrow(post_hoc_result)) {
        post_hoc_result[row,column] <- 
          fnames[as.numeric(post_hoc_result[row,column])]
      }
    }
    
    # Merge the "Group" columns into meaningful labels
    post_hoc_result$Comparison <- paste0(post_hoc_result$Comparison_Group_1,
                                         "-", 
                                         post_hoc_result$Comparison_Group_2)
    
    # Remove extraneous columns
    post_hoc_result[1:2] <- NULL
    
    # GENERATE COMPACT LETTER DISPLAYS
    # Get p values in a `multcompLetters()`-friendly format (named num vector)
    p.vals <- post_hoc_result$`p-value`
    names(p.vals) <- post_hoc_result$Comparison
    
    # Create the CLD
    cld_result <- multcompLetters(p.vals,
                                  threshold = 0.05)
    
    # WRITE THE RESULTS INTO THE PRE-BUILT LISTS
    post_hocs[[metabolite]] <- post_hoc_result
    clds[[metabolite]] <- cld_result
  }
  
  # RETURN THE RESULTS
  return(list(
    post_hocs = post_hocs,
    clds = clds
  ))
}

print("Metabolomics_Funcs.R ran successfully")