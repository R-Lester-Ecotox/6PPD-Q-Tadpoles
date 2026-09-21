
# ---- Helpers ---------------------------------------------------------

#' Report on the reliability of an NMDS stress value
check_stress <- function(stress) {
  msg <- sprintf("NMDS stress = %.4f", stress)
  if (stress >= 0.20) {
    warning(msg, " - above the unreliability threshold of 0.20.",
            call. = FALSE)
  } else if (stress >= 0.10) {
    message(msg, " - okay. Broad interpretation only.")
  } else if (stress >= 0.05) {
    message(msg, " - good. Safe to interpret clusters.")
  } else {
    message(msg, " - extremely low! Layout is basically exact.")
  }
  invisible(stress)
}

#' PERMANOVA + dispersion check + NMDS ordination plot for one grouping
#' variable.
#'
#' @param group_var Column name in `data` (character, backticks not needed).
#' @param data Sample metadata, one row per sample, in the same order as `d`.
#' @param d A `dist` object of between-sample distances.
#' @param nmds_fit A fitted `metaMDS` object built from the same samples.
#' @param perm_count Number of permutations.
#' @param out_dir Directory for the CSV and PNG outputs.
#' @param plot_title Optional plot title; defaults to "<group_var> ...".
#' @param file_prefix Optional filename stem; defaults to a sanitised
#'   `group_var`.
#' @return Invisibly, a list with the PERMANOVA table, dispersion test,
#'   variance explained and output paths.
permanova_by_group <- function(group_var,
                               data,
                               d,
                               nmds_fit,
                               perm_count = 9999,
                               out_dir     = ".",
                               plot_title  = NULL,
                               file_prefix = NULL,
                               xlim   = NULL,
                               ylim   = NULL,
                               width  = 1800,
                               height = 1800,
                               dpi    = 300,
                               palette = NULL) {
  
  stopifnot(is.character(group_var), length(group_var) == 1L)
  if (!group_var %in% names(data)) {
    stop("`", group_var, "` is not a column in `data`.", call. = FALSE)
  }
  if (attr(d, "Size") != nrow(data)) {
    stop("`d` has ", attr(d, "Size"), " samples but `data` has ",
         nrow(data), " rows.", call. = FALSE)
  }
  
  if (is.null(file_prefix)) {
    file_prefix <- gsub("^_|_$", "", gsub("[^A-Za-z0-9]+", "_", group_var))
  }
  if (is.null(plot_title)) {
    plot_title <- paste(group_var, "NMDS Ordination")
  }
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  
  groups <- factor(data[[group_var]])
  keep   <- !is.na(groups)
  if (any(!keep)) {
    warning(sum(!keep), " sample(s) dropped: NA in `", group_var, "`.",
            call. = FALSE)
  }
  groups <- droplevels(groups[keep])
  if (nlevels(groups) < 2L) {
    stop("`", group_var, "` needs at least two non-empty levels.",
         call. = FALSE)
  }
  
  # Subset the distance matrix if any samples were dropped
  d_sub <- if (all(keep)) {d} else {stats::as.dist(as.matrix(d)[keep, keep])}
  
  # -- PERMANOVA -------------------------------------------------------
  fit_data <- data.frame(grp = groups)
  permanova_model <- adonis2(d_sub ~ grp,
                             data         = fit_data,
                             permutations = perm_count,
                             #by           = "terms"
                             )
  
  res <- as.data.frame(permanova_model)
  rownames(res)[1] <- group_var
  names(res)[grep("^Pr\\(", names(res))] <- "model_p_val"
  
  var_explained <- res$R2[1] * 100
  
  # -- Dispersion check ------------------------------------------------
  bd <- permutest(betadisper(d_sub, groups), permutations = perm_count)
  
  res$dispersion_p_val <- c(bd$tab[["Pr(>F)"]][1],
                            rep(NA_real_, nrow(res) - 1L))
  
  notes <- c(
    paste0("R2 is the proportion (0-1) of variance in the distance ",
           "matrix accounted for by `", group_var, "`."),
    "model_p_val < 0.05 indicates a significant effect of the variable.",
    paste0("dispersion_p_val < 0.05 means within-group dispersions ",
           "DIFFER, i.e. the homogeneity assumption is violated and ",
           "apparent separation may be a spread artefact. p >= 0.05 ",
           "supports the assumption. Less critical in balanced designs.")
  )
  res$notes <- c(notes, rep("", max(0L, nrow(res) - length(notes))))[
    seq_len(nrow(res))]
  
  csv_path <- file.path(out_dir, paste0(file_prefix, "_PERMANOVA.csv"))
  write.csv(res, file = csv_path)
  
  # -- Ordination plot -------------------------------------------------
  scrs <- scores(nmds_fit, display = "sites")[keep, , drop = FALSE]
  
  plt <- if (is.null(palette)) {
    grDevices::hcl.colors(nlevels(groups), palette = "Dark 3")
  } else {
    rep_len(palette, nlevels(groups))
  }

  png_path <- file.path(out_dir,
                        paste0(file_prefix, "_PERMANOVA_Ordination.png"))
  
  CairoPNG(png_path, width = width, height = height, dpi = dpi)
  on.exit(if (dev.cur() > 1L) dev.off(), add = TRUE)
  
  par(
    cex      = 1.2,
    cex.axis = 1.0,
    cex.lab  = 1.4,
    cex.main = 1.5,
    cex.sub  = 1.2
  )
  
  plot(scrs, type = "n", asp = 1,
       xlim = xlim, ylim = ylim,
       xlab = "NMDS1", ylab = "NMDS2",
       main = plot_title,
       sub  = sprintf(
         "Stress = %.3f  |  R2 = %.1f%%  |  p = %.4f",
         nmds_fit$stress, var_explained, res$model_p_val[1]))
  points(scrs, col = plt[groups], pch = 19)
  ordiellipse(scrs, groups, col = plt, lwd = 2, label = TRUE)
  legend("topright", legend = levels(groups), col = plt,
         pch = 19, bty = "n", title = group_var)
  dev.off()
  
  message("Wrote: ", csv_path, " and ", png_path)
  
  invisible(list(
    group_var = group_var,
    permanova = res,
    dispersion = bd,
    var_explained = var_explained,
    csv = csv_path,
    png = png_path
  ))
}

# ---- Run -------------------------------------------------------------

perm_count  <- 9999
dist_method <- "euclidean"

comm <- long_no_control[, 3:ncol(long_no_control)]
d    <- vegdist(comm, method = dist_method)

set.seed(42)  # reproducible NMDS starts *and* permutations

nmds_fit <- metaMDS(comm, distance = dist_method, k = 2, trymax = 100)
check_stress(nmds_fit$stress)

permanova_gosner <- permanova_by_group(
  group_var = "Gosner stage",
  data = tadpoles,
  d = d,
  nmds_fit = nmds_fit,
  perm_count = perm_count,
  out_dir = "permanova"
)

permanova_treatment <- permanova_by_group(
  group_var = "Treatment",
  data = tadpoles,
  d = d,
  nmds_fit = nmds_fit,
  perm_count = perm_count,
  out_dir = "permanova",
  palette = pal[2:4]
)

# Pre-function version
# # Set global variables ----
# perm_count <- 9999
# dist_method <- "euclidean"
# 
# # Calculate distances ----
# d <- vegdist(long_no_control[3:ncol(long_no_control)], method = dist_method)
# 
# # Make sure the relevant columns in tadpoles are factors
# tadpoles <- tadpoles %>% 
#   mutate(
#     `Gosner stage` = as_factor(`Gosner stage`),
#     Treatment = as_factor(Treatment)
#   )
# 
# # Calculate PERMANOVA ----
# permanova_model <- adonis2(d ~ `Gosner stage`,
#                            data = tadpoles,
#                            permutation = perm_count)
# 
# names(permanova_model)[5] <- c("model_p_val")
# permanova_model
# 
# gosner_var_explained <- permanova_model$R2[1] * 100
# 
# # Dispersion check
# bd <- 
#   permutest(
#     betadisper(
#       d,
#       tadpoles$`Gosner stage`),
#     permutations = perm_count)
# 
# permanova_model$dispersion_p_val <- c(bd$tab$`Pr(>F)`[1], NA, NA)
# 
# permanova_model$notes <- c(
#   "Model - R2 indicates the % variance accounted for by the tested variable. Values are from 0 to 1",
#   "model_p_val indiciates significance of the model. < 0.05 is significant",
#   "dispersion_p_val indicates whether the homogeneity of dispersion assumption holds. < 0.05 means it does and the model is good. If the study is balanced then this is less important"
# )
# 
# write.csv(permanova_model, file = "Gosner_PERMANOVA.csv")
# 
# # Ordinance plot ----
# nmds_fit <- 
#   metaMDS(
#     long_no_control[3:ncol(long_no_control)],
#     distance = dist_method,
#     k = 2) # number of dimensions
# 
# stress_level <- nmds_fit$stress
# 
# # Display a warning or message indiciating the appropriateness of the stress level
# 
# if (stress_level >= 0.20) {
#   warning(paste0(
#     "MDS stress of ", 
#     stress_level, 
#     " is above the unreliability threshold of 0.20"))
# } else if (stress_level >= 0.10) {
#   message(paste0(
#     "MDS stress levels are okay (",
#     stress_level,
#     "). Broad interpretation only."))
# } else if (stress_level >= 0.05) {
#   message(paste0(
#     "MDS stress levels are good (",
#     stress_level,
#     "). Safe to interpret all clusters."))
# } else {
#   message(paste0(
#     "MDS stress levels are extremely low (",
#     stress_level,
#     ")! Layout is basically exact."))
# }
# 
# # Make the ordinance plot
# ordinance_plot <- 
#   function() {
#     ordiplot(
#       nmds_fit,
#       type = "n",
#       main = "Gosner Stage PERMANOVA Ordinance Plot",
#       xlim = c(-200, 200))
#     points(
#       nmds_fit,
#       display = "sites",
#       col = 1:6,
#       pch = 19)
#     ordiellipse(
#       nmds_fit,
#       tadpoles$`Gosner stage`,
#       col = 1:6,
#       lwd = 2,
#       label = TRUE)
#   }
# # Save the plot as PNG
# CairoPNG(
#   "Gosner_PERMANOVA_Ordinance.png",
#   width = 1800,
#   height = 1800,
#   dpi = 300
# )
# ordinance_plot()
# dev.off()