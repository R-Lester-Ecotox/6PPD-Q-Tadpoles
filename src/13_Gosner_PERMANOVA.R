#' PERMANOVA + dispersion check + PCA ordination plot for one grouping
#' variable.
#'
#' @param group_var Column name in `data` (character, backticks not needed).
#' @param data Sample metadata, one row per sample, in the same order as `d`.
#' @param d A `dist` object of between-sample distances.
#' @param pca_fit A `prcomp` object.
#' @param perm_count Number of permutations.
#' @param out_dir Directory for the CSV and PNG outputs.
#' @param plot_title Optional plot title; defaults to "<group_var> ...".
#' @param file_prefix Optional filename stem; defaults to a sanitised
#'   `group_var`.
#' @return Invisibly, a list with the PERMANOVA table, dispersion test,
#'   variance explained and output paths.
permanova_by_group <- function(
    group_var,
    data,
    d,
    pca_fit,
    perm_count = 9999,
    out_dir = ".",
    plot_title = NULL,
    file_prefix = NULL,
    xlim = NULL,
    ylim = NULL,
    width = 1800,
    height = 1800,
    dpi = 300,
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
    plot_title <- paste(group_var, "PCA Ordination")
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
  
  # PERMANOVA ----
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
  
  # Dispersion check ----
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
  
  # Ordination plot ----
  scrs <- pca_fit$x[, 1:2]
  
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
    cex = 1.2,
    cex.axis = 1.0,
    cex.lab = 1.4,
    cex.main = 1.5,
    cex.sub = 1.2
  )
  
  pc_var <- summary(pca_fit)$importance[2, 1:2] * 100
  
  plot(scrs, type = "n", asp = NULL,
       xlim = xlim, ylim = ylim,
       xlab = sprintf("PC1 (%.1f%%)", pc_var[1]), 
       ylab = sprintf("PC2 (%.1f%%)", pc_var[2]), 
       main = plot_title,
       sub = sprintf("R2 = %.1f%% | p = %.4f",
                      var_explained, res$model_p_val[1]))
  
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

# Run ----
comm <- long_no_control[, 3:ncol(long_no_control)]

perm_count <- 9999

dist_method = "euclidean"

d <- vegdist(comm, method = dist_method)

pca_fit <- prcomp(comm, center = TRUE, scale. = FALSE)

scrs <- pca_fit$x[, 1:2]

permanova_gosner <- permanova_by_group(
  group_var = "Gosner stage",
  data = tadpoles,
  d = d,
  pca_fit = pca_fit,
  perm_count = perm_count,
  out_dir = "permanova"
)

permanova_treatment <- permanova_by_group(
  group_var = "Treatment",
  data = tadpoles,
  d = d,
  pca_fit = pca_fit,
  perm_count = perm_count,
  out_dir = "permanova",
  palette = pal[2:4]
)

print("Gosner_PERMANOVA.R ran successfully")
