# Principal Component Analysis ----
# Select data
PCA_data <- long_no_control[, 3:ncol(long)]

# Run the test
PCA_test <- prcomp(PCA_data, center = TRUE, scale. = TRUE)
summary(PCA_test)

# Isolate the variance data
PCA_var <- PCA_test$sdev ^ 2
PCA_percent_var <- PCA_var / sum(PCA_var) * 100

# Convert to dataframe and include an index for the PC
PCA_var <- data.frame(
  PC = 1:length(PCA_var),
  var = PCA_var,
  per_cent_var = PCA_percent_var
)

# TO DO:
# • Theme cowplot
# • Label each bar with the actual number
# • Remove dead space underneath
# • x-axis to 100

# Create a colour palette
PCA_pal <- c(Solvent = "#bebebe", Low = "#278192", High = "#3f3f7b")

# Consitent theme elements
PCA_theme <-
  theme_cowplot(14) +
  theme(
    axis.ticks = element_blank(),
    plot.background = element_rect(
      colour = "black",
      fill = NA,
      linewidth = 0.5))

# Create a bar blot of the explained variance
PCA_bar <- PCA_var[1:10, ] %>% 
  ggplot(aes(x = as.factor(PC), y = per_cent_var)) +
  geom_col(fill = PCA_pal[[2]]) +
  scale_y_continuous(expand = c(0, 0)) +
  xlab("PC") +
  ylab("% of total variance") +
  PCA_theme

PCA_bar

# Create a scatterplot of the first two principal components ----
# Wrangle data
PCA_PCs <- PCA_test$x %>% 
  as.data.frame() %>% 
  mutate(Treatment = long_no_control$Treatment, .before = 1)

# Variables to choose PCs
PC_x <- 1
PC_y <- 2

# Generate the plot
PCA_lims <- range(with(PCA_PCs, c(PC1, PC2)))

PCA_scatter <- PCA_PCs %>% 
  ggplot(
    aes(
      x = .data[[names(PCA_PCs)[1+PC_x]]],
      y = .data[[names(PCA_PCs)[1+PC_y]]],
      colour = Treatment)) +
  geom_point(size = 5) +
  coord_cartesian(xlim = PCA_lims, ylim = PCA_lims) +
  xlab(sprintf("PC%d (%.2f%%)", PC_x, PCA_var$per_cent_var[PC_x])) +
  ylab(sprintf("PC%d (%.2f%%)", PC_y, PCA_var$per_cent_var[PC_y])) +
  scale_color_manual(values = PCA_pal, name = "Treatment") +
  PCA_theme +
  theme(
    legend.position = "top",
    legend.justification.top = "centre",
    legend.title = element_blank())

PCA_scatter

# Loadings plot ----

PCA_loadings <- PCA_test$rotation %>% as.data.frame()

PCA_ldg_breaks <-
  seq(
    signif(min(PCA_loadings$PC1), digits = 2),
    signif(max(PCA_loadings$PC1), digits = 2),
    length.out = 3
  )

PCA_ldg_plot <- PCA_loadings %>% 
  ggplot(
    aes(
      x= .data[[names(PCA_loadings)[PC_x]]],
      y = .data[[names(PCA_loadings)[PC_y]]])) +
  xlab(sprintf("PC%d (%.2f%%)", PC_x, PCA_var$per_cent_var[PC_x])) +
  ylab(sprintf("PC%d (%.2f%%)", PC_y, PCA_var$per_cent_var[PC_y])) +
  geom_point(
    size = 1,
    colour = PCA_pal[[2]]) +
  PCA_theme +
  scale_x_continuous(breaks = PCA_ldg_breaks) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1))

PCA_ldg_plot

# Merge all the plots together ----
PCA_small_plots <- plot_grid(
  PCA_bar, 
  PCA_ldg_plot, 
  ncol = 1)
PCA_small_plots

PCA_all_plots <- plot_grid(
  PCA_scatter, 
  PCA_small_plots,
  rel_widths = c(1, 0.5))
PCA_all_plots

# Save the final multipanel plot ----
px <- 1500

CairoPNG(
  filename = "PCA.png",
  width = px * 1.5,
  height = px,
  units = "px",
  dpi = 300
)
PCA_all_plots
dev.off()

print("Multivariate_6PPDQ.R ran successfully")

