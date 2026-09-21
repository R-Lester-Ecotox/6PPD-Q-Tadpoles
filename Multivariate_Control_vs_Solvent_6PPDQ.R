# Principal Component Analysis ----
# Select data
CS_PCA_data <- long_no_6ppdq[, 3:ncol(long)]

# Run the test
CS_PCA_test <- prcomp(CS_PCA_data, center = TRUE, scale. = TRUE)
summary(CS_PCA_test)

# Isolate the variance data
CS_PCA_var <- CS_PCA_test$sdev ^ 2
CS_PCA_percent_var <- CS_PCA_var / sum(CS_PCA_var) * 100

# Convert to dataframe and include an index for the PC
CS_PCA_var <- data.frame(
  PC = 1:length(CS_PCA_var),
  var = CS_PCA_var,
  per_cent_var = CS_PCA_percent_var
)

# Create a colour palette
CS_PCA_pal <- c(Solvent = "#bebebe", Control = "grey15")

# Consitent theme elements
CS_PCA_theme <-
  theme_cowplot(14) +
  theme(
    axis.ticks = element_blank(),
    plot.background = element_rect(
      colour = "black",
      fill = NA,
      linewidth = 0.5))

# Create a bar blot of the explained variance
CS_PCA_bar <- CS_PCA_var[1:10, ] %>% 
  ggplot(aes(x = as.factor(PC), y = per_cent_var)) +
  geom_col(fill = CS_PCA_pal[[2]]) +
  scale_y_continuous(expand = c(0, 0)) +
  xlab("PC") +
  ylab("% of total variance") +
  CS_PCA_theme

CS_PCA_bar

# Create a scatterplot of the first two principal components ----
# Wrangle data
CS_PCA_PCs <- CS_PCA_test$x %>% 
  as.data.frame() %>% 
  mutate(Treatment = long_no_6ppdq$Treatment, .before = 1)

# Variables to choose PCs
PC_x <- 1
PC_y <- 2

# Generate the plot
CS_PCA_lims <- range(with(CS_PCA_PCs, c(PC1, PC2)))

CS_PCA_scatter <- CS_PCA_PCs %>% 
  ggplot(
    aes(
      x = .data[[names(CS_PCA_PCs)[1+PC_x]]],
      y = .data[[names(CS_PCA_PCs)[1+PC_y]]],
      colour = Treatment)) +
  geom_point(size = 5) +
  coord_cartesian(xlim = CS_PCA_lims, ylim = CS_PCA_lims) +
  xlab(sprintf("PC%d (%.2f%%)", PC_x, CS_PCA_var$per_cent_var[PC_x])) +
  ylab(sprintf("PC%d (%.2f%%)", PC_y, CS_PCA_var$per_cent_var[PC_y])) +
  scale_color_manual(values = CS_PCA_pal, name = "Treatment") +
  CS_PCA_theme +
  theme(
    legend.position = "top",
    legend.justification.top = "centre",
    legend.title = element_blank())

CS_PCA_scatter

# Loadings plot ----

CS_PCA_loadings <- CS_PCA_test$rotation %>% as.data.frame()

CS_PCA_ldg_breaks <-
  seq(
    signif(min(PCA_loadings$PC1), digits = 1),
    signif(max(PCA_loadings$PC1), digits = 1),
    length.out = 3
  )

CS_PCA_ldg_plot <- CS_PCA_loadings %>% 
  ggplot(
    aes(
      x= .data[[names(CS_PCA_loadings)[PC_x]]],
      y = .data[[names(CS_PCA_loadings)[PC_y]]])) +
  xlab(sprintf("PC%d (%.2f%%)", PC_x, CS_PCA_var$per_cent_var[PC_x])) +
  ylab(sprintf("PC%d (%.2f%%)", PC_y, CS_PCA_var$per_cent_var[PC_y])) +
  geom_point(
    size = 1,
    colour = CS_PCA_pal[[2]]) +
  CS_PCA_theme +
  scale_x_continuous(breaks = CS_PCA_ldg_breaks) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1))

CS_PCA_ldg_plot

# Merge all the plots together ----
CS_PCA_small_plots <- plot_grid(
  CS_PCA_bar, 
  CS_PCA_ldg_plot, 
  ncol = 1)
CS_PCA_small_plots

CS_PCA_all_plots <- plot_grid(
  CS_PCA_scatter, 
  CS_PCA_small_plots,
  rel_widths = c(1, 0.5))
CS_PCA_all_plots

# Save the final multipanel plot ----
px <- 1500

CairoPNG(
  filename = "Control_vs_Solvent_PCA.png",
  width = px * 1.5,
  height = px,
  units = "px",
  dpi = 300
)
CS_PCA_all_plots
dev.off()

print("Multivariate_Control_vs_Solvent_6PPDQ.R ran successfully")