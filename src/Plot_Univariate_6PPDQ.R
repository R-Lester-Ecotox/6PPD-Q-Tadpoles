# Set the width/height (in pixels) of an individual plot in the
# final exported figures
plot_side_size <- 320
set_scale <- 5
set_dpi <- 300

# Main analysis metabolites-----
# Generate boxplots for the significant metabolites
metabolite_plots <- metabolite_boxplots(
  data = long_no_control,
  metabolites = sig_metabs,
  cld_data = clds,
  title_size = 12,
  wrap_width = 40)

# Extract the names of significant variables that are classified as metabolites
sig_metabs <- significant_anovas$variable[significant_anovas$Type == "Metabolite"]

# Arrange all metabolite boxplots into a grid with two rows
metabo_plots_merge <- 
  plot_grid(
    plotlist = metabolite_plots,
    nrow = 1)

# Combine the "Metabolites" heading with the metabolite plot grid.
# The heading gets a smaller relative height than the plots themselves.
metabo_plots_final <- plot_grid(
  # metabo_title,
  metabo_plots_merge,
  ncol = 1,
  rel_heights = c(40, plot_side_size)
)

metabo_plots_final
# Open a Cairo PNG graphics device.
# The output image is one plot wide and one plot high,
CairoPNG("plots/metabolites_6PPDQ.png",
         width = plot_side_size * set_scale,
         height = (plot_side_size + 40) * set_scale,
         dpi = set_dpi)

# Draw the assembled metabolite figure to the PNG device
metabo_plots_final

# Close the graphics device and save the PNG
dev.off()

# Morphometric indices----
# Extract the names of the morphometric indices
names(tadpoles) <- c("Tank", "Treatment", "Wet Weight (mg)", "Dry Weight (mg)",
                     "Gosner stage", "SVL (mm)", "Full length (mm)")
morphos <- names(tadpoles)[!(names(tadpoles) %in% c("Tank", "Treatment"))]

# Generate boxplots of all morphometric indices
morpho_plots <- metabolite_boxplots(
  tadpoles, 
  metabolites = morphos)

# Generate a blank spacer plot for padding
spacer = ggplot() + theme_void()

# Arrange the first three morphos into one row ...
morpho_plots_row_1 <- plot_grid(
  plotlist = morpho_plots[1:3],
  nrow = 1)

# ...and the remaining two into another row
morpho_plots_row_2 <- plot_grid(
  spacer,
  morpho_plots[[4]],
  morpho_plots[[5]],
  spacer,
  nrow = 1,
  rel_widths = c(0.5, 1, 1, 0.5))

# Combine the morpho heading with the two plot rows
morpho_plots_final <- plot_grid(
  morpho_plots_row_1,
  morpho_plots_row_2,
  ncol = 1,
  rel_heights = c(plot_side_size, plot_side_size)
)

# Open a Cairo PNG graphics device.
# The output image is three plots wide and two plots high,
# plus 40 pixels for the title row.
CairoPNG("plots/Morphometric_indices.png",
         width = plot_side_size * 3 * set_scale,
         height = plot_side_size * 2 * set_scale,
         dpi = set_dpi)

# Draw the assembled morphometric indices figure to the PNG device
morpho_plots_final

# Close the graphics device and save the PNG
dev.off()

print("Plot_Univariate_6PPDQ.R ran successfully")

