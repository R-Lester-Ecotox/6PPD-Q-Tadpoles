# Functions ----
# Function to generate a boxplot for each significant metabolite/lipid
# and add compact letter display (CLD) significance annotations
metabolite_boxplots <- function(data,
                                class_col = "Treatment",
                                metabolites,
                                cld_data = NULL,
                                palette = PCA_pal) {
  
  # Create an empty list with one element for each metabolite.
  # The list elements are named using the metabolite names.
  plotlist <- vector(mode = "list", length = length(metabolites)) %>%
    set_names(metabolites)
  
  # Loop through each metabolite/lipid
  for (i in metabolites) {
    if(!is.null(cld_data)) {
      # Create a data frame containing the information needed to
      # position the compact letter display (e.g. a, ab, b) above
      # each treatment group.
      cld_df <- data.frame(
        
        # Extract the treatment names from the CLD object and convert
        # them to factors so they correspond to the x-axis categories
        Treatment =
          names(cld_data[[i]]$Letters) %>%
          as_factor(),
        
        # Extract the significance letters for each treatment
        Letters =
          cld_data[[i]]$Letters,
        
        # Set the vertical position of the letters.
        # The maximum observed value for the metabolite is repeated
        # three times (once for each treatment).
        Height =
          rep(max(data[[i]]), 3)
      )
    }
    # Create the boxplot
    plot <- ggplot(data = data,
                   aes(x = .data[[class_col]], 
                       y = .data[[i]])) +
      
      # Draw a boxplot for each treatment and fill according to
      # the treatment/class column
      geom_boxplot(aes(fill = .data[[class_col]])) +
      
      # Apply the supplied colour palette to the treatment groups
      scale_fill_manual(values = palette) +
      
      # Add the individual observations on top of the boxplots
      geom_point() +
      
      # Add the compact letter display above the boxplots (if there is one)
      {if(!is.null(cld_data)) {
        geom_text(data = cld_df,
                  aes(x = Treatment, 
                      y = Height, 
                      label = Letters),
                  vjust = 1,
                  hjust = 1.5,
                  size = 9)
      }} +
      # Apply the cowplot theme with a relatively large base font size
      theme_cowplot(22) +
      
      # Use the metabolite/lipid name as the plot title
      labs(title = i) +
      
      # Remove unnecessary plot elements because these plots will
      # eventually be assembled into a larger figure
      theme(
        legend.position = "none", 
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        plot.title = element_text(hjust = 0.5)
      )
    
    # Store the completed plot in the list using the metabolite name
    # as the list index
    plotlist[[i]] <- plot
  }
  
  # Return the list of plots
  plotlist
}
# Function to create a row of plots containing a text label.
# Blank plots are added so that the label occupies the same width
# as the plot grid it is being used alongside.
plottext <- function(text, width) {
  
  # Create a blank ggplot to use as a spacer
  blank <- ggplot() + theme_void()
  
  # Create a plot containing only the supplied text
  textplot <- ggdraw() +
    draw_label(
      text,
      x = 0,
      y = 0.5,
      hjust = 0,
      size = 36,
      fontface = "bold"
    )
  
  # Combine the text plot with blank plots to make a row.
  # 'width' determines how many columns the resulting row occupies.
  plot_grid(
    plotlist = c(
      list(textplot),
      rep(list(blank), width - 1)),
    nrow = 1
  )
}

# Set the width/height (in pixels) of an individual plot in the
# final exported figures
plot_side_size <- 320

# Main analysis metabolites-----
# Generate boxplots for the significant metabolites
metabolite_plots <- metabolite_boxplots(
  data = long_no_control,
  metabolites = sig_metabs,
  cld_data = clds)

# Extract the names of significant variables that are classified as metabolites
sig_metabs <- significant_anovas$variable[significant_anovas$Type == "Metabolite"]

# Create the "Metabolites" heading.
# Width = 2 because the metabolite plots will be arranged in 2 columns.
metabo_title <- plottext("Metabolites", width = 2)


# Arrange all metabolite boxplots into a grid with two rows
metabo_plots_merge <- 
  plot_grid(
    plotlist = metabolite_plots,
    nrow = 2)


# Combine the "Metabolites" heading with the metabolite plot grid.
# The heading gets a smaller relative height than the plots themselves.
metabo_plots_final <- plot_grid(
  metabo_title,
  metabo_plots_merge,
  ncol = 1,
  rel_heights = c(40, plot_side_size * 2)
)


# Open a Cairo PNG graphics device.
# The output image is two plots wide and two plots high,
# plus 40 pixels for the title row.
CairoPNG("metabolites_6PPDQ.png",
         width = plot_side_size * 2,
         height = (plot_side_size * 2) + 40)

# Draw the assembled metabolite figure to the PNG device
metabo_plots_final

# Close the graphics device and save the PNG
dev.off()

# Main analysis lipids----
# Extract the names of significant variables that are classified as lipids
sig_lipids <- significant_anovas$variable[significant_anovas$Type == "Lipid"]

# Generate boxplots for the significant lipids
lipid_plots <- metabolite_boxplots(
  data = long_no_control,
  metabolites = sig_lipids,
  cld_data = clds)

# Create the "Lipids" heading.
# Width = 2 because the lipid plots will be arranged in 2 columns.
lipid_title <- plottext("Lipids", width = 2)


# Arrange all lipid boxplots into a single row
lipid_plots_merge <- plot_grid(
  plotlist = lipid_plots,
  nrow = 2)


# Combine the "Lipids" heading with the lipid plot row
lipid_plots_final <- plot_grid(
  lipid_title, 
  lipid_plots_merge, 
  ncol = 1, 
  rel_heights = c(40, plot_side_size * 2)
)

# Open a Cairo PNG graphics device.
# The output image is two plots wide and two plots high,
# plus 40 pixels for the title row.
CairoPNG("lipids_6PPDQ.png",
         width = plot_side_size * 2,
         height = (plot_side_size * 2) + 40)

# Draw the assembled lipid figure to the PNG device
lipid_plots_final

# Close the graphics device and save the PNG
dev.off()

# GS34 stratified metabolites ----
# Extract the names of significant variables that are classified as metabolites
GS34_sig_metabs <- 
  GS34_significant_anovas$variable[GS34_significant_anovas$Type == "Metabolite"]
# Generate boxplots for the significant metabolites
GS34_metabolite_plots <- metabolite_boxplots(
  data = long_no_control,
  metabolites = GS34_sig_metabs,
  cld_data = GS34_clds)


# Create the "Metabolites" heading.
# Width = 2 because the metabolite plots will be arranged in 2 columns.
GS34_metabo_title <- plottext("Metabolites", width = 2)


# Arrange all metabolite boxplots into a grid with two rows
GS34_metabo_plots_merge <- 
  plot_grid(
    plotlist = GS34_metabolite_plots,
    nrow = 2)


# Combine the "Metabolites" heading with the metabolite plot grid.
# The heading gets a smaller relative height than the plots themselves.
GS34_metabo_plots_final <- plot_grid(
  GS34_metabo_title,
  GS34_metabo_plots_merge,
  ncol = 1,
  rel_heights = c(40, plot_side_size * 2)
)


# Open a Cairo PNG graphics device.
# The output image is two plots wide and two plots high,
# plus 40 pixels for the title row.
CairoPNG("GS34_metabolites_6PPDQ.png",
         width = plot_side_size * 2,
         height = (plot_side_size * 2) + 40)

# Draw the assembled metabolite figure to the PNG device
GS34_metabo_plots_final

# Close the graphics device and save the PNG
dev.off()

# GS34 stratified lipids -----
# Extract the names of significant variables that are classified as lipids
GS34_sig_lipids <- 
  GS34_significant_anovas$variable[GS34_significant_anovas$Type == "Lipid"]

# Generate boxplots for the significant lipids
GS34_lipid_plots <- metabolite_boxplots(
  data = long_no_control,
  metabolites = GS34_sig_lipids,
  cld_data = GS34_clds)

# Create the "Lipids" heading.
# Width = 2 because the lipid plots will be arranged in 2 columns.
GS34_lipid_title <- plottext("Lipids", width = 2)


# Arrange all lipid boxplots into two rows
GS34_lipid_plots_merge <- plot_grid(
  plotlist = GS34_lipid_plots,
  nrow = 2)


# Combine the "Lipids" heading with the lipid plot row
GS34_lipid_plots_final <- plot_grid(
  GS34_lipid_title, 
  GS34_lipid_plots_merge, 
  ncol = 1, 
  rel_heights = c(40, plot_side_size * 2)
)

# Open a Cairo PNG graphics device.
# The output image is three plots wide and one plot high,
# plus 40 pixels for the title row.
CairoPNG("GS34_lipids_6PPDQ.png",
         width = plot_side_size * 2,
         height = (plot_side_size * 2) + 40)

# Draw the assembled lipid figure to the PNG device
GS34_lipid_plots_final

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

# Create the "Morphometric indices" heading
# Width = 3 as there will be 3 columns
morpho_title <- plottext("Morphometric indices", width = 3)

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
  morpho_title,
  morpho_plots_row_1,
  morpho_plots_row_2,
  ncol = 1,
  rel_heights = c(40, plot_side_size, plot_side_size)
)

# Open a Cairo PNG graphics device.
# The output image is three plots wide and two plots high,
# plus 40 pixels for the title row.
CairoPNG("Morphometric_indices.png",
         width = plot_side_size * 3,
         height = (plot_side_size * 2) + 40)

# Draw the assembled morphometri indices figure to the PNG device
morpho_plots_final

# Close the graphics device and save the PNG
dev.off()
