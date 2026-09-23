# Functions ----
# Function to generate a boxplot for each significant metabolite/lipid
# and add compact letter display (CLD) significance annotations
metabolite_boxplots <- function(data,
                                class_col = "Treatment",
                                metabolites,
                                cld_data = NULL,
                                palette = pal) {
  
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
