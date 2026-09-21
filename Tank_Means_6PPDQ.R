# Isolate the metabolites data from the `merged` data frame and take the mean of
#  each metabolite for each tank
tank_means_data <- merged[4:ncol(merged)] %>% 
  # Transpose data fro easier dplyr operations
  t() %>% as.data.frame() %>%
  # Create an index to identify each tank for grouping
  mutate(Tank_Index = rep(1:24, each = 2), .before = 1) %>% 
  group_by(Tank_Index) %>% 
  # Calculate the means
  summarise_if(is.numeric, mean, na.rm = TRUE) %>% 
  # Delete the index column
  mutate(Tank_Index = NULL) %>% 
  # Return data to its original orientation to rejoin with metadata
  t() %>% as.data.frame()

# Create a vector of new names for each rep (as the previous names included an
# index for the tank that is irrelevant now) and apply it
tank_means_names <- c(
  sprintf("Control %d", 1:6),
  sprintf("High %d", 1:6),
  sprintf("Low %d", 1:6),
  sprintf("Solvent %d", 1:6))
names(tank_means_data) <- tank_means_names

# Add in the Name column from `merged` to use as a joining index
tank_means_data$Name <- merged$Name

# Join the metadata from `merged` with the tank-meaned data
tank_means_merged <- left_join(merged[1:3], tank_means_data)

# Print debugging
print("Tank_Means_6PPDQ.R ran successfully")
