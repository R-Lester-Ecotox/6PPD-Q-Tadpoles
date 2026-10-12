# Create a vector to specify which columns are data (i.e. *not* metadata)
data_cols_merged <- data_cols_CCM[5:length(data_cols_CCM)]

# Define the groups. Necessary for `eighty_per_cent_rule`
reps <- 12
data_cols_merged_grouping <- 
  c(
    rep("Control", reps),
    rep("High", reps),
    rep("Low", reps),
    rep("Solvent", reps)) %>% 
  setNames(data_cols_merged)

# Create a named vector of the dry weights. Req. for `weight_normalise` 
dry_weights <- metadata$`Dry Weight mg` %>% 
  setNames(paste0("Tadpole ", metadata$Sample))

# Merge data and perform final preprocessing steps. These steps are the same for
# all datasets, and are done together to simplify the script and minimise compute
merged <- bind_rows(CCM_processed, QToF_Lipid_processed, QToF_Metabo_processed) %>%
  select(-matches("AAOA|Recon|PBQC|NoInj")) %>% 
  eighty_per_cent_rule(
    cols = data_cols_merged, 
    groups = data_cols_merged_grouping) %>% 
  remove_duplicate_metabolites() %>%
  weight_normalise(
    cols = data_cols_merged,
    weights = dry_weights) %>% 
  median_normalise() %>%
  pareto_scale() %>% 
  QRILC_impute_with_metadata(cols = data_cols_merged)
# Clean up the names
names(merged) <- str_remove(names(merged), "Tadpole ")

# Print debugging
print("Merge_6PPDQ.R ran successfully")


