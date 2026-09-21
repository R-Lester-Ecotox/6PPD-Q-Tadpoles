# Create a vector to specify which columns are data (i.e. *not* metadata)
data_cols_CCM_logic <- names(CCM) %>% str_detect(c("PBQC|AAOA|Name|NoInj"))

# Subset the names of CCM by this string. To pass to data pre-processing funcs.
data_cols_CCM <- names(CCM)[!data_cols_CCM_logic]

# Pre-process CCM data
CCM_processed <- CCM %>%
  mutate(baseline = calculate_baseline(., pattern = "NoInj")) %>% 
  subtract_baseline(cols = data_cols_CCM, baseline_col = "baseline") %>%
  mutate(baseline = calculate_baseline(., pattern = "Recon")) %>% 
  subtract_baseline(cols = data_cols_CCM, baseline_col = "baseline") %>%
  low_count_filter(cols = data_cols_CCM, threshold = 300) %>%
  mutate(Method = "CCM_QQQ", Type = "Metabolite", .after = 1)

# Print error checking
print("CCM_6PPDQ.R ran succesfully")
