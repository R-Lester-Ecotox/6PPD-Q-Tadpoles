# Create a vector to specify which columns are data (i.e. *not* metadata)
data_cols_qtof_lipid <- 
  get_metabolite_columns(
    neg_mode = QToF_Lipid_NEG, 
    pos_mode = QToF_Lipid_POS, 
    pattern = "Tadpole")

# Prepare the data for analysis using the nice tidy functions from QToF_Funcs.R
QToF_Lipid_processed <-
  POS_NEG_Merge(
    pos_mode_data = QToF_Lipid_POS,
    neg_mode_data = QToF_Lipid_NEG) %>%
  NA_2_Unknown() %>%
  remove_low_scores(scores = .$`Score(MFE,max)`, threshold = 80) %>%
  remove_high_mass_error(masserror = .$`RSD(Mass, ppm)`) %>%
  mutate(baseline1 = calculate_baseline(., pattern = "NoInj")) %>%
  subtract_baseline(cols = data_cols_qtof_metabo, baseline_col = "baseline1") %>%
  mutate(baseline2 = calculate_baseline(., pattern = "Recon")) %>% 
  subtract_baseline(cols = data_cols_qtof_metabo, baseline_col = "baseline2") %>% 
  sum_same_metabolites(cols = data_cols_qtof_lipid) %>%
  low_count_filter(cols = data_cols_qtof_lipid) %>%
  dplyr::select(contains(c("Name", "Method", "Tadpole"))) %>%
  mutate(Type = "Lipid", .after = 2)

# Print debugging
print("QToF_Lipids_6PPDQ.R ran successfully")