# Stratify the dataset and remove controls
GS34_metadata <- metadata[metadata$`Gosner stage` == 34 
                          & ! metadata$Treatment == "Control",]

# Subset the merged dataframe to include metadata columns and the samples from
# GS34_metadata$Sample
GS34_merged <- merged[ , c(names(merged)[1:3],
                           GS34_metadata$Sample)]

GS34_grouping <- GS34_metadata$Treatment
names(GS34_grouping) <- GS34_metadata$Sample

# Get GS34 data organised into long format
GS34_long <- GS34_merged[, 4: ncol(GS34_merged)] %>% 
  t() %>% as.data.frame() %>%
  set_names(merged_metadata$Name) %>% 
  add_column(GS34_grouping, .before = 1) %>% 
  mutate(Treatment = as_factor(GS34_grouping),
         GS34_grouping = NULL,
         Tank = c(1:length(GS34_grouping))) %>% 
  relocate(Treatment, Tank, .before = 1)

# Isolate data columns from `long`
GS34_long_data_cols <- GS34_long %>%
  dplyr::select(-any_of(c("Treatment", "Tank"))) %>%
  names()

# # Run a 20% bootstrapped trimmed mean 1-way ANOVA test for each metabolite. MODIFY BEFORE FINAL
# GS34_anovas_trimmed_bt <- trimmed_mean_1way_aov_bt_parallel(
#   data = long_no_control,
#   cols = long_data_cols,
#   class_col = "Treatment",
#   nboot = boots
# )
# 
# write_csv(GS34_anovas_trimmed_bt, "trimmed_anova_gosner_stratified.csv")
GS34_anovas_trimmed_bt <- read_csv("results/trimmed_anova_gosner_stratified.csv")

# Sort for significant ANOVAs
GS34_significant_anovas <- GS34_anovas_trimmed_bt %>% 
  filter(p_val_adj <= 0.05) %>% 
  mutate(
    Method = NA_character_,
    Type = NA_character_
  )
# NOT RUN: No significant ANOVAs so post-hoc testing is unnecessary
# # Lookup and attach the method (i.e. QToF/QQQ) and type (Metabolite/Lipid)
# for(x in 1:nrow(GS34_significant_anovas)) {
#   GS34_significant_anovas$Method[x] <- 
#     merged$Method[merged$Name == GS34_significant_anovas$variable[x]]
#   
#   GS34_significant_anovas$Type[x]  <- 
#     merged$Type[merged$Name == GS34_significant_anovas$variable[x]]
# }
# 
# # Post hoc testing ----
# GS34_sig_metabs <- GS34_significant_anovas$variable
# 
# 
# # Run WRS2::mcppb20 on all significant ANOVAs
# GS34_post_hoc_results <- one_way_mcppb20(
#   metabolites = GS34_sig_metabs,
#   data = long_no_control
# )
# 
# # Isolate the results
# GS34_post_hocs <- GS34_post_hoc_results$post_hocs
# GS34_clds <- GS34_post_hoc_results$clds

print("Unviariate_Gosner_Stratified_6PPDQ.R ran succesfully")