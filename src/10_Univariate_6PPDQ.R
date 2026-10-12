# Retain metadata
merged_metadata <- tank_means_merged[, 1:3]

# Create named vector with class information
groupings <- c(
    rep("Control", times = 6),
    rep("High", times = 6),
    rep("Low", times = 6),
    rep("Solvent", times = 6)
  )
names(groupings) <- tank_means_names

# Transpose data an add in class information
long <- tank_means_merged[, 4: ncol(tank_means_merged)] %>% 
  t() %>% as.data.frame() %>%
  set_names(merged_metadata$Name) %>% 
  add_column(groupings, .before = 1) %>% 
  mutate(Treatment = as_factor(groupings),
         groupings = NULL,
         Tank = c(1:24)) %>% 
  relocate(Treatment, Tank, .before = 1)

# # Subset long for testing DELETE BEFORE FINAL
# subset_val <- 100
# 
# long <- long[, 1:(subset_val + 2)]

# Isolate data columns from `long`
long_data_cols <- long %>%
  dplyr::select(-any_of(c("Treatment", "Tank"))) %>%
  names()

# Remove "Control" treatment, to use "Solvent" as point of comparison
long_no_control <- long %>%
  filter(! Treatment == "Control")%>% 
  mutate(Treatment = as.character(Treatment),
         Treatment = as_factor(Treatment),
         Treatment = fct_relevel(Treatment, "Solvent", "Low", "High"))

# Set bootstrap count
boots = 5000

# NOT RUN - bootstrapped tests take a lot of compute. I have computed this once
# and load from a CSV normally.

# # Run a 20% bootstrapped trimmed mean 1-way ANOVA test for each metabolite
# anovas_trimmed_bt <- trimmed_mean_1way_aov_bt_parallel(
#   data = long_no_control,
#   cols = long_data_cols,
#   class_col = "Treatment",
#   nboot = boots
# )
# 
# # Save trimmed mean results as CSV
# write.csv(anovas_trimmed_bt, "trimmed_anova_results_1.csv")

anovas_trimmed_bt <- read.csv("results/trimmed_anova_results_1.csv")

# Sort for significant ANOVAs
significant_anovas <- anovas_trimmed_bt %>% 
  filter(p_val_adj <= 0.05) %>% 
  mutate(
    Method = NA_character_,
    Type = NA_character_
  )

# Lookup and attach the method (i.e. QToF/QQQ) and type (Metabolite/Lipid)
for(x in 1:nrow(significant_anovas)) {
  significant_anovas$Method[x] <- 
    merged$Method[merged$Name == significant_anovas$variable[x]]
  
  significant_anovas$Type[x]  <- 
    merged$Type[merged$Name == significant_anovas$variable[x]]
}

# Post hoc testing ----
sig_metabs <- significant_anovas$variable

# Run WRS2::mcppb20 on all significant ANOVAs
post_hoc_results <- one_way_mcppb20(
  metabolites = sig_metabs,
  data = long_no_control
)

# Isolate the results
post_hocs <- post_hoc_results$post_hocs
clds <- post_hoc_results$clds


print("Univariate_6PPDQ.R ran successfully")

# NOT RUN: Testing script to optimise iterations of `median_aov`
# ~ 3500 seems to be an ideal amount for reaching an accurate p_value. 
# This takes ~ 2.5 mins for 100 metabolites which is incredibly time consuming. 

# nboots <- c(500, 1000, 2000, 3000, 4000, 5000)
# p_vals <- data.frame(matrix(ncol = length(nboots), nrow = 100))
# names(p_vals) <- as.character(nboots)
# 
# for (nboot in nboots) {
#   start <- Sys.time()
#   
#   anovas <- median_aov_parallel(
#     data = long,
#     cols = long_data_cols,
#     class_col = "Treatment",
#     nboot = nboot
#   )
#   
#   end <- Sys.time()
#   
#   time_taken <- as.numeric(end - start)
#   
#   print(sprintf("%d iterations took %f", nboot, time_taken))
#   
#   p_vals[[as.character(nboot)]] <- anovas$p_val
# }
# 
# rownames(p_vals) <- anovas$variable
# 
# p_vals_rounded <- round(p_vals, 3)
# 
# write_csv(p_vals, "med1way_iteration_test.csv")
#  