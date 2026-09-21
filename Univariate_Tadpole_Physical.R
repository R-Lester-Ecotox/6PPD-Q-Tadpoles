# Data wrangling ----
# Sort out factors nicely
tadpoles <- metadata %>%
  mutate(
    Treatment =
      str_replace_all(Treatment, c("1 ug/L" = "Low", "10 ug/L" = "High")),
    Treatment = as_factor(Treatment) %>%
      fct_relevel(factor_order),
    Tank = as_factor(Tank)) %>%
  filter(! Treatment == "Control") %>% 
  # Then take tank means
  group_by(Tank, Treatment) %>%
  summarise_if(is.numeric, mean, na.rm = TRUE) %>%
  arrange(Treatment)

# Kruskal-Wallis rank sums ----
# Specify the columns to test with KW
physical_data_cols <- names(tadpoles)[3:ncol(tadpoles)]

# Initialise a list to receive the data
physical_kw <- vector(
  mode = "list",
  length = length(physical_data_cols)
) %>% setNames(physical_data_cols)

# KW test for each variable. Storing data neatly in dataframes
for(var in physical_data_cols) {
  formula <- as.formula(paste0("`", var, "`~Treatment"))
  
  kw <- kruskal.test(
    formula = formula,
    data = tadpoles)
  
  physical_kw[[var]] <- data.frame(
    chi.sq = kw$statistic[[1]],
    df     = kw$parameter[[1]],
    p.value = kw$p.value
  )
}

# Combine dataframes and arrange by p.value
physical_kw <- do.call(rbind, physical_kw) %>%
  rownames_to_column(var = "Parameter") %>% 
  arrange(p.value)

# Summary statistics----
# Initialise a list to receive the data
physical_summary_stats <- vector(
  mode = "list",
  length = length(physical_data_cols)
) %>% setNames(physical_data_cols)

# Generate summary stats for each variable. Storing data neatly in 1-row dataframes
for(var in physical_data_cols) {
  physical_summary_stats[[var]] <- 
    summary(tadpoles[[var]]) %>% 
    as.matrix() %>% 
    t() %>% 
    as.data.frame() %>% 
    mutate(Variable = c(var),
           SD = sd(tadpoles[[var]]),
           .before = 1) %>% 
    relocate(Mean, .before = 2)
}

# Join all the 1-row dataframes together
physical_summary_stats <- bind_rows(physical_summary_stats)

# Save the results to CSV
write.csv(physical_kw, "physical_kruskal-wallis.csv")
write.csv(physical_summary_stats, "physical_summary_stats.csv")

# # Post-hoc testing ----
# # <NOT RUN> as Gosner stage isn't significant when "Control" is filtered out
# # Post-hoc test for Gosner stage, which is nearly significant (p ~ 0.06)
# gosner_posthoc_model <- mcppb20(
#   `Gosner stage` ~ `Treatment`,
#   data = tadpoles)
# 
# # Compact letter display ----
# # Need the treatment names in order for compact letter display
# gosner_fnames <- gosner_posthoc_model$fnames
# names(gosner_fnames) <- as.character(c(1:length(gosner_fnames)))
# 
# # Isolate the useful data
# gosner_posthoc_result <- gosner_posthoc_model$comp %>% 
#   as.data.frame()
# 
# # Condense the grouping columns into the right format for `multCompLetters`
# names(gosner_posthoc_result)[1:2] <- 
#   c("Comparison_Group_1", "Comparison_Group_2")
# 
# gosner_posthoc_result <- gosner_posthoc_result %>%
#   mutate(
#     Comparison_Group_1 = as.character(Comparison_Group_1),
#     Comparison_Group_1 = str_replace_all(Comparison_Group_1, gosner_fnames),
#     Comparison_Group_2 = as.character(Comparison_Group_2),
#     Comparison_Group_2 = str_replace_all(Comparison_Group_2, gosner_fnames),
#     Comparison = paste0(Comparison_Group_1, "-", Comparison_Group_2),
#     Comparison_Group_1 = NULL,
#     Comparison_Group_2 = NULL
#   )
# 
# # Named vector of p-values, with Comparisons as names for multCompLetters
# gosner_cld_data <- gosner_posthoc_result$`p-value`
# names(gosner_cld_data) <- gosner_posthoc_result$Comparison
# 
# # Caclualte the compact letter display
# gosner_cld <- multcompLetters(
#   gosner_cld_data, threshold = 0.05)

# Print debugging ----
print("Univariate_Tadpole_Physical.R ran successfully")