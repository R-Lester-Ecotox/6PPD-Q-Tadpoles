# Subset data to include only "Control" and "Solvent" groups
long_no_6ppdq <- long %>%
  filter(! (Treatment == "Low" | Treatment == "High") )%>% 
  mutate(Treatment = as.character(Treatment),
         Treatment = as_factor(Treatment))
long_no_6ppdq_data_cols <- names(long_no_6ppdq)[3:length(names(long_no_6ppdq))]

# Perform Yuen's trimmed means t-test on each metabolite
t_test_no_6ppdq <- yuen_trimmed_t_test(
  long_no_6ppdq,
  cols = long_no_6ppdq_data_cols,
  class_col = "Treatment"
)

write_csv(t_test_no_6ppdq, file = "results/control_vs_solvent_yuen.csv")

# Isolate significant results
t_test_no_6ppdq_significant <- t_test_no_6ppdq %>% 
  filter(t_test_no_6ppdq$p_val_adj < 0.05)

print("Univariate_Control_vs_Solvent_6PPDQ.R ran succesfully")