# Load and wrangle data
read_qtof_tsv <- function(filename) {
  result <- 
    read_tsv(
      file = filename,
      locale = locale(encoding="latin1")) %>%
    dplyr::select(-contains("injection"))
}

QToF_Metabo_NEG <- read_qtof_tsv("Tadpole_QToF_NEG.tsv")
QToF_Metabo_POS <- read_qtof_tsv("Tadpole_QToF_POS.tsv")
QToF_Lipid_NEG <- read_qtof_tsv("Tadpole_Lipid_NEG.tsv")
QToF_Lipid_POS <- read_qtof_tsv("Tadpole_Lipid_POS.tsv")

names <- read_csv("metabolites_harmonised.csv")

CCM <- read_csv("Tadpoles_QQQ_Metabolites.csv") %>% 
  t() %>% 
  as.data.frame() %>% 
  row_to_names(1) %>%
  mutate(across(everything(), as.numeric)) %>% 
  rownames_to_column(var = "Name") %>% 
  left_join(names[1:2], by = c("Name" = "CCM")) %>% 
  mutate(Name = CCM_harmonised) %>% 
  dplyr::select(-CCM_harmonised)

factor_order <- c("Control", "Solvent", "Low", "High")

metadata <- read_csv("tadpole.csv")

print("Import_6PPDQ.R ran succesfully")
