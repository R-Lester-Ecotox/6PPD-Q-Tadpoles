# Make results reproducible
set.seed("13")
# ---- Prepare data for mixOmics ----
# X matrix: samples in rows, features (metabolites/lipids) in columns
mixomics_merged <- tank_means_merged %>%
  # Remove any existing row names so column_to_rownames() doesn't clash
  `rownames<-`(NULL) %>%
  # Use the feature name column as row names so the data are purely numeric
  column_to_rownames("Name") %>%
  # Drop non-numeric annotation columns (e.g. Type, Method)
  select(-contains(c("Type", "Method"))) %>%
  # Transpose so samples are rows and features are columns, as mixOmics expects
  t()

# Y vector: class label for each sample (6 tanks per treatment).
# The order MUST match the row order of mixomics_merged.
mixomics_class <- c(
  rep("C", 6),
  rep("H", 6),
  rep("L", 6),
  rep("S", 6)
)

# ---- Unsupervised overview: PCA ----
# Explore the main sources of variation without using the class labels.
# scale = FALSE keeps features on their original scale (assumes data are
# already scaled/transformed upstream).
pca <- mixOmics::pca(mixomics_merged, ncomp = 3, scale = FALSE)

# Plot samples on the first two PCs, coloured by treatment
plotIndiv(pca, group = mixomics_class, ind.names = F, legend = T)

# ---- Supervised model: PLS-DA ----
# Fit an initial PLS-DA with extra components so we can assess how many
# are actually needed
plsda <- mixOmics::plsda(
  mixomics_merged,
  mixomics_class,
  ncomp = 5)

# Cross-validate the PLS-DA (5-fold, repeated 100 times) to estimate
# classification error as a function of the number of components
perf.plsda <- perf(
  plsda,
  validation = "Mfold",
  folds = 5,
  progressBar = T,
  nrepeat = 100)

# Plot error rates (with SD) to choose the optimal number of components
plot(perf.plsda, sd = TRUE, legend.position = "horizontal")

# ---- Sparse PLS-DA: feature selection ----
# Candidate numbers of features to keep per component:
# 1-10 in steps of 1, then 20-100 in steps of 10
list.keepX <- c(1:10, seq(20, 100, 10))

# Tune the number of features to retain on each component using
# cross-validation (5-fold, repeated 50 times) and the max distance
# classifier. This overwrites the tune.splsda() function name with its result.
tune.splsda <- tune.splsda(
  mixomics_merged,
  mixomics_class,
  ncomp = 3,
  validation = "Mfold",
  folds = 5,
  dist = "max.dist",
  test.keepX = list.keepX,
  nrepeat = 50,
  progressBar = TRUE
)

# Shows the optimum keepX values according to tuning
plot(tune.splsda, sd = T, title)

# Extract the optimal number of features per component
select.keepX <- tune.splsda$choice.keepX

# Fit the final sparse PLS-DA using the tuned keepX values
splsda <- mixOmics::splsda(
  mixomics_merged,
  mixomics_class,
  ncomp = 2,
  keepX = select.keepX[1:2]
)

# ---- Evaluate the final sPLS-DA ----
# Cross-validate the final model to get classification error rates
# and feature selection stability
perf.splsda <- perf(
  splsda,
  validation = "Mfold",
  folds = 5,
  dist = "max.dist",
  nrepeat = 50,
  progressBar = TRUE
)

# Per-class error rates for each component
perf.splsda$error.rate.class

# ---- Selected features and their stability ----
# Names of the features selected on components 1 and 2
select.name1 <- selectVar(splsda, comp = 1)$name
select.name2 <- selectVar(splsda, comp = 2)$name

# Stability = proportion of cross-validation folds in which each selected
# feature was chosen (closer to 1 = more robust)
stability1 <- perf.splsda$features$stable$comp1[select.name1]
stability2 <- perf.splsda$features$stable$comp2[select.name2]

# Combine loading values with stability for each component
vips_comp1 <- cbind(selectVar(splsda, comp = 1)$value, stability1)
vips_comp2 <- cbind(selectVar(splsda, comp = 2)$value, stability2)

# ---- Visualise the final model ----
# Sample plot on components 1 and 2 with 95% confidence ellipses per group
plotIndiv(
  splsda,
  comp = c(1,2),
  ind.names = FALSE,
  ellipse = TRUE,
  legend = TRUE,
  star = FALSE,
  title = "sPLSDA scores plot"
)
