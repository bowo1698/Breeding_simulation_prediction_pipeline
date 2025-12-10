## Preprocessing data

# 1.  Handle missing values ​​in genotype test population (imputation or remove)
# 2.  Genotype matrix standardization (center & scale)
# 3.  Quality control: MAF filtering, call rate filtering

#install.packages(c("tidyverse", "randomForest", "xgboost", 
#                   "e1071", "rrBLUP", "gridExtra", "hibayes"))

# Load libraries
suppressPackageStartupMessages({
  library(tidyverse)
  library(randomForest)
  library(xgboost)
  library(e1071)
  library(hibayes)
  library(rrBLUP)
  library(gridExtra)
  library(yaml)
})

# load config and extract
config <- read_yaml("config.yaml")
# Extract parameters
sim_output_dir <- config$output$output_dir
test_gens <- config$breeding$test_generations

# Create output directory 
run_name <- basename(sim_output_dir)
output_dir <- file.path("prediction_output", run_name)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

message("Output directory created:", output_dir, "\n\n")

# Load data
genotype_ref <- read_csv(file.path(sim_output_dir, "genotype_reference.csv"), show_col_types = FALSE)
phenotype_ref <- read_csv(file.path(sim_output_dir, "phenotype_reference.csv"), show_col_types = FALSE)
genotype_test <- read_csv(file.path(sim_output_dir, paste0("genotype_test_gen", test_gens[1], ".csv")), show_col_types = FALSE)
phenotype_test <- read_csv(file.path(sim_output_dir, paste0("phenotype_test_gen", test_gens[1], ".csv")), show_col_types = FALSE)
genotype_test_gen17 <- read_csv(file.path(sim_output_dir, paste0("genotype_test_gen", test_gens[2], ".csv")), show_col_types = FALSE)
phenotype_test_gen17 <- read_csv(file.path(sim_output_dir, paste0("phenotype_test_gen", test_gens[2], ".csv")), show_col_types = FALSE)
genotype_test_gen18 <- read_csv(file.path(sim_output_dir, paste0("genotype_test_gen", test_gens[3], ".csv")), show_col_types = FALSE)
phenotype_test_gen18 <- read_csv(file.path(sim_output_dir, paste0("phenotype_test_gen", test_gens[3], ".csv")), show_col_types = FALSE)

# Separate ID columns
ref_ids <- genotype_ref[[1]]
test_ids <- genotype_test[[1]]
test_ids_gen17 <- genotype_test_gen17[[1]]
test_ids_gen18 <- genotype_test_gen18[[1]]

# Extract genotype matrices
geno_ref <- data.matrix(genotype_ref[, -1])
geno_test <- data.matrix(genotype_test[, -1])
geno_test_gen17 <- data.matrix(genotype_test_gen17[, -1])
geno_test_gen18 <- data.matrix(genotype_test_gen18[, -1])

# 1. Quality control
maf_threshold <- 0.05
callrate_threshold <- 0.95

# Calculate MAF from reference
allele_freq <- colMeans(geno_ref, na.rm = TRUE) / 2
maf <- pmin(allele_freq, 1 - allele_freq)

# Calculate call rate from test
callrate <- colMeans(!is.na(geno_test))

# Filter SNPs
keep_snps <- (maf >= maf_threshold) & (callrate >= callrate_threshold)

# Apply filter
geno_ref_filtered <- geno_ref[, keep_snps]
geno_test_filtered <- geno_test[, keep_snps]
geno_test_gen17_filtered <- geno_test_gen17[, keep_snps]
geno_test_gen18_filtered <- geno_test_gen18[, keep_snps]

# 2. Missing value imputation filtering
# Impute missing values in REFERENCE population
for(i in 1:ncol(geno_ref_filtered)) {
  if(any(is.na(geno_ref_filtered[, i]))) {
    geno_ref_filtered[is.na(geno_ref_filtered[, i]), i] <- mean(geno_ref_filtered[, i], na.rm = TRUE)
  }
}

# Impute missing values in TEST population
for(i in 1:ncol(geno_test_filtered)) {
  if(any(is.na(geno_test_filtered[, i]))) {
    geno_test_filtered[is.na(geno_test_filtered[, i]), i] <- mean(geno_test_filtered[, i], na.rm = TRUE)
  }
}

# Impute missing values in TEST Gen 17
for(i in 1:ncol(geno_test_gen17_filtered)) {
  if(any(is.na(geno_test_gen17_filtered[, i]))) {
    geno_test_gen17_filtered[is.na(geno_test_gen17_filtered[, i]), i] <- mean(geno_test_gen17_filtered[, i], na.rm = TRUE)
  }
}

# Impute missing values in TEST Gen 18
for(i in 1:ncol(geno_test_gen18_filtered)) {
  if(any(is.na(geno_test_gen18_filtered[, i]))) {
    geno_test_gen18_filtered[is.na(geno_test_gen18_filtered[, i]), i] <- mean(geno_test_gen18_filtered[, i], na.rm = TRUE)
  }
}

# 3. Standardization (center and scale)
# Calculate mean and sd from reference population
snp_means <- colMeans(geno_ref_filtered)
snp_sds <- apply(geno_ref_filtered, 2, sd)

# Standardize both populations
geno_ref_std <- scale(geno_ref_filtered, center = snp_means, scale = snp_sds)
geno_test_std <- scale(geno_test_filtered, center = snp_means, scale = snp_sds)
geno_test_gen17_std <- scale(geno_test_gen17_filtered, center = snp_means, scale = snp_sds)
geno_test_gen18_std <- scale(geno_test_gen18_filtered, center = snp_means, scale = snp_sds)

# Prepare phenotype vectors
y_ref <- phenotype_ref$phenotype
y_test <- phenotype_test$phenotype
tbv_test <- phenotype_test$tbv
y_test_gen17 <- phenotype_test_gen17$phenotype
tbv_test_gen17 <- phenotype_test_gen17$tbv
y_test_gen18 <- phenotype_test_gen18$phenotype
tbv_test_gen18 <- phenotype_test_gen18$tbv

cat("Preprocessing complete:\n")
cat("SNPs retained:", sum(keep_snps), "out of", ncol(geno_ref), "\n")
cat("Reference:", nrow(geno_ref_std), "individuals\n")
cat("Test:", nrow(geno_test_std), "individuals\n")

## Feature Engineering

#1.  Calculate the kinship/genomic relationship matrix (GRM) for GBLUP
#2.  Standardize the genotype matrix for ML models
#3.  CV fold for ML models

# Genomic Relationship Matrix (GRM) for GBLUP
# Using VanRaden method: G = ZZ'/k where k = 2*sum(p(1-p))
# Recalculate from geno_ref_filtered
allele_freq_filtered <- colMeans(geno_ref_filtered, na.rm = TRUE) / 2

Z <- geno_ref_filtered - matrix(2 * allele_freq_filtered, 
                                  nrow = nrow(geno_ref_filtered), 
                                  ncol = ncol(geno_ref_filtered), byrow = TRUE)
k <- 2 * sum(allele_freq_filtered * (1 - allele_freq_filtered))
GRM <- tcrossprod(Z) / k

# GRM evaluation
# Diagonal (self-relationship)
diag_mean <- mean(diag(GRM))
diag_range <- range(diag(GRM))

# Off-diagonal (relatedness)
offdiag <- GRM[lower.tri(GRM)]
offdiag_mean <- mean(offdiag)
offdiag_range <- range(offdiag)

# Create 5-fold CV indices for ML models
set.seed(123)
n_ref <- nrow(geno_ref_std)
fold_indices <- sample(rep(1:5, length.out = n_ref))

# Create CV fold list
cv_folds <- list()
for(fold in 1:5) {
  cv_folds[[fold]] <- list(
    train_idx = which(fold_indices != fold),
    val_idx = which(fold_indices == fold)
  )
}

cat("Feature engineering complete:\n")
cat("GRM dimensions:", dim(GRM), "\n")
cat("Diagonal mean:", diag_mean, "range:", diag_range, "\n")
cat("Off-diagonal mean:", offdiag_mean, "range:", offdiag_range, "\n")
cat("Standardized genotype:", dim(geno_ref_std), "\n")
cat("CV folds created: 5 folds\n")
cat("Fold sizes:", sapply(cv_folds, function(x) length(x$val_idx)), "\n")

# ============================================
# 3. MODEL TRAINING - STATISTICAL METHODS
# ============================================

# ============================================
# 3.A. GBLUP using rrBLUP
# ============================================

# Prepare data for GBLUP
gblup_time <- system.time({
  gblup_fit <- mixed.solve(y = y_ref, K = GRM)
})

# Extract results
gblup_intercept <- as.numeric(gblup_fit$beta)
gblup_gebv_ref <- gblup_fit$u
gblup_Vg <- gblup_fit$Vu
gblup_Ve <- gblup_fit$Ve
gblup_h2 <- gblup_Vg / (gblup_Vg + gblup_Ve)

cat("GBLUP completed - Time:", round(gblup_time[3], 2), "sec\n")
cat("h2:", round(gblup_h2, 3), "\n")

# ============================================
# 3.B. Bayesian Methods using HIBAYES
# ============================================

set.seed(123)

# Prepare data frame for HIBAYES
pheno_data <- data.frame(
  id = ref_ids,
  phenotype = y_ref
)

# BayesA
bayesA_time <- system.time({
fit_bayesA <- ibrm(
  formula = phenotype ~ 1,
  data = pheno_data,
  M = geno_ref_filtered,
  M.id = ref_ids,
  method = "BayesA",
  niter = 30000,
  nburn = 10000,
  thin = 5,
  threads = 4,
  verbose = FALSE
)
})

cat("BayesA completed - Time:", round(bayesA_time[3], 2), "sec\n")
cat("h2:", round(fit_bayesA$h2, 3), "\n")

# BayesR
bayesR_time <- system.time({
fit_bayesR <- ibrm(
  formula = phenotype ~ 1,
  data = pheno_data,
  M = geno_ref_filtered,
  M.id = ref_ids,
  method = "BayesR",
  fold = c(0, 0.001, 0.01, 0.1),  
  vg = gblup_Vg,                   # use Vg from GBLUP
  ve = gblup_Ve,                   # use Ve from GBLUP
  niter = 30000,
  nburn = 10000,
  thin = 5,
  threads = 4,
  verbose = FALSE
)
})

cat("BayesR completed - Time:", round(bayesR_time[3], 2), "sec\n")
cat("h2:", round(fit_bayesR$h2, 3), "\n")

# ============================================
# 3.B. MACHINE LEARNING MODELS WITH TUNING
# ============================================

# Clean data
na_cols <- colSums(is.na(geno_ref_std)) > 0
geno_ref_std_clean <- geno_ref_std[, !na_cols]
geno_test_std_clean <- geno_test_std[, !na_cols]
geno_test_gen17_std_clean <- geno_test_gen17_std[, !na_cols]
geno_test_gen18_std_clean <- geno_test_gen18_std[, !na_cols]

cat("Removed", sum(na_cols), "SNPs with zero variance\n")
cat("Clean matrix:", ncol(geno_ref_std_clean), "SNPs\n\n")

# ============================================
# 3.B.1. Random Forest with CV tuning
# ============================================
if(FALSE){ # If you want to do tuning via CV, just remove the if(FALSE) condition
rf_grid <- expand.grid(
  ntree = c(300, 500, 700),
  mtry = c(50, 100, 200),
  nodesize = c(5, 10)
)

rf_cv_results <- matrix(NA, nrow = nrow(rf_grid), ncol = 5)

rf_tuning_time <- system.time({
for(i in seq_len(nrow(rf_grid))) {
  for(fold in 1:5) {
    train_idx <- cv_folds[[fold]]$train_idx
    val_idx <- cv_folds[[fold]]$val_idx
    
    rf_temp <- randomForest(
      x = geno_ref_std_clean[train_idx, ],
      y = y_ref[train_idx],
      xtest = geno_ref_std_clean[val_idx, ],
      ytest = y_ref[val_idx],
      ntree = rf_grid$ntree[i],
      mtry = rf_grid$mtry[i],
      nodesize = rf_grid$nodesize[i]
    )
    
    rf_cv_results[i, fold] <- rf_temp$test$mse[rf_grid$ntree[i]]
  }
}
})

rf_cv_mean <- rowMeans(rf_cv_results)
best_rf <- which.min(rf_cv_mean)

cat("RF tuning time:", round(rf_tuning_time[3], 2), "sec\n")
cat("RF best params:\n")
cat("  ntree =", rf_grid$ntree[best_rf], "\n")
cat("  mtry =", rf_grid$mtry[best_rf], "\n")
cat("  nodesize =", rf_grid$nodesize[best_rf], "\n")
}

rf_train_time <- system.time({
set.seed(123)
rf_model <- randomForest(
  x = geno_ref_std_clean,
  y = y_ref,
  ntree = 500, #rf_grid$ntree[best_rf], # un-comment this for tuning via CV
  mtry = 200, #rf_grid$mtry[best_rf], # # un-comment this for tuning via CV
  nodesize = 10, #rf_grid$nodesize[best_rf], # un-comment this for tuning via CV
  importance = TRUE
)
})

cat("RF training complete with time of:", round(rf_train_time[3], 2), "sec\n\n")

# ============================================
# 3.B.2. XGBoost with CV tuning
# ============================================
if(FALSE){
xgb_grid <- expand.grid(
  eta = c(0.01, 0.1),
  max_depth = c(4, 6),
  subsample = c(0.8, 1.0),
  colsample_bytree = c(0.6, 0.8),
  lambda = c(0, 0.1, 1, 10)
)

xgb_cv_results <- matrix(NA, nrow = nrow(xgb_grid), ncol = 5)

xgb_tuning_time <- system.time({
for(i in seq_len(nrow(xgb_grid))) {
  for(fold in 1:5) {
    train_idx <- cv_folds[[fold]]$train_idx
    val_idx <- cv_folds[[fold]]$val_idx
    
    dtrain <- xgb.DMatrix(data = geno_ref_std_clean[train_idx, ], 
                          label = y_ref[train_idx])
    dval <- xgb.DMatrix(data = geno_ref_std_clean[val_idx, ], 
                        label = y_ref[val_idx])
    
    xgb_temp <- xgb.train(
      params = list(
        objective = "reg:squarederror",
        eta = xgb_grid$eta[i],
        max_depth = xgb_grid$max_depth[i],
        subsample = xgb_grid$subsample[i],
        colsample_bytree = xgb_grid$colsample_bytree[i],
        lambda = xgb_grid$lambda[i]
      ),
      data = dtrain,
      nrounds = 200,
      verbose = 0
    )
    
    pred <- predict(xgb_temp, dval)
    xgb_cv_results[i, fold] <- mean((pred - y_ref[val_idx])^2)
  }
}
})

xgb_cv_mean <- rowMeans(xgb_cv_results)
best_xgb <- which.min(xgb_cv_mean)

cat("XGB tuning time:", round(xgb_tuning_time[3], 2), "sec\n")
cat("XGB best params:\n")
cat("  eta =", xgb_grid$eta[best_xgb], "\n")
cat("  max_depth =", xgb_grid$max_depth[best_xgb], "\n")
cat("  subsample =", xgb_grid$subsample[best_xgb], "\n")
cat("  colsample_bytree =", xgb_grid$colsample_bytree[best_xgb], "\n")
cat("  lambda =", xgb_grid$lambda[best_xgb], "\n")
}

xgb_train_time <- system.time({
xgb_train <- xgb.DMatrix(data = geno_ref_std_clean, label = y_ref)
xgb_model <- xgb.train(
  params = list(
    objective = "reg:squarederror",
    eta = 0.1, #xgb_grid$eta[best_xgb], # un-comment this for tuning via CV
    max_depth = 4, #xgb_grid$max_depth[best_xgb], # un-comment this for tuning via CV
    subsample = 1, #xgb_grid$subsample[best_xgb], # un-comment this for tuning via CV
    colsample_bytree = 0.6, #xgb_grid$colsample_bytree[best_xgb], # un-comment this for tuning via CV
    lambda = 1 #xgb_grid$lambda[best_xgb] # un-comment this for tuning via CV
  ),
  data = xgb_train,
  nrounds = 200,
  verbose = 0
)
})

cat("XGB training time:", round(xgb_train_time[3], 2), "sec\n\n")


# ============================================
# 3.B.3. SVR with CV tuning
# ============================================
if(FALSE){
svr_grid <- expand.grid(
  cost = c(0.1, 1, 10), #0.001 - 100
  epsilon = c(0.01, 0.1)
)

svr_cv_results <- matrix(NA, nrow = nrow(svr_grid), ncol = 5)

svr_tuning_time <- system.time({
for(i in seq_len(nrow(svr_grid))) {
  for(fold in 1:5) {
    train_idx <- cv_folds[[fold]]$train_idx
    val_idx <- cv_folds[[fold]]$val_idx
    
    svr_temp <- svm(
      x = geno_ref_std_clean[train_idx, ],
      y = y_ref[train_idx],
      kernel = "radial",
      cost = svr_grid$cost[i],
      epsilon = svr_grid$epsilon[i]
    )
    
    pred <- predict(svr_temp, geno_ref_std_clean[val_idx, ])
    svr_cv_results[i, fold] <- mean((pred - y_ref[val_idx])^2)
  }
}
})

svr_cv_mean <- rowMeans(svr_cv_results)
best_svr <- which.min(svr_cv_mean)

cat("SVR tuning time:", round(svr_tuning_time[3], 2), "sec\n")
cat("SVR best params: cost =", svr_grid$cost[best_svr],
    "epsilon =", svr_grid$epsilon[best_svr], "\n")
}  
#gamma

svr_train_time <- system.time({
svr_model <- svm(
  x = geno_ref_std_clean,
  y = y_ref,
  kernel = "radial", #linear, polynomial, sigmoid
  cost = 1, # svr_grid$cost[best_svr], # un-comment this for tuning via CV
  epsilon = 0.1 # svr_grid$epsilon[best_svr] # un-comment this for tuning via CV
)
})

cat("SVR training time:", round(svr_train_time[3], 2), "sec\n\n")

# ============================================
# 4. PREDICTION ON TEST POPULATION
# ============================================

# ============================================
# 4.A. GBLUP Prediction
# ============================================

# Compute cross-relationship between test and reference
allele_freq_filtered <- colMeans(geno_ref_filtered, na.rm = TRUE) / 2
Z_ref <- geno_ref_filtered - matrix(2 * allele_freq_filtered, 
                                     nrow = nrow(geno_ref_filtered), 
                                     ncol = ncol(geno_ref_filtered), 
                                     byrow = TRUE)
Z_test <- geno_test_filtered - matrix(2 * allele_freq_filtered, 
                                       nrow = nrow(geno_test_filtered), 
                                       ncol = ncol(geno_test_filtered), 
                                       byrow = TRUE)
k <- 2 * sum(allele_freq_filtered * (1 - allele_freq_filtered))
G_test_ref <- tcrossprod(Z_test, Z_ref) / k

lambda_reg <- 1e-6  # Regularization parameter
n_ref <- nrow(GRM)
n_test <- nrow(Z_test)
stopifnot(ncol(G_test_ref) == n_ref)
GRM_reg <- GRM + diag(lambda_reg, nrow(GRM))
gblup_pred_time <- system.time({
  pred_gblup <- gblup_intercept + G_test_ref %*% solve(GRM_reg) %*% (y_ref - gblup_intercept)
})

# GBLUP Gen 17
Z_test_gen17 <- geno_test_gen17_filtered - matrix(2 * allele_freq_filtered, 
                                                    nrow = nrow(geno_test_gen17_filtered), 
                                                    ncol = ncol(geno_test_gen17_filtered), 
                                                    byrow = TRUE)
G_test_gen17_ref <- tcrossprod(Z_test_gen17, Z_ref) / k
pred_gblup_gen17 <- gblup_intercept + G_test_gen17_ref %*% solve(GRM_reg) %*% (y_ref - gblup_intercept)

# GBLUP Gen 18
Z_test_gen18 <- geno_test_gen18_filtered - matrix(2 * allele_freq_filtered, 
                                                    nrow = nrow(geno_test_gen18_filtered), 
                                                    ncol = ncol(geno_test_gen18_filtered), 
                                                    byrow = TRUE)
G_test_gen18_ref <- tcrossprod(Z_test_gen18, Z_ref) / k
pred_gblup_gen18 <- gblup_intercept + G_test_gen18_ref %*% solve(GRM_reg) %*% (y_ref - gblup_intercept)

cat("GBLUP prediction time:", round(gblup_pred_time[3], 2), "sec\n")

# ============================================
# 4.B. Bayesian Methods Prediction
# ============================================

bayesA_pred_time <- system.time({
  pred_bayesA <- fit_bayesA$mu + geno_test_filtered %*% fit_bayesA$alpha
})

bayesR_pred_time <- system.time({
  pred_bayesR <- fit_bayesR$mu + geno_test_filtered %*% fit_bayesR$alpha
})

# Bayesian Methods Gen 12
pred_bayesA_gen17 <- fit_bayesA$mu + geno_test_gen17_filtered %*% fit_bayesA$alpha
pred_bayesR_gen17 <- fit_bayesR$mu + geno_test_gen17_filtered %*% fit_bayesR$alpha

# Bayesian Methods Gen 13
pred_bayesA_gen18 <- fit_bayesA$mu + geno_test_gen18_filtered %*% fit_bayesA$alpha
pred_bayesR_gen18 <- fit_bayesR$mu + geno_test_gen18_filtered %*% fit_bayesR$alpha

cat("BayesA prediction time:", round(bayesA_pred_time[3], 2), "sec\n")
cat("BayesR prediction time:", round(bayesR_pred_time[3], 2), "sec\n")

# ============================================
# 4.C. Machine Learning Predictions
# ============================================

rf_pred_time <- system.time({
  pred_rf <- predict(rf_model, geno_test_std_clean)
})

xgb_pred_time <- system.time({
  dtest <- xgb.DMatrix(data = geno_test_std_clean) #
  pred_xgb <- predict(xgb_model, dtest) #
})

svr_pred_time <- system.time({
  pred_svr <- predict(svr_model, geno_test_std_clean)
})

# ML Methods Gen 17
pred_rf_gen17 <- predict(rf_model, geno_test_gen17_std_clean)
dtest_gen17 <- xgb.DMatrix(data = geno_test_gen17_std_clean)
pred_xgb_gen17 <- predict(xgb_model, dtest_gen17)
pred_svr_gen17 <- predict(svr_model, geno_test_gen17_std_clean)

# ML Methods Gen 18
pred_rf_gen18 <- predict(rf_model, geno_test_gen18_std_clean)
dtest_gen18 <- xgb.DMatrix(data = geno_test_gen18_std_clean)
pred_xgb_gen18 <- predict(xgb_model, dtest_gen18)
pred_svr_gen18 <- predict(svr_model, geno_test_gen18_std_clean)

cat("RF prediction time:", round(rf_pred_time[3], 2), "sec\n")
cat("XGB prediction time:", round(xgb_pred_time[3], 2), "sec\n")
cat("SVR prediction time:", round(svr_pred_time[3], 2), "sec\n")

# ============================================
# 5. MODEL EVALUATION
# ============================================

# Function to calculate evaluation metrics
evaluate_model <- function(predicted, true_tbv, model_name) {
  # Correlation
  r <- cor(predicted, true_tbv, use = "complete.obs")
  
  # R-squared
  r2 <- r^2
  
  # RMSE
  rmse <- sqrt(mean((predicted - true_tbv)^2, , na.rm = TRUE))
  
  # MAE
  mae <- mean(abs(predicted - true_tbv), , na.rm = TRUE)
  
  # Bias (regression slope)
  bias_lm <- lm(predicted ~ true_tbv)
  bias_slope <- coef(bias_lm)[2]
  bias_intercept <- coef(bias_lm)[1]
  
  # Ranking correlation (Spearman)
  rank_cor <- cor(predicted, true_tbv, method = "spearman", use = "complete.obs")
  
  return(data.frame(
    Model = model_name,
    Correlation = r,
    R2 = r2,
    RMSE = rmse,
    MAE = mae,
    Bias_Slope = bias_slope,
    Bias_Intercept = bias_intercept,
    Rank_Correlation = rank_cor
  ))
}

# Evaluate Gen 11
results_gen16 <- rbind(
  evaluate_model(as.vector(pred_gblup), tbv_test, "GBLUP"),
  evaluate_model(as.vector(pred_bayesA), tbv_test, "BayesA"),
  evaluate_model(as.vector(pred_bayesR), tbv_test, "BayesR"),
  evaluate_model(pred_rf, tbv_test, "Random_Forest"),
  evaluate_model(pred_xgb, tbv_test, "XGBoost"),
  evaluate_model(pred_svr, tbv_test, "SVR")
)

# Evaluate Gen 12
results_gen17 <- rbind(
  evaluate_model(as.vector(pred_gblup_gen17), tbv_test_gen17, "GBLUP"),
  evaluate_model(as.vector(pred_bayesA_gen17), tbv_test_gen17, "BayesA"),
  evaluate_model(as.vector(pred_bayesR_gen17), tbv_test_gen17, "BayesR"),
  evaluate_model(pred_rf_gen17, tbv_test_gen17, "Random_Forest"),
  evaluate_model(pred_xgb_gen17, tbv_test_gen17, "XGBoost"),
  evaluate_model(pred_svr_gen17, tbv_test_gen17, "SVR")
)

# Evaluate Gen 13
results_gen18 <- rbind(
  evaluate_model(as.vector(pred_gblup_gen18), tbv_test_gen18, "GBLUP"),
  evaluate_model(as.vector(pred_bayesA_gen18), tbv_test_gen18, "BayesA"),
  evaluate_model(as.vector(pred_bayesR_gen18), tbv_test_gen18, "BayesR"),
  evaluate_model(pred_rf_gen18, tbv_test_gen18, "Random_Forest"),
  evaluate_model(pred_xgb_gen18, tbv_test_gen18, "XGBoost"),
  evaluate_model(pred_svr_gen18, tbv_test_gen18, "SVR")
)

# Combine results
results_gen16$Generation <- "Gen_16"
results_gen17$Generation <- "Gen_17"
results_gen18$Generation <- "Gen_18"
results <- rbind(results_gen16, results_gen17, results_gen18)

# Add computational efficiency metrics
results$Train_Time <- rep(c(
  gblup_time[3],
  bayesA_time[3],
  bayesR_time[3],
  rf_train_time[3],
  xgb_train_time[3],
  svr_train_time[3]
), 3)

results$Pred_Time <- rep(c(
  gblup_pred_time[3],
  bayesA_pred_time[3],
  bayesR_pred_time[3],
  rf_pred_time[3],
  xgb_pred_time[3],
  svr_pred_time[3]
), 3)

# Display results
cat("\n=== MODEL COMPARISON RESULTS ===\n\n")
print(results, row.names = FALSE, digits = 4)

# Comparison table
comparison <- data.frame(
  Model = results_gen16$Model,
  Corr_Gen16 = results_gen16$Correlation,
  Corr_Gen17 = results_gen17$Correlation,
  Corr_Gen18 = results_gen18$Correlation,
  Diff = results_gen18$Correlation - results_gen16$Correlation,
  Percent_Decline = 100 * (results_gen18$Correlation - results_gen16$Correlation) / results_gen16$Correlation
)

# Rank models by correlation
results_ranked <- results[order(-results$Correlation), ]

cat("\n=== MODELS RANKED BY PREDICTION ACCURACY ===\n\n")
print(results_ranked[, c("Model", "Correlation", "R2", "RMSE")], row.names = FALSE, digits = 4)

cat("\n=== ACCURACY COMPARISON: GEN 11 vs GEN 12 vs GEN 13 ===\n\n")
print(comparison, row.names = FALSE, digits = 4)

# Save performance metrics
write.csv(results, 
          file.path(output_dir, "model_performance_metrics.csv"), 
          row.names = FALSE)

# Save ranked results
write.csv(results_ranked, 
          file.path(output_dir, "models_ranked_by_accuracy.csv"), 
          row.names = FALSE)

write.csv(comparison, 
          file.path(output_dir, "comparison_among_generations.csv"), 
          row.names = FALSE)

cat("Saved: model_performance_metrics.csv\n")
cat("Saved: models_ranked_by_accuracy.csv\n\n")

# ============================================
# 6. MODEL COMPARISON & VISUALIZATION
# ============================================

# ============================================
# 6.1. Extract CV Performance for Boxplot
# ============================================

if(FALSE){

# Convert CV results to correlation per fold
cv_performance <- data.frame(
  Model = character(),
  Fold = integer(),
  Correlation = numeric()
)

# Random Forest CV
for(fold in 1:5) {
  val_idx <- cv_folds[[fold]]$val_idx
  train_idx <- cv_folds[[fold]]$train_idx
  
  rf_temp <- randomForest(
    x = geno_ref_std_clean[train_idx, ],
    y = y_ref[train_idx],
    ntree = 500, #rf_grid$ntree[best_rf],
    mtry = 200, #rf_grid$mtry[best_rf],
    nodesize = 10 #rf_grid$nodesize[best_rf]
  )
  pred <- predict(rf_temp, geno_ref_std_clean[val_idx, ])
  cv_performance <- rbind(cv_performance, 
                          data.frame(Model = "Random_Forest", 
                                   Fold = fold, 
                                   Correlation = cor(pred, y_ref[val_idx])))
}

# XGBoost CV
for(fold in 1:5) {
  val_idx <- cv_folds[[fold]]$val_idx
  train_idx <- cv_folds[[fold]]$train_idx
  
  dtrain <- xgb.DMatrix(data = geno_ref_std_clean[train_idx, ], 
                        label = y_ref[train_idx])
  dval <- xgb.DMatrix(data = geno_ref_std_clean[val_idx, ])
  
  xgb_temp <- xgb.train(
    params = list(
      objective = "reg:squarederror",
      eta = 0.1, #xgb_grid$eta[best_xgb],
      max_depth = 4, #xgb_grid$max_depth[best_xgb],
      subsample = 1, #xgb_grid$subsample[best_xgb],
      colsample_bytree = 0.6, #xgb_grid$colsample_bytree[best_xgb],
      lambda = 1 #xgb_grid$lambda[best_xgb]
    ),
    data = dtrain,
    nrounds = 200,
    verbose = 0
  )
  pred <- predict(xgb_temp, dval)
  cv_performance <- rbind(cv_performance, 
                          data.frame(Model = "XGBoost", 
                                   Fold = fold, 
                                   Correlation = cor(pred, y_ref[val_idx])))
}

# ============================================
# 6.2. Boxplot of CV Performance
# ============================================

p_boxplot <- ggplot(cv_performance, aes(x = Model, y = Correlation, fill = Model)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, alpha = 0.5) +
  theme_minimal() +
  labs(title = "Cross-Validation Performance Across ML Models",
       y = "Correlation (5-fold CV)",
       x = "Model") +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(output_dir, "cv_performance_boxplot.png"), 
       plot = p_boxplot,
       width = 10, height = 6, dpi = 300, bg = "white")

cat("Saved: cv_performance_boxplot.png\n")

write.csv(cv_performance, 
          file.path(output_dir, "cv_performance_5fold.csv"), 
          row.names = FALSE)

cat("Saved: cv_performance_5fold.csv\n\n")

}

# ============================================
# 6.3. Scatter Plots: Predicted vs True TBV
# ============================================

scatter_data <- data.frame(
  True_TBV = c(rep(tbv_test, 6), rep(tbv_test_gen17, 6), rep(tbv_test_gen18, 6)),
  Predicted = c(as.vector(pred_gblup), as.vector(pred_bayesA), 
                as.vector(pred_bayesR), pred_rf, pred_xgb, pred_svr,
                as.vector(pred_gblup_gen17), as.vector(pred_bayesA_gen17),
                as.vector(pred_bayesR_gen17), pred_rf_gen17, pred_xgb_gen17,
                pred_svr_gen17, as.vector(pred_gblup_gen18), as.vector(pred_bayesA_gen18),
                as.vector(pred_bayesR_gen18), pred_rf_gen18, pred_xgb_gen18,
                pred_svr_gen18),
  Model = rep(c("GBLUP", "BayesA", "BayesR", "Random_Forest", "XGBoost", "SVR"), 3),
  Generation = rep(c("Gen_16", "Gen_17", "Gen_18"), each = length(tbv_test) * 6)
)

p_scatter <- ggplot(scatter_data, aes(x = True_TBV, y = Predicted, color = Generation)) +
  geom_point(alpha = 0.5, size = 1.5) +
  geom_abline(intercept = 0, slope = 1, color = "black", linetype = "dashed") +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.8) +
  facet_wrap(~ Model, scales = "free", ncol = 3) +
  theme_minimal() +
  labs(title = "Predicted vs true breeding values",
       x = "True breeding value (TBV)",
       y = "Predicted breeding value") +
  theme(strip.text = element_text(face = "bold"),
        legend.position = "bottom") +
  scale_color_manual(values = c("Gen_16" = "#2E86AB", "Gen_17" = "#A23B72", "Gen_18" = "#8e341f"))

ggsave(file.path(output_dir, "predicted_vs_true_scatter_both_gen.png"), 
       plot = p_scatter,
       width = 18, height = 8, dpi = 300, bg = "white")

cat("Saved: predicted_vs_true_scatter_both_gen.png\n")

# ============================================
# 6.4. Performance Comparison Barplot
# ============================================

p_corr <- ggplot(results, aes(x = Model, y = Correlation, fill = Generation)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() +
  labs(title = "Model prediction accuracy",
       x = "Model",
       y = "Correlation with TBV") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom") +
  scale_fill_manual(values = c("Gen_16" = "#2E86AB", "Gen_17" = "#A23B72", "Gen_18" = "#8e341f")) +
  geom_text(aes(label = round(Correlation, 3)), 
            position = position_dodge(width = 0.9), 
            vjust = -0.5, size = 3)

p_rmse <- ggplot(results, aes(x = Model, y = RMSE, fill = Generation)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() +
  labs(title = "Model prediction error",
       x = "Model",
       y = "RMSE") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom") +
  scale_fill_manual(values = c("Gen_16" = "#2E86AB", "Gen_17" = "#A23B72", "Gen_18" = "#8e341f")) +
  geom_text(aes(label = round(RMSE, 3)), 
            position = position_dodge(width = 0.9), 
            vjust = -0.5, size = 3)

# Combine plots
p_combined <- grid.arrange(p_corr, p_rmse, ncol = 2)

# Save to file
ggsave(
  filename = file.path(output_dir, "performance_comparison_barplots.png"),
  plot = p_combined,
  width = 18,
  height = 6,
  dpi = 300,
  bg = "white"
)

cat("Saved: performance_comparison_barplots.png\n")

# ============================================
# 6.5. Accuracy Decline Visualization
# ============================================

p_decline <- ggplot(comparison, aes(x = reorder(Model, -Percent_Decline), 
                                     y = Percent_Decline, fill = Model)) +
  geom_bar(stat = "identity") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  coord_flip() +
  theme_minimal() +
  labs(title = "Prediction accuracy trend from gen 11 to gen 13",
       x = "Model",
       y = "% Change in correlation") +
  theme(legend.position = "none") +
  geom_text(aes(label = paste0(round(Percent_Decline, 1), "%")), hjust = 1.2)

ggsave(file.path(output_dir, "accuracy_decline_gen16_to_gen18.png"), 
       plot = p_decline,
       width = 18, height = 6, dpi = 300, bg = "white")

cat("Saved: accuracy_decline_gen16_to_gen18.png\n")

# ============================================
# 6.6. Computational Efficiency Comparison
# ============================================
train_time_data <- results[results$Generation == "Gen_16", c("Model", "Train_Time")]

p_time <- ggplot(train_time_data, aes(x = reorder(Model, Train_Time), y = Train_Time, fill = Model)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  theme_minimal() +
  labs(title = "Training time comparison",
       x = "Model",
       y = "Training time (seconds)") +
  theme(legend.position = "none") +
  geom_text(aes(label = round(Train_Time, 2)), hjust = -0.1)

ggsave(file.path(output_dir, "training_time_comparison.png"), 
       plot = p_time,
       width = 15, height = 6, dpi = 300, bg = "white")

cat("Saved: training_time_comparison.png\n\n")

# ============================================
# 7. OUTPUT RESULTS
# ============================================

# ============================================
# 7.1. Save Prediction Results
# ============================================

# Gen 11 predictions
predictions_gen16 <- data.frame(
  individual_id = test_ids,
  True_TBV = tbv_test,
  True_Phenotype = y_test,
  GBLUP = as.vector(pred_gblup),
  BayesA = as.vector(pred_bayesA),
  BayesR = as.vector(pred_bayesR),
  Random_Forest = pred_rf,
  XGBoost = pred_xgb,
  SVR = pred_svr,
  Generation = paste0("Gen_", test_gens[1])
)

# Gen 12 predictions
predictions_gen17 <- data.frame(
  individual_id = test_ids_gen17,
  True_TBV = tbv_test_gen17,
  True_Phenotype = y_test_gen17,
  GBLUP = as.vector(pred_gblup_gen17),
  BayesA = as.vector(pred_bayesA_gen17),
  BayesR = as.vector(pred_bayesR_gen17),
  Random_Forest = pred_rf_gen17,
  XGBoost = pred_xgb_gen17,
  SVR = pred_svr_gen17,
  Generation = paste0("Gen_", test_gens[2])
)

# Gen 13 predictions
predictions_gen18 <- data.frame(
  individual_id = test_ids_gen18,
  True_TBV = tbv_test_gen18,
  True_Phenotype = y_test_gen18,
  GBLUP = as.vector(pred_gblup_gen18),
  BayesA = as.vector(pred_bayesA_gen18),
  BayesR = as.vector(pred_bayesR_gen18),
  Random_Forest = pred_rf_gen18,
  XGBoost = pred_xgb_gen18,
  SVR = pred_svr_gen18,
  Generation = paste0("Gen_", test_gens[3])
)

# Save 
write.csv(predictions_gen16, 
          file.path(output_dir, "predictions_gen16.csv"), 
          row.names = FALSE)

write.csv(predictions_gen17, 
          file.path(output_dir, "predictions_gen17.csv"), 
          row.names = FALSE)

write.csv(predictions_gen18, 
          file.path(output_dir, "predictions_gen18.csv"), 
          row.names = FALSE)

cat("\n=== PREDICTION RESULTS", paste0("Gen ", test_gens[1]), "(First 10) ===\n")
print(head(predictions_gen16, 10))

cat("\n=== PREDICTION RESULTS", paste0("Gen ", test_gens[2]), "(First 10) ===\n")
print(head(predictions_gen17, 10))

cat("\n=== PREDICTION RESULTS", paste0("Gen ", test_gens[3]), "(First 10) ===\n")
print(head(predictions_gen18, 10))

cat("Saved: predictions_gen16.csv\n")
cat("Saved: predictions_gen17.csv\n")

# ============================================
# 7.2. Feature Importance (RF & XGBoost)
# ============================================

# Random Forest Variable Importance
rf_importance <- importance(rf_model)
rf_imp_df <- data.frame(
  SNP = rownames(rf_importance),
  IncMSE = rf_importance[, "%IncMSE"],
  IncNodePurity = rf_importance[, "IncNodePurity"]
)
rf_imp_top <- rf_imp_df[order(-rf_imp_df$IncMSE), ][1:20, ]

write.csv(rf_imp_top, 
          file.path(output_dir, "rf_top20_important_snps.csv"), 
          row.names = FALSE)

cat("\n=== TOP 20 IMPORTANT SNPs (Random Forest) ===\n")
print(rf_imp_top, row.names = FALSE)
cat("Saved: rf_top20_important_snps.csv\n")

# XGBoost Feature Importance
xgb_importance <- xgb.importance(model = xgb_model)
xgb_imp_top <- head(xgb_importance, 20)

write.csv(xgb_imp_top, 
          file.path(output_dir, "xgb_top20_important_features.csv"), 
          row.names = FALSE)

cat("\n=== TOP 20 IMPORTANT FEATURES (XGBoost) ===\n")
print(xgb_imp_top, row.names = FALSE)
cat("Saved: xgb_top20_important_features.csv\n")

# Plot importance
p_rf_imp <- ggplot(rf_imp_top, aes(x = reorder(SNP, IncMSE), y = IncMSE)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  theme_minimal() +
  labs(title = "Top 20 SNPs - Random Forest",
       x = "SNP",
       y = "% Increase in MSE")

p_xgb_imp <- ggplot(xgb_imp_top, aes(x = reorder(Feature, Gain), y = Gain)) +
  geom_bar(stat = "identity", fill = "coral") +
  coord_flip() +
  theme_minimal() +
  labs(title = "Top 20 Features - XGBoost",
       x = "Feature",
       y = "Gain")

# Combine plots
p_ml_importance_combined <- grid.arrange(p_rf_imp, p_xgb_imp, ncol = 2)

# Save to file
ggsave(
  filename = file.path(output_dir, "feature_importance_comparison.png"),
  plot = p_ml_importance_combined,
  width = 12,
  height = 6,
  dpi = 300,
  bg = "white"
)

cat("Saved: feature_importance_comparison.png\n")

# ============================================
# 7.3. SNP Effects (Bayesian Methods)
# ============================================

snp_effects <- data.frame(
  SNP = colnames(geno_ref_filtered),
  BayesA_Effect = fit_bayesA$alpha,
  BayesR_Effect = fit_bayesR$alpha,
  BayesA_PIP = fit_bayesA$pip,
  BayesR_PIP = fit_bayesR$pip
)

# Top SNPs by absolute effect size
snp_effects$BayesA_AbsEffect <- abs(snp_effects$BayesA_Effect)
snp_top_bayesA <- snp_effects[order(-snp_effects$BayesA_AbsEffect), ][1:20, ]

write.csv(snp_effects, 
          file.path(output_dir, "bayesian_snp_effects_all.csv"), 
          row.names = FALSE)

write.csv(snp_top_bayesA, 
          file.path(output_dir, "bayesA_top20_snps.csv"), 
          row.names = FALSE)

cat("\n=== TOP 20 SNPs BY EFFECT SIZE (BayesA) ===\n")
print(snp_top_bayesA[, c("SNP", "BayesA_Effect", "BayesA_PIP")], row.names = FALSE)
cat("Saved: bayesian_snp_effects_all.csv\n")
cat("Saved: bayesA_top20_snps.csv\n")

# Plot Manhattan-style SNP effects
p_snp_effects <- ggplot(snp_effects, aes(x = seq_along(SNP), y = BayesA_Effect)) +
  geom_point(aes(color = BayesA_PIP), alpha = 0.6, size = 1) +
  scale_color_gradient(low = "grey", high = "red") +
  theme_minimal() +
  labs(title = "SNP Effects - BayesA",
       x = "SNP Index",
       y = "Effect Size",
       color = "PIP") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black")

ggsave(file.path(output_dir, "bayesA_snp_effects_manhattan.png"), 
       plot = p_snp_effects,
       width = 10, height = 6, dpi = 300, bg = "white")

cat("Saved: bayesA_snp_effects_manhattan.png\n\n")

# ============================================
# 7.4. Summary Report
# ============================================

report_output <- capture.output({
  cat(strrep("=", 60), "\n")
  cat("         GENOMIC PREDICTION MODEL COMPARISON REPORT\n")
  cat(strrep("=", 60), "\n\n")
  
  cat("Dataset Information:\n")
  cat("  Reference Population:", nrow(geno_ref_std), "individuals\n")
  cat("  Test Population:", nrow(geno_test_std), "individuals\n")
  cat("  SNPs after QC: ", ncol(geno_ref_std_clean), "\n")
  cat("  True heritability (GBLUP): ", round(gblup_h2, 3), "\n\n")
  
  cat("Model Performance Summary:\n")
  cat("  Best Model (Correlation): ", results_ranked$Model[1], 
      " (r = ", round(results_ranked$Correlation[1], 3), ")\n", sep = "")
  cat("  Best Model (RMSE): ", results_gen16[which.min(results_gen16$RMSE), "Model"], 
      " (RMSE = ", round(min(results_gen16$RMSE), 3), ")\n", sep = "")
  cat("  Fastest Training: ", results_gen16[which.min(results_gen16$Train_Time), "Model"], 
      " (", round(min(results_gen16$Train_Time), 2), " sec)\n", sep = "")
  cat("  Fastest Prediction: ", results_gen16[which.min(results_gen16$Pred_Time), "Model"], 
      " (", round(min(results_gen16$Pred_Time), 4), " sec)\n\n", sep = "")
  cat("\nAccuracy Decline Summary:\n")
  cat("  Average decline: ", round(mean(comparison$Percent_Decline), 2), "%\n")
  cat("  Most robust model: ", comparison$Model[which.min(abs(comparison$Percent_Decline))], "\n")
  
  cat("\nDetailed Model Performance:\n")
  print(results)
  
  cat("\n\nKey Findings:\n")
  cat("  - Statistical methods (GBLUP) vs ML methods comparison\n")
  cat("  - Bayesian methods struggled with sparse genetic architecture\n")
  cat("  - ML methods show competitive performance with proper tuning\n")
  cat("  - Computational efficiency varies significantly across methods\n\n")
  
  cat(strrep("=", 60), "\n")
})

# Write ke file
writeLines(report_output, file.path(output_dir, "summary_report.txt"))

# Progress message (muncul di HPC log)
cat("Saved: summary_report.txt")

# ============================================
# FINAL MESSAGE
# ============================================

cat(strrep("=", 60))
cat("ALL RESULTS SAVED SUCCESSFULLY!")
cat(strrep("=", 60))
cat("Output directory: ", output_dir)
cat("\nFiles saved:")
cat("  CSV files: 10")
cat("  PNG images: 6")
cat("  Text report: 1")
cat(strrep("=", 60))
