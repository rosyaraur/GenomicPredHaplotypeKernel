# Load necessary library
if (!require("BGLR")) install.packages("BGLR")
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("tidyr")) install.packages("tidyr")
library(ggplot2)
library(tidyr)
library(BGLR)


## =================================================================
## Module 1: Data Simulation
## =================================================================
simulate_genomic_data <- function(n_ind = 300, n_markers = 60000, 
                                  block_size = 3, n_qtl = 10, h2 = 0.4) {
  # 1. Simulate Marker Matrix (M)
  M <- matrix(sample(c(-1, 0, 1), size = n_ind * n_markers, replace = TRUE), 
              nrow = n_ind, ncol = n_markers)
  colnames(M) <- paste0("M", 1:n_markers)
  
  # 2. Construct Haplotype Matrix (H)
  n_blocks <- n_markers / block_size
  hap_list <- list()
  for (b in 1:n_blocks) {
    snps_in_block <- M[, ((b-1)*block_size + 1):(b*block_size)]
    hap_strings <- apply(snps_in_block, 1, paste, collapse = "_")
    hap_dummies <- model.matrix(~ as.factor(hap_strings) - 1)
    colnames(hap_dummies) <- paste0("B", b, "_", colnames(hap_dummies))
    hap_list[[b]] <- hap_dummies
  }
  H <- do.call(cbind, hap_list)
  
  # 3. Simulate Trait based on Haplotype QTLs
  qtl_idx <- sample(1:ncol(H), n_qtl)
  qtl_effects <- rnorm(n_qtl, mean = 0, sd = 4)
  TBV <- H[, qtl_idx] %*% qtl_effects
  
  var_g <- var(TBV)
  var_e <- var_g * (1 - h2) / h2
  Y <- as.numeric(TBV + rnorm(n_ind, mean = 0, sd = sqrt(var_e)))
  
  return(list(M = M, H = H, Y = Y))
}

## =================================================================
## Module 2: Utilities (Kernel Computation)
## =================================================================
compute_gaussian_kernel <- function(M) {
  D2 <- as.matrix(dist(M))^2
  theta <- 1 / median(D2)
  K <- exp(-theta * D2)
  return(K)
}

## =================================================================
## Module 3: Model Fitting (Standard)
## =================================================================
fit_standard_model <- function(Y_masked, matrix_input, model_type = "BRR", 
                               nIter = 2000, burnIn = 500) {
  # model_type: "BRR" for GBLUP, "RKHS" for Gaussian Kernel
  if (model_type == "RKHS") {
    ETA <- list(list(K = matrix_input, model = "RKHS"))
  } else {
    ETA <- list(list(X = matrix_input, model = model_type))
  }
  
  fit <- BGLR(y = Y_masked, ETA = ETA, 
              nIter = nIter, burnIn = burnIn, verbose = FALSE)
  return(fit$yHat)
}

## =================================================================
## Module 4: Model Fitting (GWAS-Assisted)
## =================================================================
fit_gwas_assisted_model <- function(Y, Y_masked, train_idx, fixed_matrix, 
                                    random_matrix = NULL, random_kernel = NULL, 
                                    n_fixed = 5, nIter = 2000, burnIn = 500) {
  
  # 1. Perform GWAS strictly on the training set to prevent data leakage
  p_vals <- numeric(ncol(fixed_matrix))
  for (i in 1:ncol(fixed_matrix)) {
    if(var(fixed_matrix[train_idx, i]) > 0) {
      fit <- lm(Y[train_idx] ~ fixed_matrix[train_idx, i])
      p_vals[i] <- summary(fit)$coefficients[2, 4]
    } else { 
      p_vals[i] <- 1 
    }
  }
  
  # 2. Identify top loci
  top_idx <- order(p_vals)[1:n_fixed]
  
  # 3. Build the BGLR ETA list
  ETA <- list(list(X = fixed_matrix[, top_idx, drop=FALSE], model = "FIXED"))
  
  # Add the random polygenic background (either Marker Matrix or Kernel)
  if (!is.null(random_kernel)) {
    ETA[[2]] <- list(K = random_kernel, model = "RKHS")
  } else if (!is.null(random_matrix)) {
    # Remove fixed markers from random matrix to avoid double-fitting
    ETA[[2]] <- list(X = random_matrix[, -top_idx], model = "BRR")
  }
  
  # 4. Fit the model
  fit <- BGLR(y = Y_masked, ETA = ETA, 
              nIter = nIter, burnIn = burnIn, verbose = FALSE)
  return(fit$yHat)
}

## =================================================================
## Updated Module 5: Cross-Validation Orchestrator
## =================================================================
run_modular_cv <- function(k_folds = 5, n_fixed = 5) {
  
  cat("1. Simulating Data & Computing Kernel...\n")
  data <- simulate_genomic_data()
  K <- compute_gaussian_kernel(data$M)
  
  folds <- sample(rep(1:k_folds, length.out = length(data$Y)))
  
  # Matrices to store results
  results <- matrix(NA, nrow = k_folds, ncol = 4)
  colnames(results) <- c("Standard_GBLUP", "Standard_RKHS", 
                         "GWAS_GBLUP", "GWAS_Hap_RKHS")
  
  # Vector to store the unbiased GEBVs for the best model
  best_model_gebvs <- numeric(length(data$Y))
  
  cat("2. Running", k_folds, "- Fold Cross Validation...\n")
  for (f in 1:k_folds) {
    cat(sprintf("   Processing Fold %d...\n", f))
    
    test_idx <- which(folds == f)
    train_idx <- which(folds != f)
    
    Y_masked <- data$Y
    Y_masked[test_idx] <- NA 
    
    # Fit Models (using functions from the previous script)
    pred_gblup <- fit_standard_model(Y_masked, data$M, "BRR")
    pred_rkhs  <- fit_standard_model(Y_masked, K, "RKHS")
    
    pred_gwas_gblup <- fit_gwas_assisted_model(
      data$Y, Y_masked, train_idx, data$M, random_matrix = data$M, n_fixed = n_fixed
    )
    
    pred_gwas_rkhs <- fit_gwas_assisted_model(
      data$Y, Y_masked, train_idx, data$H, random_kernel = K, n_fixed = n_fixed
    )
    
    # Store unbiased test-set GEBVs for the population
    best_model_gebvs[test_idx] <- pred_gwas_rkhs[test_idx]
    
    # Record Accuracies
    results[f, "Standard_GBLUP"] <- cor(pred_gblup[test_idx], data$Y[test_idx])
    results[f, "Standard_RKHS"]  <- cor(pred_rkhs[test_idx], data$Y[test_idx])
    results[f, "GWAS_GBLUP"]     <- cor(pred_gwas_gblup[test_idx], data$Y[test_idx])
    results[f, "GWAS_Hap_RKHS"]  <- cor(pred_gwas_rkhs[test_idx], data$Y[test_idx])
  }
  
  cat("\n=== Final Mean Accuracies ===\n")
  print(round(colMeans(results), 4))
  
  return(list(
    Accuracies = results, 
    True_Phenotypes = data$Y, 
    Predicted_GEBVs = best_model_gebvs
  ))
}

## =================================================================
## Module 6: Extract Breeding Values (GEBVs)
## =================================================================
#' Extract and rank the top performing lines based on GEBVs
#'
#' @param gebv_vector A vector of predicted phenotypes (GEBVs)
#' @param selection_intensity The percentage of top lines to select (e.g., 0.10 for top 10%)
#' @return A dataframe of the selected lines sorted by GEBV
extract_top_lines <- function(gebv_vector, selection_intensity = 0.10) {
  
  # Create a dataframe pairing individual IDs with their GEBVs
  n_ind <- length(gebv_vector)
  df <- data.frame(
    Line_ID = paste0("Line_", 1:n_ind),
    GEBV = gebv_vector
  )
  
  # Sort descending (assuming higher values are better, like Grain Yield)
  df_sorted <- df[order(-df$GEBV), ]
  
  # Calculate how many lines to keep
  n_select <- max(1, round(n_ind * selection_intensity))
  
  # Extract the top tier
  selected_lines <- head(df_sorted, n_select)
  
  cat(sprintf("\nSelected the top %d lines (%.0f%% intensity).\n", 
              n_select, selection_intensity * 100))
  
  return(selected_lines)
}


## =================================================================
## Module 7: Plot Results
## =================================================================
#' Plot cross-validation accuracies using ggplot2
#'
#' @param cv_matrix The matrix of accuracies output by the orchestrator
#' @return A ggplot object
plot_cv_results <- function(cv_matrix) {
  
  # Convert the matrix to a dataframe
  df <- as.data.frame(cv_matrix)
  df$Fold <- 1:nrow(df)
  
  # Reshape data from wide to long format for ggplot
  df_long <- pivot_longer(df, 
                          cols = -Fold, 
                          names_to = "Model", 
                          values_to = "Accuracy")
  
  # Reorder factor levels for better visual progression on the X-axis
  df_long$Model <- factor(df_long$Model, levels = c("Standard_GBLUP", 
                                                    "Standard_RKHS", 
                                                    "GWAS_GBLUP", 
                                                    "GWAS_Hap_RKHS"))
  
  # Create the boxplot
  p <- ggplot(df_long, aes(x = Model, y = Accuracy, fill = Model)) +
    geom_boxplot(alpha = 0.8, outlier.shape = 16, outlier.size = 2) +
    scale_fill_viridis_d() + # Uses a colorblind-friendly palette
    theme_minimal(base_size = 14) +
    labs(
      title = "Genomic Prediction Accuracies",
      subtitle = "Comparing Standard vs. GWAS-Assisted Models",
      x = "Prediction Model",
      y = "Pearson Correlation (Accuracy)"
    ) +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 15, hjust = 1, face = "bold"),
      panel.grid.major.x = element_blank()
    )
  
  print(p)
  return(p)
}

#####################
# --- Execution --- #
# 1. Run the pipeline
final_output <- run_modular_cv()

# 2. Visualize the cross-validation accuracies
accuracy_plot <- plot_cv_results(final_output$Accuracies)

# 3. Extract the top 10% of lines based on the best model's predictions
elite_lines <- extract_top_lines(final_output$Predicted_GEBVs, selection_intensity = 0.10)
print(head(elite_lines))
