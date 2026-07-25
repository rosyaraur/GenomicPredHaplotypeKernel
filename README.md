# GWAS-Assisted Genomic Prediction Pipeline

A modular R framework for incorporating Genome-Wide Association Study (GWAS) results into Genomic Selection (GS) models. 

This repository implements the proof-of-concept methodology outlined in [Sehgal et al. (2020) *Frontiers in Plant Science*](https://www.frontiersin.org/journals/plant-science/articles/10.3389/fpls.2020.00197/full). It demonstrates how modeling major loci as fixed effects—while capturing the polygenic background using non-linear Gaussian kernels (RKHS)—can significantly improve prediction accuracies for complex traits like grain yield.

## 🧬 Key Features

* **Two-Part Mixed Architecture:** Leverages the `BGLR` package to model top loci as unpenalized fixed effects while treating the remaining genome as a penalized random effect.
* **Haplotype Block Construction:** Groups adjacent SNPs into localized blocks to capture regional epistatic interactions.
* **Reproducing Kernel Hilbert Space (RKHS):** Uses a non-linear Gaussian kernel to capture non-additive genetic variance (cryptic epistasis) missed by standard GBLUP.
* **Strict Cross-Validation:** The pipeline strictly encapsulates the GWAS step within the training loop of each fold, ensuring zero data leakage between training and testing sets.
* **Highly Modular:** Cleanly separates data simulation, kernel computation, model fitting, and visualization so you can easily swap in your own empirical datasets.

## 📋 Prerequisites

This pipeline requires R (>= 4.0.0) and the following packages:

```R
install.packages(c("BGLR", "ggplot2", "tidyr"))

```

## 🚀 Quick Start

To run the full simulation, cross-validation, and selection pipeline, execute the orchestrator script.

```R
# 1. Load the scripts (assuming they are saved in your working directory)
source("1_data_simulation.R")
source("2_model_fitting.R")
source("3_orchestrator.R")

# 2. Run the 5-fold cross-validation pipeline (defaults to 5 fixed effects)
final_output <- run_modular_cv(k_folds = 5, n_fixed = 5)

# 3. Visualize the accuracies across the folds
accuracy_plot <- plot_cv_results(final_output$Accuracies)

# 4. Extract the top 10% of lines for a simulated breeding program
elite_lines <- extract_top_lines(final_output$Predicted_GEBVs, selection_intensity = 0.10)
print(head(elite_lines))

```

## 📂 Repository Structure

* **`1_data_simulation.R`**: Contains `simulate_genomic_data()` to generate synthetic marker matrices (`M`), haplotype matrices (`H`), and phenotypic vectors (`Y`) based on mixed genetic architectures.
* **`2_model_fitting.R`**: Contains the core logic for the models.
* `compute_gaussian_kernel()`: Generates the genetic distance matrix.
* `fit_standard_model()`: Fits standard GBLUP and RKHS models.
* `fit_gwas_assisted_model()`: Performs internal GWAS and fits the hybrid fixed/random models.


* **`3_orchestrator.R`**: Contains `run_modular_cv()`, which manages data splitting, phenotype masking, and logging model performance. Includes utility functions for plotting and extracting Genomic Estimated Breeding Values (GEBVs).

## 📊 Models Compared

The pipeline runs and compares four distinct models to demonstrate the stepwise improvement in accuracy:

1. **Standard GBLUP:** Traditional Bayesian Ridge Regression; strictly additive, all markers penalized equally.
2. **Standard RKHS:** Non-linear Gaussian kernel; captures non-additive variance, equal penalization.
3. **GWAS-Assisted GBLUP:** Top SNPs treated as fixed effects; remaining SNPs treated as random (additive).
4. **Haplotype-GWAS-Assisted RKHS (Best):** Top haplotype blocks treated as fixed effects; polygenic background handled by the non-linear Gaussian kernel.

## 💡 Using Your Own Data

To use empirical data instead of simulated data, bypass `simulate_genomic_data()`. Simply format your genotype data as a numeric matrix `M` (e.g., -1, 0, 1), compute your haplotype matrix `H`, and provide a phenotypic vector `Y`. Feed these directly into the `run_modular_cv()` loop.

## 📜 License

This project is licensed under the MIT License 

```

```
