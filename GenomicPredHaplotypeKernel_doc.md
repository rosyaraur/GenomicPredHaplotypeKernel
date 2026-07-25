# Documentation: GWAS-Assisted Genomic Prediction Pipeline

This documentation outlines the methodology and R implementation for a modular genomic prediction pipeline. The framework is designed to integrate Genome-Wide Association Study (GWAS) results into Genomic Selection (GS) models, specifically utilizing haplotype blocks and Reproducing Kernel Hilbert Space (RKHS) regressions, as described by Sehgal et al. (2020).

---

## 1. Methodological Overview

Standard genomic prediction models, such as Ridge Regression BLUP (rrBLUP/GBLUP), assume all genetic markers contribute equally to a trait's variance. They apply a uniform penalty, shrinking all marker effects toward zero to prevent overfitting. However, for traits governed by a few major quantitative trait loci (QTLs) alongside a highly polygenic background (e.g., grain yield), this uniform shrinkage artificially suppresses the true impact of the major loci.

This pipeline resolves this by adopting a **two-part mixed model architecture**:

1. **Fixed Effects (Major Loci):** The most significant markers or haplotype blocks identified via GWAS are modeled without shrinkage.
2. **Random Effects (Polygenic Background):** The remaining genetic variance is captured using penalized regression or non-linear kernel methods.

### 1.1 Haplotype Construction

Single Nucleotide Polymorphisms (SNPs) are often inherited together in blocks due to linkage disequilibrium. By grouping adjacent SNPs into haplotypes (e.g., converting independent alleles into a single regional state like `-1_0_1`), the model can capture localized epistatic interactions. The pipeline converts these blocks into a design matrix of dummy variables representing the presence or absence of specific regional alleles.

### 1.2 Gaussian Kernel (RKHS)

Standard GBLUP models strictly additive genetic variance. To capture non-linear, non-additive genetic variance (cryptic epistasis), the pipeline employs RKHS regression using a Gaussian kernel. The kernel defines the genetic similarity matrix $K$ based on the squared Euclidean distance $d_{ij}^2$ between individuals $i$ and $j$:

$$K_{ij} = \exp(-\theta d_{ij}^2)$$

Where $\theta$ is a bandwidth parameter, typically defined as the inverse of the median squared distance across the population.

### 1.3 Strict Cross-Validation (Preventing Data Leakage)

A critical methodological constraint is that GWAS must be performed **exclusively on the training fold** during cross-validation. Identifying top markers using the entire dataset prior to splitting allows the fixed effects to "see" the test data, severely inflating prediction accuracy estimates. This pipeline rigorously encapsulates the GWAS step within the training loop of each fold.

---

## 2. Dependencies

The pipeline requires the following R packages:

* **`BGLR`**: Bayesian Generalized Linear Regression (handles mixed fixed/random components and RKHS kernels natively).
* **`ggplot2`**: For visualizing cross-validation accuracy distributions.
* **`tidyr`**: For data reshaping prior to plotting.

---

## 3. Function Reference

### Module 1: Data Simulation

**`simulate_genomic_data(n_ind = 300, n_markers = 600, block_size = 3, n_qtl = 10, h2 = 0.5)`**
Simulates a synthetic population with a complex trait architecture.

* **`n_ind`**: Number of individuals.
* **`n_markers`**: Total number of simulated SNPs.
* **`block_size`**: Number of contiguous SNPs to group into a single haplotype block.
* **`n_qtl`**: Number of causal haplotype blocks driving the true breeding value.
* **`h2`**: Target narrow-sense heritability.
* **Returns**: A list containing `M` (SNP matrix), `H` (Haplotype dummy matrix), and `Y` (Phenotypic vector).

### Module 2: Utilities

**`compute_gaussian_kernel(M)`**
Calculates the non-linear genetic similarity matrix.

* **`M`**: The numeric matrix of SNP genotypes.
* **Returns**: A square symmetric kernel matrix $K$.

### Module 3: Standard Model Fitting

**`fit_standard_model(Y_masked, matrix_input, model_type = "BRR", nIter = 2000, burnIn = 500)`**
Fits baseline genomic prediction models treating all loci equally.

* **`Y_masked`**: Phenotypic vector with testing set indices replaced by `NA`.
* **`matrix_input`**: The genetic data (either `M` for GBLUP or `K` for RKHS).
* **`model_type`**: String specifying the BGLR model (`"BRR"` for Bayesian Ridge Regression, `"RKHS"` for Gaussian Kernel).
* **Returns**: A vector of predicted phenotypes (`yHat`).

### Module 4: GWAS-Assisted Model Fitting

**`fit_gwas_assisted_model(Y, Y_masked, train_idx, fixed_matrix, random_matrix = NULL, random_kernel = NULL, n_fixed = 5, nIter = 2000, burnIn = 500)`**
The core analytical function. Performs training-set-exclusive GWAS, identifies top loci, and solves the mixed architecture.

* **`Y`**: The complete phenotypic vector (used only at `train_idx` for GWAS).
* **`Y_masked`**: Phenotypic vector with testing set indices set to `NA` (for BGLR solving).
* **`train_idx`**: Vector of integer indices representing the training set.
* **`fixed_matrix`**: The matrix (SNP or Haplotype) to screen for top GWAS hits.
* **`random_matrix` / `random_kernel**`: The polygenic background matrix. One must be provided. If `random_matrix` is used, top loci are removed from it to prevent double-fitting.
* **`n_fixed`**: The number of top loci to extract and model as fixed effects.
* **Returns**: A vector of predicted phenotypes (`yHat`).

### Module 5: Orchestration

**`run_modular_cv(k_folds = 5, n_fixed = 5)`**
Manages the $k$-fold cross-validation loop, coordinates data masking, and logs performance across four distinct models: Standard GBLUP, Standard RKHS, GWAS-Assisted GBLUP, and Haplotype-GWAS RKHS.

* **`k_folds`**: Number of cross-validation folds.
* **`n_fixed`**: Number of major loci to fix in the assisted models.
* **Returns**: A list containing `Accuracies` (matrix of fold results), `True_Phenotypes`, and `Predicted_GEBVs` (unbiased out-of-sample predictions for the entire population based on the best model).

### Module 6: Selection Extraction

**`extract_top_lines(gebv_vector, selection_intensity = 0.10)`**
Ranks the population based on Genomic Estimated Breeding Values (GEBVs) and extracts the elite tier.

* **`gebv_vector`**: The out-of-sample predictions returned by the orchestrator.
* **`selection_intensity`**: Proportion of the top individuals to retain (e.g., `0.10` for top 10%).
* **Returns**: A sorted data frame mapping line IDs to their GEBVs.

### Module 7: Visualization

**`plot_cv_results(cv_matrix)`**
Generates publication-ready boxplots comparing the distribution of prediction accuracies across the tested methodologies.

* **`cv_matrix`**: The `Accuracies` matrix returned by `run_modular_cv`.
* **Returns**: A `ggplot2` object.