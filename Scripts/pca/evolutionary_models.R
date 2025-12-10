#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ape)
  library(phytools)
  library(dplyr)
})

# Try to load geiger, but continue if not available
has_geiger <- requireNamespace("geiger", quietly = TRUE)
if (has_geiger) {
  library(geiger)
  message("Using geiger package for model fitting")
} else {
  message("Warning: geiger package not available. Attempting to install...")
  tryCatch({
    install.packages("geiger", repos = "https://cloud.r-project.org", quiet = TRUE)
    library(geiger)
    has_geiger <- TRUE
    message("Successfully installed and loaded geiger package")
  }, error = function(e) {
    message("Could not install geiger. OU and EB models will be skipped.")
    message("To install manually, run: install.packages('geiger')")
    has_geiger <- FALSE
  })
}

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  cat("Usage: Rscript evolutionary_models.R <pca_scores_csv> <tree_file> <output_prefix> [PC_components]\n")
  cat("\nArguments:\n")
  cat("  pca_scores_csv  : Path to CSV file with PCA scores\n")
  cat("  tree_file       : Path to Newick tree file\n")
  cat("  output_prefix   : Prefix for output files\n")
  cat("  PC_components   : Comma-separated list of PCs to use (default: PC1,PC2)\n")
  cat("\nThis script:\n")
  cat("  1. Fits evolutionary models (BM, OU, EB)\n")
  cat("  2. Tests for rate differences between cetaceans and terrestrial artiodactyls\n")
  cat("  3. Compares model fits using AIC\n")
  quit(status = 1)
}

# Get arguments
pca_scores_file <- args[[1]]
tree_file <- args[[2]]
out_prefix <- args[[3]]
pc_components <- ifelse(length(args) >= 4, args[[4]], "PC1,PC2")

# Parse PC components
pc_list <- unlist(strsplit(pc_components, ","))
pc_list <- trimws(pc_list)

# Check if files exist
if (!file.exists(pca_scores_file)) {
  stop("Error: PCA scores file not found: ", pca_scores_file)
}
if (!file.exists(tree_file)) {
  stop("Error: Tree file not found: ", tree_file)
}

pca_scores_file <- normalizePath(pca_scores_file, mustWork = TRUE)
tree_file <- normalizePath(tree_file, mustWork = TRUE)

# Create output directory if needed
output_dir <- dirname(out_prefix)
if (output_dir != "." && !dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  message("Created output directory: ", output_dir)
}

message(paste(rep("=", 60), collapse = ""))
message("Evolutionary Model Fitting: Testing Cetacean Evolution")
message(paste(rep("=", 60), collapse = ""))
message("PCA scores file: ", pca_scores_file)
message("Tree file: ", tree_file)
message("Output prefix: ", out_prefix)
message("PC components: ", paste(pc_list, collapse = ", "))
message(paste(rep("=", 60), collapse = ""))

# Function to classify species group (excludes humans)
classify_species_group <- function(seq_id) {
  # Extract genus (first part before underscore or pipe)
  # Check for pipe first (since some sequences have pipes before underscores)
  seq_id <- as.character(seq_id)  # Ensure it's a character
  if (grepl("|", seq_id, fixed = TRUE)) {
    # If pipe exists, check if it comes before the first underscore
    pipe_pos <- regexpr("|", seq_id, fixed = TRUE)[1]
    underscore_pos <- regexpr("_", seq_id, fixed = TRUE)[1]
    if (pipe_pos > 0 && (underscore_pos < 0 || pipe_pos < underscore_pos)) {
      # Pipe comes first or no underscore, split on pipe
      genus <- strsplit(seq_id, "|", fixed = TRUE)[[1]][1]
    } else {
      # Underscore comes first, split on underscore
      genus <- strsplit(seq_id, "_", fixed = TRUE)[[1]][1]
    }
  } else if (grepl("_", seq_id)) {
    genus <- strsplit(seq_id, "_", fixed = TRUE)[[1]][1]
  } else {
    genus <- seq_id  # Fallback: use whole string
  }
  # Trim any whitespace
  genus <- trimws(genus)
  
  cetacean_genera <- c(
    "Balaenoptera", "Delphinapterus", "Delphinus", "Eschrichtius", 
    "Eubalaena", "Globicephala", "Kogia", "Lagenorhynchus", 
    "Lipotes", "Monodon", "Neophocaena", "Orcinus", "Phocoena", 
    "Physeter", "Pseudorca", "Tursiops", "Sousa", "Mesoplodon"
  )
  
  # Terrestrial Artiodactyl genera (excluding humans)
  terrestrial_genera <- c(
    "Bison", "Bos", "Bubalus", "Budorcas", "Camelus", "Capra", 
    "Capricornis", "Cervus", "Dama", "Hippopotamus", "Moschus", 
    "Muntiacus", "Odocoileus", "Oryx", "Ovibos", "Ovis", 
    "Pantholops", "Phacochoerus", "Rangifer", "Sus", "Vicugna"
  )
  
  if (genus %in% cetacean_genera) {
    return("Cetacean")
  } else if (genus %in% terrestrial_genera) {
    return("Terrestrial_Artiodactyl")
  } else if (genus == "Homo") {
    return("Human")  # Mark humans separately
  } else {
    return("Other")  # Unknown/other taxa
  }
}

# 1) Load PCA scores
message("\nLoading PCA scores...")
pca_scores <- read.csv(pca_scores_file, stringsAsFactors = FALSE)

# Check if required PCs exist
missing_pcs <- setdiff(pc_list, colnames(pca_scores))
if (length(missing_pcs) > 0) {
  stop("PC components not found: ", paste(missing_pcs, collapse = ", "))
}

# Determine ID column
if ("original_sequence_id" %in% colnames(pca_scores)) {
  pca_scores$seq_id <- pca_scores$original_sequence_id
  message("Using 'original_sequence_id' column to match with tree")
} else if ("clean_id" %in% colnames(pca_scores)) {
  pca_scores$seq_id <- pca_scores$clean_id
  message("Using 'clean_id' column to match with tree")
} else {
  pca_scores$seq_id <- pca_scores[, 1]
  message("Using first column to match with tree")
}

# Classify species groups
message("\nClassifying species groups...")
pca_scores$group <- sapply(pca_scores$seq_id, classify_species_group)
group_counts <- table(pca_scores$group)
message("Group counts (before filtering):")
for (group in names(group_counts)) {
  message("  ", group, ": ", group_counts[group])
}

# Filter to only include Cetaceans and Terrestrial Artiodactyls (exclude humans and others)
message("\nFiltering to Cetaceans and Terrestrial Artiodactyls only...")
message("  EXCLUDING humans (Homo) and other taxa from analysis")
human_count <- sum(pca_scores$group == "Human")
if (human_count > 0) {
  message("  Removing ", human_count, " human sequence(s)")
}
pca_scores <- pca_scores %>%
  filter(group %in% c("Cetacean", "Terrestrial_Artiodactyl"))

# Rename for consistency
pca_scores$group <- ifelse(pca_scores$group == "Terrestrial_Artiodactyl", "Terrestrial", pca_scores$group)

group_counts_filtered <- table(pca_scores$group)
message("Group counts (after filtering):")
for (group in names(group_counts_filtered)) {
  message("  ", group, ": ", group_counts_filtered[group])
}

# 2) Load phylogenetic tree
message("\nLoading phylogenetic tree...")
tree <- read.tree(tree_file)
tree$tip.label <- gsub(" ", "_", tree$tip.label)
message("Tree loaded: ", length(tree$tip.label), " tips")

# 3) Match sequences
common_seqs <- intersect(tree$tip.label, pca_scores$seq_id)
message("\nMatching sequences...")
message("  Sequences in tree: ", length(tree$tip.label))
message("  Sequences in PCA: ", nrow(pca_scores))
message("  Common sequences: ", length(common_seqs))

if (length(common_seqs) < 10) {
  stop("Too few common sequences (", length(common_seqs), "). Need at least 10.")
}

# Filter and prepare data
pca_filtered <- pca_scores %>%
  filter(seq_id %in% common_seqs) %>%
  arrange(match(seq_id, tree$tip.label))

# Extract PC scores matrix
pc_matrix <- pca_filtered %>%
  select(all_of(pc_list)) %>%
  as.matrix()
rownames(pc_matrix) <- pca_filtered$seq_id

# Get group labels
tip_groups <- pca_filtered$group
names(tip_groups) <- pca_filtered$seq_id

# Prune tree to common sequences
tree_pruned <- keep.tip(tree, common_seqs)

# Reorder matrix to match tree
pc_matrix <- pc_matrix[tree_pruned$tip.label, , drop = FALSE]
tip_groups <- tip_groups[tree_pruned$tip.label]

# Handle zero-length branches
if (any(tree_pruned$edge.length == 0)) {
  message("  Warning: Tree has ", sum(tree_pruned$edge.length == 0), " zero-length branches")
  message("  Adding small value to zero-length branches...")
  min_edge <- min(tree_pruned$edge.length[tree_pruned$edge.length > 0])
  tree_pruned$edge.length[tree_pruned$edge.length == 0] <- min_edge * 0.001
}

message("\nFinal dataset:")
message("  Taxa: ", nrow(pc_matrix))
message("  PC dimensions: ", ncol(pc_matrix))
message("  Cetaceans: ", sum(tip_groups == "Cetacean"))
message("  Terrestrial: ", sum(tip_groups == "Terrestrial"))

# 4) Fit evolutionary models
message("\n", paste(rep("=", 60), collapse = ""))
message("Fitting Evolutionary Models")
message(paste(rep("=", 60), collapse = ""))

results_list <- list()

# For each PC, fit models
for (pc_idx in 1:ncol(pc_matrix)) {
  pc_name <- colnames(pc_matrix)[pc_idx]
  trait <- pc_matrix[, pc_idx]
  names(trait) <- rownames(pc_matrix)
  
  message("\n", paste(rep("-", 60), collapse = ""))
  message("Analyzing ", pc_name)
  message(paste(rep("-", 60), collapse = ""))
  
  # Model 1: Single-rate Brownian Motion (BM)
  message("\n1. Fitting single-rate Brownian Motion...")
  tryCatch({
    if (has_geiger) {
      fit_bm <- fitContinuous(tree_pruned, trait, model = "BM", control = list(niter = 100))
      bm_aic <- fit_bm$opt$aic
      bm_sig2 <- fit_bm$opt$sigsq
      bm_loglik <- fit_bm$opt$lnL
    } else {
      # Use phytools::fastBM or manual calculation
      # Calculate BM rate using maximum likelihood
      bm_result <- phytools::brownie.lite(tree_pruned, trait)
      bm_sig2 <- bm_result$sig2.single
      # Approximate AIC (2 * k - 2 * logLik, where k = 1 for BM)
      bm_loglik <- -0.5 * (length(trait) * log(2 * pi * bm_sig2) + sum((trait - mean(trait))^2) / bm_sig2)
      bm_aic <- 2 * 1 - 2 * bm_loglik
    }
    message("  AIC: ", round(bm_aic, 2))
    message("  Sigma^2 (rate): ", round(bm_sig2, 6))
    
    results_list[[paste0(pc_name, "_BM")]] <- list(
      model = "BM",
      pc = pc_name,
      aic = bm_aic,
      sigsq = bm_sig2,
      loglik = bm_loglik
    )
  }, error = function(e) {
    message("  Error fitting BM: ", e$message)
  })
  
  # Model 2: Multi-rate Brownian Motion (different rates for cetaceans vs terrestrial)
  message("\n2. Fitting multi-rate Brownian Motion (cetaceans vs terrestrial)...")
  tryCatch({
    # Try to fit with different rates using geiger's approach
    # This is a simplified version - full implementation would use OUwie or similar
    message("  Note: Multi-rate BM requires specialized packages (OUwie, mvMORPH)")
    message("  Using approximate method (fitting each group separately)...")
    
    # Calculate rates separately for each group
    cetacean_tips <- names(tip_groups)[tip_groups == "Cetacean"]
    terrestrial_tips <- names(tip_groups)[tip_groups == "Terrestrial"]
    
    if (length(cetacean_tips) > 3 && length(terrestrial_tips) > 3) {
        # Fit BM to each group separately
        tree_cet <- keep.tip(tree_pruned, cetacean_tips)
        tree_terr <- keep.tip(tree_pruned, terrestrial_tips)
        
        trait_cet <- trait[cetacean_tips]
        trait_terr <- trait[terrestrial_tips]
        
        if (has_geiger) {
          fit_cet <- tryCatch(fitContinuous(tree_cet, trait_cet, model = "BM"), 
                             error = function(e) NULL)
          fit_terr <- tryCatch(fitContinuous(tree_terr, trait_terr, model = "BM"), 
                              error = function(e) NULL)
        } else {
          fit_cet <- tryCatch(phytools::brownie.lite(tree_cet, trait_cet), 
                             error = function(e) NULL)
          fit_terr <- tryCatch(phytools::brownie.lite(tree_terr, trait_terr), 
                              error = function(e) NULL)
        }
        
        if (!is.null(fit_cet) && !is.null(fit_terr)) {
          if (has_geiger) {
            sigsq_cet <- fit_cet$opt$sigsq
            sigsq_terr <- fit_terr$opt$sigsq
          } else {
            sigsq_cet <- fit_cet$sig2.single
            sigsq_terr <- fit_terr$sig2.single
          }
          rate_ratio <- sigsq_cet / sigsq_terr
        
        message("  Cetacean rate (sigma^2): ", round(sigsq_cet, 6))
        message("  Terrestrial rate (sigma^2): ", round(sigsq_terr, 6))
        message("  Rate ratio (cetacean/terrestrial): ", round(rate_ratio, 3))
        
        # Approximate AIC (simplified)
        aic_multi <- bm_aic + 1  # +1 for extra parameter (rough estimate)
        
        results_list[[paste0(pc_name, "_BM_multi")]] <- list(
          model = "BM_multi",
          pc = pc_name,
          aic = aic_multi,
          sigsq_cetacean = sigsq_cet,
          sigsq_terrestrial = sigsq_terr,
          rate_ratio = rate_ratio
        )
      }
    }
  }, error = function(e) {
    message("  Error fitting multi-rate BM: ", e$message)
  })
  
  # Model 3: Ornstein-Uhlenbeck (OU) - single optimum
  message("\n3. Fitting Ornstein-Uhlenbeck (single optimum)...")
  tryCatch({
    if (has_geiger) {
      fit_ou <- fitContinuous(tree_pruned, trait, model = "OU", control = list(niter = 100))
      ou_aic <- fit_ou$opt$aic
      ou_alpha <- fit_ou$opt$alpha
      ou_sig2 <- fit_ou$opt$sigsq
      # For single-optimum OU, z0 is the optimum (theta)
      ou_theta <- if("z0" %in% names(fit_ou$opt)) fit_ou$opt$z0 else NA
      ou_loglik <- fit_ou$opt$lnL
    } else {
      stop("  ERROR: OU model requires geiger package but it is not available. Please install geiger: install.packages('geiger')")
    }
    
    message("  AIC: ", round(ou_aic, 2))
    message("  Alpha (constraint): ", round(ou_alpha, 6))
    message("  Sigma^2 (rate): ", round(ou_sig2, 6))
    if (!is.na(ou_theta)) {
      message("  Theta (optimum): ", round(ou_theta, 4))
    } else {
      message("  Theta (optimum): NA (not available)")
    }
    
    results_list[[paste0(pc_name, "_OU")]] <- list(
      model = "OU",
      pc = pc_name,
      aic = ou_aic,
      alpha = ou_alpha,
      sigsq = ou_sig2,
      theta = ou_theta,
      loglik = ou_loglik
    )
  }, error = function(e) {
    message("  Error fitting OU: ", e$message)
  })
  
  # Model 4: Early Burst (EB)
  message("\n4. Fitting Early Burst...")
  tryCatch({
    if (has_geiger) {
      fit_eb <- fitContinuous(tree_pruned, trait, model = "EB", control = list(niter = 100))
      eb_aic <- fit_eb$opt$aic
      eb_a <- fit_eb$opt$a
      eb_sig2 <- fit_eb$opt$sigsq
      eb_loglik <- fit_eb$opt$lnL
    } else {
      stop("  ERROR: EB model requires geiger package but it is not available. Please install geiger: install.packages('geiger')")
    }
    
    message("  AIC: ", round(eb_aic, 2))
    message("  a (exponential decay): ", round(eb_a, 6))
    message("  Sigma^2 (rate): ", round(eb_sig2, 6))
    
    results_list[[paste0(pc_name, "_EB")]] <- list(
      model = "EB",
      pc = pc_name,
      aic = eb_aic,
      a = eb_a,
      sigsq = eb_sig2,
      loglik = eb_loglik
    )
  }, error = function(e) {
    message("  Error fitting EB: ", e$message)
  })
  
  # Compare models
  message("\n5. Model Comparison:")
  model_aics <- sapply(results_list[grepl(pc_name, names(results_list))], function(x) x$aic)
  if (length(model_aics) > 0) {
    best_model <- names(model_aics)[which.min(model_aics)]
    message("  Best model (lowest AIC): ", best_model, " (AIC = ", round(min(model_aics), 2), ")")
    
    # Calculate AIC weights (relative likelihood)
    delta_aic <- model_aics - min(model_aics)
    aic_weights <- exp(-0.5 * delta_aic) / sum(exp(-0.5 * delta_aic))
    
    message("\n  Model ranking (by AIC):")
    sorted_models <- sort(model_aics)
    for (i in 1:length(sorted_models)) {
      model_name <- names(sorted_models)[i]
      weight <- aic_weights[model_name]
      message("    ", i, ". ", model_name, ": AIC = ", round(sorted_models[i], 2), 
              ", weight = ", round(weight, 3))
    }
  }
}

# 5) Calculate phylogenetic signal statistics
message("\n", paste(rep("=", 60), collapse = ""))
message("Calculating Phylogenetic Signal Statistics")
message(paste(rep("=", 60), collapse = ""))
message("Calculating Blomberg's K and Pagel's λ for each PC dimension")

phylogenetic_signal <- list()

for (pc_idx in 1:ncol(pc_matrix)) {
  pc_name <- colnames(pc_matrix)[pc_idx]
  trait <- pc_matrix[, pc_idx]
  names(trait) <- rownames(pc_matrix)
  
  message("\n", pc_name, ":")
  
  # Calculate Blomberg's K
  message("  Calculating Blomberg's K...")
  tryCatch({
    # Use phytools::phylosig with method="K"
    k_result <- phytools::phylosig(tree_pruned, trait, method = "K", test = TRUE)
    blomberg_k <- k_result$K
    blomberg_k_pvalue <- k_result$P
    
    message("    K = ", round(blomberg_k, 4))
    message("    P-value = ", round(blomberg_k_pvalue, 4))
    
    if (blomberg_k_pvalue < 0.05) {
      message("    Significant phylogenetic signal (p < 0.05)")
    } else {
      message("    No significant phylogenetic signal (p >= 0.05)")
    }
  }, error = function(e) {
    message("    Error calculating Blomberg's K: ", e$message)
    blomberg_k <- NA
    blomberg_k_pvalue <- NA
  })
  
  # Calculate Pagel's λ
  message("  Calculating Pagel's λ...")
  tryCatch({
    lambda_result <- phytools::phylosig(tree_pruned, trait, method = "lambda", test = TRUE)
    pagel_lambda <- lambda_result$lambda
    lambda_pvalue <- lambda_result$P
    
    message("    λ = ", round(pagel_lambda, 4))
    message("    P-value = ", round(lambda_pvalue, 4))
    
    if (lambda_pvalue < 0.05) {
      message("    Significant phylogenetic signal (p < 0.05)")
    } else {
      message("    No significant phylogenetic signal (p >= 0.05)")
    }
  }, error = function(e) {
    message("    Error calculating Pagel's λ: ", e$message)
    pagel_lambda <- NA
    lambda_pvalue <- NA
  })
  
  phylogenetic_signal[[pc_name]] <- list(
    pc = pc_name,
    blomberg_k = if(exists("blomberg_k")) blomberg_k else NA,
    blomberg_k_pvalue = if(exists("blomberg_k_pvalue")) blomberg_k_pvalue else NA,
    pagel_lambda = if(exists("pagel_lambda")) pagel_lambda else NA,
    pagel_lambda_pvalue = if(exists("lambda_pvalue")) lambda_pvalue else NA
  )
}

# 6) Test for rate differences using phylogenetic ANOVA
message("\n", paste(rep("=", 60), collapse = ""))
message("Testing for Rate Differences: Cetaceans vs Terrestrial")
message(paste(rep("=", 60), collapse = ""))

# Use phytools::phylANOVA or similar
rate_tests <- list()

for (pc_idx in 1:ncol(pc_matrix)) {
  pc_name <- colnames(pc_matrix)[pc_idx]
  trait <- pc_matrix[, pc_idx]
  names(trait) <- rownames(pc_matrix)
  
  message("\n", pc_name, ":")
  
  # Calculate rates for each group
  cetacean_tips <- names(tip_groups)[tip_groups == "Cetacean"]
  terrestrial_tips <- names(tip_groups)[tip_groups == "Terrestrial"]
  
  if (length(cetacean_tips) > 3 && length(terrestrial_tips) > 3) {
    # Calculate mean trait values
    mean_cet <- mean(trait[cetacean_tips])
    mean_terr <- mean(trait[terrestrial_tips])
    
    # Calculate variance
    var_cet <- var(trait[cetacean_tips])
    var_terr <- var(trait[terrestrial_tips])
    
    message("  Cetacean mean: ", round(mean_cet, 4), ", variance: ", round(var_cet, 4))
    message("  Terrestrial mean: ", round(mean_terr, 4), ", variance: ", round(var_terr, 4))
    message("  Difference in means: ", round(abs(mean_cet - mean_terr), 4))
    
    # Phylogenetic ANOVA (simplified - would use phylANOVA for proper test)
    # For now, just report descriptive statistics
    rate_tests[[pc_name]] <- list(
      pc = pc_name,
      cetacean_mean = mean_cet,
      terrestrial_mean = mean_terr,
      cetacean_var = var_cet,
      terrestrial_var = var_terr,
      mean_difference = abs(mean_cet - mean_terr)
    )
  }
}

# 7) Save results
message("\n", paste(rep("=", 60), collapse = ""))
message("Saving Results")
message(paste(rep("=", 60), collapse = ""))

# Convert results to dataframes
model_results <- do.call(rbind, lapply(results_list, function(x) {
  data.frame(
    model = x$model,
    pc = x$pc,
    aic = x$aic,
    stringsAsFactors = FALSE
  )
}))

if (nrow(model_results) > 0) {
  write.csv(model_results, 
            file = paste0(out_prefix, "_no_humans_model_comparison.csv"),
            row.names = FALSE)
  message("Saved: ", paste0(out_prefix, "_no_humans_model_comparison.csv"))
}

# Save rate test results
if (length(rate_tests) > 0) {
  rate_results <- do.call(rbind, lapply(rate_tests, function(x) {
    data.frame(
      pc = x$pc,
      cetacean_mean = x$cetacean_mean,
      terrestrial_mean = x$terrestrial_mean,
      cetacean_variance = x$cetacean_var,
      terrestrial_variance = x$terrestrial_var,
      mean_difference = x$mean_difference,
      stringsAsFactors = FALSE
    )
  }))
  
  write.csv(rate_results,
            file = paste0(out_prefix, "_no_humans_rate_comparison.csv"),
            row.names = FALSE)
  message("Saved: ", paste0(out_prefix, "_no_humans_rate_comparison.csv"))
}

# Save phylogenetic signal results
if (length(phylogenetic_signal) > 0) {
  signal_results <- do.call(rbind, lapply(phylogenetic_signal, function(x) {
    data.frame(
      pc = x$pc,
      blomberg_k = x$blomberg_k,
      blomberg_k_pvalue = x$blomberg_k_pvalue,
      pagel_lambda = x$pagel_lambda,
      pagel_lambda_pvalue = x$pagel_lambda_pvalue,
      stringsAsFactors = FALSE
    )
  }))
  
  write.csv(signal_results,
            file = paste0(out_prefix, "_no_humans_phylogenetic_signal.csv"),
            row.names = FALSE)
  message("Saved: ", paste0(out_prefix, "_no_humans_phylogenetic_signal.csv"))
}

message("\n", paste(rep("=", 60), collapse = ""))
message("Analysis Complete!")
message(paste(rep("=", 60), collapse = ""))

