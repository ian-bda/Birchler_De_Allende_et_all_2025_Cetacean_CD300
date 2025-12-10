#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ape)
  library(phytools)
  library(dplyr)
  library(ggplot2)
  library(sp)
})

# Try to load geiger
has_geiger <- requireNamespace("geiger", quietly = TRUE)
if (has_geiger) {
  library(geiger)
  message("Using geiger package for EB model fitting and simulation")
} else {
  stop("ERROR: geiger package is required for EB simulation analysis. Please install: install.packages('geiger')")
}

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 4) {
  cat("Usage: Rscript eb_simulation_analysis.R <pca_scores_csv> <tree_file> <output_prefix> <n_simulations> [PC_x] [PC_y]\n")
  cat("\nArguments:\n")
  cat("  pca_scores_csv  : Path to CSV file with PCA scores\n")
  cat("  tree_file       : Path to ultrametric Newick tree file\n")
  cat("  output_prefix   : Prefix for output files\n")
  cat("  n_simulations   : Number of simulations (default: 1000)\n")
  cat("  PC_x            : Principal component for X-axis (default: PC1)\n")
  cat("  PC_y            : Principal component for Y-axis (default: PC2)\n")
  cat("\nThis script:\n")
  cat("  1. Fits EB models to PC1 and PC2\n")
  cat("  2. Simulates PC1 and PC2 under EB model 1000x\n")
  cat("  3. Calculates convex hull metrics for each simulation\n")
  cat("  4. Compares observed metrics to simulated distribution\n")
  cat("  5. Calculates p-values to test if observed is outside EB model predictions\n")
  quit(status = 1)
}

# Get arguments
pca_scores_file <- args[[1]]
tree_file <- args[[2]]
out_prefix <- args[[3]]
n_simulations <- ifelse(length(args) >= 4, as.integer(args[[4]]), 1000)
pc_x <- ifelse(length(args) >= 5, args[[5]], "PC1")
pc_y <- ifelse(length(args) >= 6, args[[6]], "PC2")

# Check if files exist
if (!file.exists(pca_scores_file)) {
  stop("Error: PCA scores file not found: ", pca_scores_file)
}
if (!file.exists(tree_file)) {
  stop("Error: Tree file not found: ", tree_file)
}

message(paste(rep("=", 60), collapse = ""))
message("EB Model Simulation Analysis")
message(paste(rep("=", 60), collapse = ""))
message("PCA scores file: ", pca_scores_file)
message("Tree file: ", tree_file)
message("Output prefix: ", out_prefix)
message("Number of simulations: ", n_simulations)
message("X-axis (PC): ", pc_x)
message("Y-axis (PC): ", pc_y)
message(paste(rep("=", 60), collapse = ""))

# Function to classify species group
classify_species_group <- function(seq_id) {
  if (is.na(seq_id) || seq_id == "") {
    return("Other")
  }
  
  seq_id_str <- as.character(seq_id)
  
  # Extract genus (first part before underscore or pipe)
  if (grepl("_", seq_id_str)) {
    genus <- strsplit(seq_id_str, "_")[[1]][1]
  } else if (grepl("\\|", seq_id_str)) {
    genus <- strsplit(seq_id_str, "\\|")[[1]][1]
  } else {
    genus <- seq_id_str
  }
  
  # Cetacean genera
  cetacean_genera <- c(
    "Balaenoptera", "Delphinapterus", "Delphinus", "Eschrichtius", 
    "Eubalaena", "Globicephala", "Kogia", "Lagenorhynchus", 
    "Lipotes", "Monodon", "Neophocaena", "Orcinus", "Phocoena", 
    "Physeter", "Pseudorca", "Tursiops", "Sousa", "Mesoplodon"
  )
  
  # Terrestrial Artiodactyl genera
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
    return("Human")
  } else {
    return("Other")
  }
}

# Function to calculate convex hull area (2D)
calculate_hull_area <- function(points) {
  if (nrow(points) < 3) {
    return(0)
  }
  
  points <- unique(points)
  if (nrow(points) < 3) {
    return(0)
  }
  
  hull <- chull(points[, 1], points[, 2])
  hull_points <- points[hull, ]
  
  n <- nrow(hull_points)
  if (n < 3) {
    return(0)
  }
  
  area <- 0
  for (i in 1:n) {
    j <- ifelse(i == n, 1, i + 1)
    area <- area + (hull_points[i, 1] * hull_points[j, 2] - hull_points[j, 1] * hull_points[i, 2])
  }
  area <- abs(area) / 2
  
  return(area)
}

# Function to calculate overlap between two convex hulls
calculate_hull_overlap <- function(hull1_points, hull2_points) {
  if (nrow(hull1_points) < 3 || nrow(hull2_points) < 3) {
    return(list(overlap_area = 0, overlap_percentage_hull1 = 0, overlap_percentage_hull2 = 0))
  }
  
  hull1 <- chull(hull1_points[, 1], hull1_points[, 2])
  hull2 <- chull(hull2_points[, 1], hull2_points[, 2])
  
  hull1_poly <- hull1_points[hull1, ]
  hull2_poly <- hull2_points[hull2, ]
  
  all_points <- rbind(hull1_points, hull2_points)
  combined_hull <- chull(all_points[, 1], all_points[, 2])
  combined_area <- calculate_hull_area(all_points[combined_hull, ])
  
  area1 <- calculate_hull_area(hull1_points)
  area2 <- calculate_hull_area(hull2_points)
  
  overlap_area <- max(0, area1 + area2 - combined_area)
  
  return(list(
    overlap_area = overlap_area,
    overlap_percentage_hull1 = ifelse(area1 > 0, (overlap_area / area1) * 100, 0),
    overlap_percentage_hull2 = ifelse(area2 > 0, (overlap_area / area2) * 100, 0)
  ))
}

# 1) Load PCA scores
message("\n1. Loading PCA scores...")
pca_scores <- read.csv(pca_scores_file, stringsAsFactors = FALSE)

# Determine ID column
if ("original_sequence_id" %in% colnames(pca_scores)) {
  pca_scores$seq_id <- pca_scores$original_sequence_id
  message("Using 'original_sequence_id' column")
} else if ("clean_id" %in% colnames(pca_scores)) {
  pca_scores$seq_id <- pca_scores$clean_id
  message("Using 'clean_id' column")
} else {
  pca_scores$seq_id <- pca_scores[, 1]
  message("Using first column as sequence ID")
}

# Check if required PCs exist
if (!pc_x %in% colnames(pca_scores)) {
  stop("PC '", pc_x, "' not found in PCA scores")
}
if (!pc_y %in% colnames(pca_scores)) {
  stop("PC '", pc_y, "' not found in PCA scores")
}

# Classify species groups
message("\n2. Classifying species groups...")
pca_scores$group <- sapply(pca_scores$seq_id, classify_species_group)

# Filter to cetaceans and terrestrial artiodactyls only (exclude humans)
message("\n3. Filtering data (excluding humans)...")
pca_scores <- pca_scores %>%
  filter(group %in% c("Cetacean", "Terrestrial_Artiodactyl"))

# Rename for consistency
pca_scores$group <- ifelse(pca_scores$group == "Terrestrial_Artiodactyl", 
                           "Non_Cetacean", 
                           pca_scores$group)

# 2) Load phylogenetic tree
message("\n4. Loading phylogenetic tree...")
tree <- read.tree(tree_file)
tree$tip.label <- gsub(" ", "_", tree$tip.label)
message("Tree loaded: ", length(tree$tip.label), " tips")

# Check if tree is ultrametric
if (!is.ultrametric(tree, tol = 1e-6)) {
  warning("WARNING: Tree may not be ultrametric. Simulation results may be inaccurate.")
}

# 3) Match sequences
common_seqs <- intersect(tree$tip.label, pca_scores$seq_id)
message("\n5. Matching sequences...")
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

# Extract PC scores
pc1_observed <- pca_filtered[[pc_x]]
pc2_observed <- pca_filtered[[pc_y]]
names(pc1_observed) <- pca_filtered$seq_id
names(pc2_observed) <- pca_filtered$seq_id

# Get group labels
tip_groups <- pca_filtered$group
names(tip_groups) <- pca_filtered$seq_id

# Prune tree to common sequences
tree_pruned <- keep.tip(tree, common_seqs)

# Reorder to match tree
pc1_observed <- pc1_observed[tree_pruned$tip.label]
pc2_observed <- pc2_observed[tree_pruned$tip.label]
tip_groups <- tip_groups[tree_pruned$tip.label]

# Handle zero-length branches
if (any(tree_pruned$edge.length == 0)) {
  message("  Warning: Tree has ", sum(tree_pruned$edge.length == 0), " zero-length branches")
  message("  Adding small value to zero-length branches...")
  min_edge <- min(tree_pruned$edge.length[tree_pruned$edge.length > 0])
  tree_pruned$edge.length[tree_pruned$edge.length == 0] <- min_edge * 0.001
}

message("\n6. Final dataset:")
message("  Taxa: ", length(tree_pruned$tip.label))
message("  Cetaceans: ", sum(tip_groups == "Cetacean"))
message("  Non-cetaceans: ", sum(tip_groups == "Non_Cetacean"))

# 4) Fit EB models to PC1 and PC2
message("\n", paste(rep("=", 60), collapse = ""))
message("Fitting EB Models to Observed Data")
message(paste(rep("=", 60), collapse = ""))

message("\nFitting EB model to ", pc_x, "...")
fit_eb_pc1 <- fitContinuous(tree_pruned, pc1_observed, model = "EB", control = list(niter = 100))
eb_pc1_a <- fit_eb_pc1$opt$a
eb_pc1_sigsq <- fit_eb_pc1$opt$sigsq
eb_pc1_aic <- fit_eb_pc1$opt$aic
message("  a (exponential decay): ", round(eb_pc1_a, 6))
message("  Sigma^2 (rate): ", round(eb_pc1_sigsq, 6))
message("  AIC: ", round(eb_pc1_aic, 2))

message("\nFitting EB model to ", pc_y, "...")
fit_eb_pc2 <- fitContinuous(tree_pruned, pc2_observed, model = "EB", control = list(niter = 100))
eb_pc2_a <- fit_eb_pc2$opt$a
eb_pc2_sigsq <- fit_eb_pc2$opt$sigsq
eb_pc2_aic <- fit_eb_pc2$opt$aic
message("  a (exponential decay): ", round(eb_pc2_a, 6))
message("  Sigma^2 (rate): ", round(eb_pc2_sigsq, 6))
message("  AIC: ", round(eb_pc2_aic, 2))

# 5) Calculate observed convex hull metrics
message("\n", paste(rep("=", 60), collapse = ""))
message("Calculating Observed Convex Hull Metrics")
message(paste(rep("=", 60), collapse = ""))

cetacean_tips <- names(tip_groups)[tip_groups == "Cetacean"]
non_cetacean_tips <- names(tip_groups)[tip_groups == "Non_Cetacean"]

observed_cetacean_points <- data.frame(
  x = pc1_observed[cetacean_tips],
  y = pc2_observed[cetacean_tips]
)

observed_non_cetacean_points <- data.frame(
  x = pc1_observed[non_cetacean_tips],
  y = pc2_observed[non_cetacean_tips]
)

observed_cetacean_area <- calculate_hull_area(observed_cetacean_points)
observed_non_cetacean_area <- calculate_hull_area(observed_non_cetacean_points)
observed_overlap <- calculate_hull_overlap(observed_cetacean_points, observed_non_cetacean_points)
observed_area_ratio <- ifelse(observed_non_cetacean_area > 0, 
                              observed_cetacean_area / observed_non_cetacean_area, 
                              NA)

message("\nObserved metrics:")
message("  Cetacean hull area: ", round(observed_cetacean_area, 6))
message("  Non-cetacean hull area: ", round(observed_non_cetacean_area, 6))
message("  Overlap area: ", round(observed_overlap$overlap_area, 6))
message("  Area ratio (cetacean/non-cetacean): ", round(observed_area_ratio, 4))

# 6) Simulate under EB model
message("\n", paste(rep("=", 60), collapse = ""))
message("Simulating Under EB Model (", n_simulations, " simulations)")
message(paste(rep("=", 60), collapse = ""))

# Initialize storage for simulation results
sim_results <- data.frame(
  simulation = 1:n_simulations,
  cetacean_area = numeric(n_simulations),
  non_cetacean_area = numeric(n_simulations),
  overlap_area = numeric(n_simulations),
  area_ratio = numeric(n_simulations),
  stringsAsFactors = FALSE
)

message("\nRunning simulations...")
pb <- txtProgressBar(min = 0, max = n_simulations, style = 3)

for (sim in 1:n_simulations) {
  setTxtProgressBar(pb, sim)
  
  # Simulate PC1 under EB model
  # EB model: rate(t) = sigsq * exp(a * t)
  # We transform the tree to account for rate decay, then simulate with BM
  tree_eb_pc1 <- tree_pruned
  node_heights_pc1 <- nodeHeights(tree_eb_pc1)
  
  # Transform edge lengths to account for EB rate decay
  # The variance accumulated along an edge from t_start to t_end is:
  # V = sigsq * integral from t_start to t_end of exp(a*t) dt
  #   = sigsq * (exp(a*t_end) - exp(a*t_start)) / a
  # For BM simulation: variance = sigsq * time
  # So we set: time_equivalent = V / sigsq = (exp(a*t_end) - exp(a*t_start)) / a
  for (i in 1:nrow(tree_eb_pc1$edge)) {
    t_start <- node_heights_pc1[i, 1]
    t_end <- node_heights_pc1[i, 2]
    
    if (abs(eb_pc1_a) > 1e-10) {
      # Calculate equivalent time for BM simulation
      # This accounts for the exponential rate decay in EB
      variance_accumulated <- (exp(eb_pc1_a * t_end) - exp(eb_pc1_a * t_start)) / eb_pc1_a
      # For BM: variance = sigsq * time, so equivalent time = variance / sigsq
      # But we want to preserve the total variance, so we use the variance directly
      # and simulate with sigsq = 1
      tree_eb_pc1$edge.length[i] <- variance_accumulated
    } else {
      # If a is very close to 0, EB reduces to BM (rate is constant)
      tree_eb_pc1$edge.length[i] <- tree_eb_pc1$edge.length[i]
    }
  }
  
  # Simulate with BM on transformed tree using sigsq = 1
  # (variance is already encoded in the transformed edge lengths)
  sim_pc1 <- fastBM(tree_eb_pc1, sig2 = eb_pc1_sigsq, internal = FALSE)
  names(sim_pc1) <- tree_pruned$tip.label
  
  # Simulate PC2 under EB model (same approach)
  tree_eb_pc2 <- tree_pruned
  node_heights_pc2 <- nodeHeights(tree_eb_pc2)
  
  for (i in 1:nrow(tree_eb_pc2$edge)) {
    t_start <- node_heights_pc2[i, 1]
    t_end <- node_heights_pc2[i, 2]
    
    if (abs(eb_pc2_a) > 1e-10) {
      variance_accumulated <- (exp(eb_pc2_a * t_end) - exp(eb_pc2_a * t_start)) / eb_pc2_a
      tree_eb_pc2$edge.length[i] <- variance_accumulated
    } else {
      tree_eb_pc2$edge.length[i] <- tree_eb_pc2$edge.length[i]
    }
  }
  
  sim_pc2 <- fastBM(tree_eb_pc2, sig2 = eb_pc2_sigsq, internal = FALSE)
  names(sim_pc2) <- tree_pruned$tip.label
  
  # Calculate convex hull metrics for this simulation
  sim_cetacean_points <- data.frame(
    x = sim_pc1[cetacean_tips],
    y = sim_pc2[cetacean_tips]
  )
  
  sim_non_cetacean_points <- data.frame(
    x = sim_pc1[non_cetacean_tips],
    y = sim_pc2[non_cetacean_tips]
  )
  
  sim_cetacean_area <- calculate_hull_area(sim_cetacean_points)
  sim_non_cetacean_area <- calculate_hull_area(sim_non_cetacean_points)
  sim_overlap <- calculate_hull_overlap(sim_cetacean_points, sim_non_cetacean_points)
  sim_area_ratio <- ifelse(sim_non_cetacean_area > 0, 
                           sim_cetacean_area / sim_non_cetacean_area, 
                           NA)
  
  sim_results$cetacean_area[sim] <- sim_cetacean_area
  sim_results$non_cetacean_area[sim] <- sim_non_cetacean_area
  sim_results$overlap_area[sim] <- sim_overlap$overlap_area
  sim_results$area_ratio[sim] <- sim_area_ratio
}

close(pb)

# 7) Calculate p-values
message("\n", paste(rep("=", 60), collapse = ""))
message("Statistical Testing: Comparing Observed to EB Model Predictions")
message(paste(rep("=", 60), collapse = ""))

# P-value = proportion of simulations that are more extreme than observed
# For one-tailed tests (observed > expected)
p_cetacean_area_high <- mean(sim_results$cetacean_area >= observed_cetacean_area)
p_cetacean_area_low <- mean(sim_results$cetacean_area <= observed_cetacean_area)
p_cetacean_area_two_tailed <- 2 * min(p_cetacean_area_high, p_cetacean_area_low)

p_non_cetacean_area_high <- mean(sim_results$non_cetacean_area >= observed_non_cetacean_area)
p_non_cetacean_area_low <- mean(sim_results$non_cetacean_area <= observed_non_cetacean_area)
p_non_cetacean_area_two_tailed <- 2 * min(p_non_cetacean_area_high, p_non_cetacean_area_low)

p_overlap_high <- mean(sim_results$overlap_area >= observed_overlap$overlap_area)
p_overlap_low <- mean(sim_results$overlap_area <= observed_overlap$overlap_area)
p_overlap_two_tailed <- 2 * min(p_overlap_high, p_overlap_low)

p_area_ratio_high <- mean(sim_results$area_ratio >= observed_area_ratio, na.rm = TRUE)
p_area_ratio_low <- mean(sim_results$area_ratio <= observed_area_ratio, na.rm = TRUE)
p_area_ratio_two_tailed <- 2 * min(p_area_ratio_high, p_area_ratio_low)

message("\nP-values (two-tailed):")
message("  Cetacean hull area: p = ", round(p_cetacean_area_two_tailed, 4))
message("  Non-cetacean hull area: p = ", round(p_non_cetacean_area_two_tailed, 4))
message("  Overlap area: p = ", round(p_overlap_two_tailed, 4))
message("  Area ratio: p = ", round(p_area_ratio_two_tailed, 4))

# 8) Create summary statistics
message("\n", paste(rep("=", 60), collapse = ""))
message("Summary Statistics")
message(paste(rep("=", 60), collapse = ""))

summary_stats <- data.frame(
  metric = c("Cetacean_area", "Non_cetacean_area", "Overlap_area", "Area_ratio"),
  observed = c(observed_cetacean_area, observed_non_cetacean_area, 
               observed_overlap$overlap_area, observed_area_ratio),
  mean_simulated = c(mean(sim_results$cetacean_area), 
                     mean(sim_results$non_cetacean_area),
                     mean(sim_results$overlap_area),
                     mean(sim_results$area_ratio, na.rm = TRUE)),
  sd_simulated = c(sd(sim_results$cetacean_area),
                   sd(sim_results$non_cetacean_area),
                   sd(sim_results$overlap_area),
                   sd(sim_results$area_ratio, na.rm = TRUE)),
  p_value_two_tailed = c(p_cetacean_area_two_tailed,
                        p_non_cetacean_area_two_tailed,
                        p_overlap_two_tailed,
                        p_area_ratio_two_tailed),
  stringsAsFactors = FALSE
)

message("\nObserved vs Simulated (mean ± SD):")
for (i in 1:nrow(summary_stats)) {
  message("\n  ", summary_stats$metric[i], ":")
  message("    Observed: ", round(summary_stats$observed[i], 6))
  message("    Simulated: ", round(summary_stats$mean_simulated[i], 6), 
          " ± ", round(summary_stats$sd_simulated[i], 6))
  message("    P-value: ", round(summary_stats$p_value_two_tailed[i], 4))
  
  # Calculate z-score
  z_score <- (summary_stats$observed[i] - summary_stats$mean_simulated[i]) / 
             summary_stats$sd_simulated[i]
  message("    Z-score: ", round(z_score, 4))
}

# 9) Create visualizations
message("\n", paste(rep("=", 60), collapse = ""))
message("Creating Visualizations")
message(paste(rep("=", 60), collapse = ""))

# Plot 1: Distribution of simulated metrics with observed overlaid
pdf(paste0(out_prefix, "_eb_simulation_distributions.pdf"), width = 12, height = 10)

par(mfrow = c(2, 2))

# Cetacean area
hist(sim_results$cetacean_area, breaks = 50, 
     main = "Cetacean Hull Area\n(EB Model Simulations)",
     xlab = "Hull Area", col = "lightblue", border = "black")
abline(v = observed_cetacean_area, col = "red", lwd = 2, lty = 2)
legend("topright", legend = c("Observed", "Simulated"), 
       col = c("red", "lightblue"), lty = c(2, 1), lwd = c(2, 5))

# Non-cetacean area
hist(sim_results$non_cetacean_area, breaks = 50,
     main = "Non-Cetacean Hull Area\n(EB Model Simulations)",
     xlab = "Hull Area", col = "lightcoral", border = "black")
abline(v = observed_non_cetacean_area, col = "red", lwd = 2, lty = 2)
legend("topright", legend = c("Observed", "Simulated"),
       col = c("red", "lightcoral"), lty = c(2, 1), lwd = c(2, 5))

# Overlap area
hist(sim_results$overlap_area, breaks = 50,
     main = "Overlap Area\n(EB Model Simulations)",
     xlab = "Overlap Area", col = "lightgreen", border = "black")
abline(v = observed_overlap$overlap_area, col = "red", lwd = 2, lty = 2)
legend("topright", legend = c("Observed", "Simulated"),
       col = c("red", "lightgreen"), lty = c(2, 1), lwd = c(2, 5))

# Area ratio
hist(sim_results$area_ratio, breaks = 50,
     main = "Area Ratio (Cetacean/Non-Cetacean)\n(EB Model Simulations)",
     xlab = "Area Ratio", col = "lightyellow", border = "black")
abline(v = observed_area_ratio, col = "red", lwd = 2, lty = 2)
legend("topright", legend = c("Observed", "Simulated"),
       col = c("red", "lightyellow"), lty = c(2, 1), lwd = c(2, 5))

dev.off()
message("Saved: ", paste0(out_prefix, "_eb_simulation_distributions.pdf"))

# Plot 2: Observed vs simulated scatter
pdf(paste0(out_prefix, "_eb_simulation_scatter.pdf"), width = 10, height = 8)

# Sample a subset of simulations for visualization
n_plot <- min(100, n_simulations)
plot_indices <- sample(1:n_simulations, n_plot)

plot(sim_results$cetacean_area[plot_indices], 
     sim_results$non_cetacean_area[plot_indices],
     xlab = "Cetacean Hull Area", ylab = "Non-Cetacean Hull Area",
     main = "Simulated Hull Areas Under EB Model",
     pch = 16, col = rgb(0, 0, 1, 0.3), cex = 0.8)
points(observed_cetacean_area, observed_non_cetacean_area,
       pch = 17, col = "red", cex = 2, lwd = 2)
legend("topright", legend = c("Observed", "Simulated (sample)"),
       pch = c(17, 16), col = c("red", "blue"), pt.cex = c(2, 0.8))

dev.off()
message("Saved: ", paste0(out_prefix, "_eb_simulation_scatter.pdf"))

# 10) Save results
message("\n", paste(rep("=", 60), collapse = ""))
message("Saving Results")
message(paste(rep("=", 60), collapse = ""))

# Save summary statistics
write.csv(summary_stats, 
          file = paste0(out_prefix, "_eb_simulation_summary.csv"),
          row.names = FALSE)
message("Saved: ", paste0(out_prefix, "_eb_simulation_summary.csv"))

# Save all simulation results
write.csv(sim_results,
          file = paste0(out_prefix, "_eb_simulation_results.csv"),
          row.names = FALSE)
message("Saved: ", paste0(out_prefix, "_eb_simulation_results.csv"))

# Save EB model parameters
eb_params <- data.frame(
  pc = c(pc_x, pc_y),
  a = c(eb_pc1_a, eb_pc2_a),
  sigsq = c(eb_pc1_sigsq, eb_pc2_sigsq),
  aic = c(eb_pc1_aic, eb_pc2_aic),
  stringsAsFactors = FALSE
)
write.csv(eb_params,
          file = paste0(out_prefix, "_eb_model_parameters.csv"),
          row.names = FALSE)
message("Saved: ", paste0(out_prefix, "_eb_model_parameters.csv"))

message("\n", paste(rep("=", 60), collapse = ""))
message("Analysis Complete!")
message(paste(rep("=", 60), collapse = ""))
message("\nInterpretation:")
message("  - P-value < 0.05: Observed metric is significantly different from EB model predictions")
message("  - P-value >= 0.05: Observed metric is consistent with EB model predictions")
message("  - This tests whether the observed morphological divergence is consistent with")
message("    an Early Burst model of evolution, or if additional processes are needed.")

