#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ape)
  library(phytools)
  library(dplyr)
  library(ggplot2)
  library(sp)
})

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  cat("Usage: Rscript convex_hull_analysis.R <pca_scores_csv> <output_prefix> [PC_x] [PC_y] [exclude_humans]\n")
  cat("\nArguments:\n")
  cat("  pca_scores_csv  : Path to CSV file with PCA scores\n")
  cat("  output_prefix   : Prefix for output files\n")
  cat("  PC_x            : Principal component for X-axis (default: PC1)\n")
  cat("  PC_y            : Principal component for Y-axis (default: PC2)\n")
  cat("  exclude_humans  : Exclude humans from analysis (default: TRUE)\n")
  cat("\nThis script calculates convex hulls for:\n")
  cat("  - Cetaceans (blue)\n")
  cat("  - Non-cetaceans/Terrestrial Artiodactyls (orange)\n")
  quit(status = 1)
}

pca_scores_file <- args[[1]]
out_prefix <- args[[2]]
pc_x <- ifelse(length(args) >= 3, args[[3]], "PC1")
pc_y <- ifelse(length(args) >= 4, args[[4]], "PC2")
exclude_humans <- ifelse(length(args) >= 5, as.logical(args[[5]]), TRUE)

# Check if files exist
if (!file.exists(pca_scores_file)) {
  stop("Error: PCA scores file not found: ", pca_scores_file)
}

message(paste(rep("=", 60), collapse = ""))
message("Convex Hull Analysis: Cetaceans vs Non-Cetaceans")
message(paste(rep("=", 60), collapse = ""))
message("PCA scores file: ", pca_scores_file)
message("Output prefix: ", out_prefix)
message("X-axis (PC): ", pc_x)
message("Y-axis (PC): ", pc_y)
message("Exclude humans: ", exclude_humans)
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
  
  # Remove duplicates
  points <- unique(points)
  if (nrow(points) < 3) {
    return(0)
  }
  
  # Calculate convex hull
  hull <- chull(points[, 1], points[, 2])
  hull_points <- points[hull, ]
  
  # Calculate area using shoelace formula
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

# Function to calculate convex hull perimeter (2D)
calculate_hull_perimeter <- function(points) {
  if (nrow(points) < 3) {
    return(0)
  }
  
  points <- unique(points)
  if (nrow(points) < 3) {
    return(0)
  }
  
  hull <- chull(points[, 1], points[, 2])
  hull_points <- points[hull, ]
  
  # Calculate perimeter
  n <- nrow(hull_points)
  perimeter <- 0
  for (i in 1:n) {
    j <- ifelse(i == n, 1, i + 1)
    dx <- hull_points[j, 1] - hull_points[i, 1]
    dy <- hull_points[j, 2] - hull_points[i, 2]
    perimeter <- perimeter + sqrt(dx^2 + dy^2)
  }
  
  return(perimeter)
}

# Function to calculate overlap between two convex hulls
calculate_hull_overlap <- function(hull1_points, hull2_points) {
  if (nrow(hull1_points) < 3 || nrow(hull2_points) < 3) {
    return(0)
  }
  
  # Get convex hulls
  hull1 <- chull(hull1_points[, 1], hull1_points[, 2])
  hull2 <- chull(hull2_points[, 1], hull2_points[, 2])
  
  hull1_poly <- hull1_points[hull1, ]
  hull2_poly <- hull2_points[hull2, ]
  
  # Check if any points from hull1 are inside hull2
  points_inside <- 0
  for (i in 1:nrow(hull1_points)) {
    if (point.in.polygon(hull1_points[i, 1], hull1_points[i, 2],
                         hull2_poly[, 1], hull2_poly[, 2]) > 0) {
      points_inside <- points_inside + 1
    }
  }
  
  # Check if any points from hull2 are inside hull1
  for (i in 1:nrow(hull2_points)) {
    if (point.in.polygon(hull2_points[i, 1], hull2_points[i, 2],
                         hull1_poly[, 1], hull1_poly[, 2]) > 0) {
      points_inside <- points_inside + 1
    }
  }
  
  # Calculate overlap area (simplified - using intersection of hulls)
  # This is an approximation
  all_points <- rbind(hull1_points, hull2_points)
  combined_hull <- chull(all_points[, 1], all_points[, 2])
  combined_area <- calculate_hull_area(all_points[combined_hull, ])
  
  area1 <- calculate_hull_area(hull1_points)
  area2 <- calculate_hull_area(hull2_points)
  
  # Overlap is the difference between combined area and sum of individual areas
  overlap_area <- max(0, area1 + area2 - combined_area)
  
  return(list(
    overlap_area = overlap_area,
    overlap_percentage_hull1 = ifelse(area1 > 0, (overlap_area / area1) * 100, 0),
    overlap_percentage_hull2 = ifelse(area2 > 0, (overlap_area / area2) * 100, 0),
    points_inside = points_inside
  ))
}

# Load PCA scores
message("\nLoading PCA scores...")
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
message("\nClassifying species groups...")
pca_scores$group <- sapply(pca_scores$seq_id, classify_species_group)
group_counts <- table(pca_scores$group)
message("Group counts:")
for (group in names(group_counts)) {
  message("  ", group, ": ", group_counts[group])
}

# Filter data
if (exclude_humans) {
  message("\nExcluding humans from analysis...")
  pca_scores <- pca_scores %>%
    filter(group != "Human")
}

# Filter to cetaceans and terrestrial artiodactyls only
pca_scores <- pca_scores %>%
  filter(group %in% c("Cetacean", "Terrestrial_Artiodactyl", "Other"))

# Rename for clarity
pca_scores$group <- ifelse(pca_scores$group == "Terrestrial_Artiodactyl", 
                           "Non_Cetacean", 
                           pca_scores$group)

# Extract PC coordinates
pca_scores$x <- pca_scores[[pc_x]]
pca_scores$y <- pca_scores[[pc_y]]

# Remove any missing values
pca_scores <- pca_scores %>%
  filter(!is.na(x) & !is.na(y))

message("\nFinal dataset:")
message("  Total sequences: ", nrow(pca_scores))
message("  Cetaceans: ", sum(pca_scores$group == "Cetacean"))
message("  Non-cetaceans: ", sum(pca_scores$group == "Non_Cetacean"))

# Separate groups
cetacean_data <- pca_scores %>%
  filter(group == "Cetacean") %>%
  select(x, y)

non_cetacean_data <- pca_scores %>%
  filter(group == "Non_Cetacean") %>%
  select(x, y)

# Check if we have enough points
if (nrow(cetacean_data) < 3) {
  stop("Not enough cetacean points (", nrow(cetacean_data), ") to calculate convex hull. Need at least 3.")
}
if (nrow(non_cetacean_data) < 3) {
  stop("Not enough non-cetacean points (", nrow(non_cetacean_data), ") to calculate convex hull. Need at least 3.")
}

# Calculate convex hulls
message("\n", paste(rep("=", 60), collapse = ""))
message("Calculating Convex Hulls")
message(paste(rep("=", 60), collapse = ""))

# Cetacean hull
cetacean_hull_idx <- chull(cetacean_data$x, cetacean_data$y)
cetacean_hull <- cetacean_data[cetacean_hull_idx, ]
cetacean_hull <- rbind(cetacean_hull, cetacean_hull[1, ])  # Close the polygon

# Non-cetacean hull
non_cetacean_hull_idx <- chull(non_cetacean_data$x, non_cetacean_data$y)
non_cetacean_hull <- non_cetacean_data[non_cetacean_hull_idx, ]
non_cetacean_hull <- rbind(non_cetacean_hull, non_cetacean_hull[1, ])  # Close the polygon

# Calculate metrics
cetacean_area <- calculate_hull_area(cetacean_data)
cetacean_perimeter <- calculate_hull_perimeter(cetacean_data)

non_cetacean_area <- calculate_hull_area(non_cetacean_data)
non_cetacean_perimeter <- calculate_hull_perimeter(non_cetacean_data)

# Calculate overlap
overlap_results <- calculate_hull_overlap(cetacean_data, non_cetacean_data)

message("\nCetacean convex hull:")
message("  Area: ", round(cetacean_area, 4))
message("  Perimeter: ", round(cetacean_perimeter, 4))
message("  Number of points: ", nrow(cetacean_data))
message("  Hull vertices: ", length(cetacean_hull_idx))

message("\nNon-cetacean convex hull:")
message("  Area: ", round(non_cetacean_area, 4))
message("  Perimeter: ", round(non_cetacean_perimeter, 4))
message("  Number of points: ", nrow(non_cetacean_data))
message("  Hull vertices: ", length(non_cetacean_hull_idx))

message("\nOverlap:")
message("  Overlap area: ", round(overlap_results$overlap_area, 4))
message("  Overlap as % of cetacean hull: ", round(overlap_results$overlap_percentage_hull1, 2), "%")
message("  Overlap as % of non-cetacean hull: ", round(overlap_results$overlap_percentage_hull2, 2), "%")
message("  Area ratio (cetacean/non-cetacean): ", round(cetacean_area / non_cetacean_area, 4))

# Create visualization
message("\n", paste(rep("=", 60), collapse = ""))
message("Creating Visualization")
message(paste(rep("=", 60), collapse = ""))

# Prepare data for plotting
plot_data <- pca_scores %>%
  filter(group %in% c("Cetacean", "Non_Cetacean"))

# Create plot
p <- ggplot(plot_data, aes(x = x, y = y, color = group, fill = group)) +
  # Plot convex hulls
  geom_polygon(data = cetacean_hull, aes(x = x, y = y), 
               fill = "steelblue", alpha = 0.2, color = "steelblue", linewidth = 1.2) +
  geom_polygon(data = non_cetacean_hull, aes(x = x, y = y), 
               fill = "darkorange", alpha = 0.2, color = "darkorange", linewidth = 1.2) +
  # Plot points
  geom_point(aes(color = group), size = 2, alpha = 0.7) +
  scale_color_manual(
    values = c("Cetacean" = "steelblue", "Non_Cetacean" = "darkorange"),
    labels = c("Cetacean" = "Cetaceans", "Non_Cetacean" = "Non-Cetaceans")
  ) +
  scale_fill_manual(
    values = c("Cetacean" = "steelblue", "Non_Cetacean" = "darkorange"),
    labels = c("Cetacean" = "Cetaceans", "Non_Cetacean" = "Non-Cetaceans")
  ) +
  labs(
    x = paste0(pc_x, " (PCA axis)"),
    y = paste0(pc_y, " (PCA axis)"),
    title = "Convex Hull Analysis: Cetaceans vs Non-Cetaceans",
    color = "Group",
    fill = "Group"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    legend.position = "right",
    panel.grid.minor = element_blank()
  )

# Save plot
plot_file <- paste0(out_prefix, "_convex_hull_", pc_x, "_", pc_y, ".pdf")
ggsave(plot_file, plot = p, width = 10, height = 8)
message("Saved plot: ", plot_file)

# Save summary statistics
summary_data <- data.frame(
  metric = c(
    "cetacean_hull_area",
    "cetacean_hull_perimeter",
    "cetacean_n_points",
    "cetacean_n_hull_vertices",
    "non_cetacean_hull_area",
    "non_cetacean_hull_perimeter",
    "non_cetacean_n_points",
    "non_cetacean_n_hull_vertices",
    "overlap_area",
    "overlap_percentage_cetacean",
    "overlap_percentage_non_cetacean",
    "area_ratio_cetacean_to_non_cetacean"
  ),
  value = c(
    cetacean_area,
    cetacean_perimeter,
    nrow(cetacean_data),
    length(cetacean_hull_idx),
    non_cetacean_area,
    non_cetacean_perimeter,
    nrow(non_cetacean_data),
    length(non_cetacean_hull_idx),
    overlap_results$overlap_area,
    overlap_results$overlap_percentage_hull1,
    overlap_results$overlap_percentage_hull2,
    cetacean_area / non_cetacean_area
  )
)

summary_file <- paste0(out_prefix, "_convex_hull_summary.csv")
write.csv(summary_data, file = summary_file, row.names = FALSE)
message("Saved summary: ", summary_file)

# Save hull coordinates
hull_coords <- list(
  cetacean = cetacean_hull,
  non_cetacean = non_cetacean_hull
)

hull_file <- paste0(out_prefix, "_convex_hull_coordinates.csv")
hull_df <- data.frame(
  group = c(rep("Cetacean", nrow(cetacean_hull)), rep("Non_Cetacean", nrow(non_cetacean_hull))),
  x = c(cetacean_hull$x, non_cetacean_hull$x),
  y = c(cetacean_hull$y, non_cetacean_hull$y)
)
write.csv(hull_df, file = hull_file, row.names = FALSE)
message("Saved hull coordinates: ", hull_file)

message("\n", paste(rep("=", 60), collapse = ""))
message("Analysis Complete!")
message(paste(rep("=", 60), collapse = ""))
message("\nOutput files:")
message("  - ", plot_file)
message("  - ", summary_file)
message("  - ", hull_file)

