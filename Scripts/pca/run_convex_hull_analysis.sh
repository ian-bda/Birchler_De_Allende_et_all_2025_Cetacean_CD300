#!/bin/bash
#SBATCH --job-name=convex_hull
#SBATCH --output=/home5/ibirchl/Yoder_Lab/CD300_Project/Cetacean/Evolutionary_Diversification/ESMFold_analysis/ESMFold_output/convex_hull_%j.out
#SBATCH --error=/home5/ibirchl/Yoder_Lab/CD300_Project/Cetacean/Evolutionary_Diversification/ESMFold_analysis/ESMFold_output/convex_hull_%j.err
#SBATCH --partition=standard
#SBATCH --mem=8G
#SBATCH --cpus-per-task=2
#SBATCH --time=1:00:00

# Change to project directory
cd /home5/ibirchl/Yoder_Lab/CD300_Project/Cetacean/Evolutionary_Diversification

# Run convex hull analysis
/opt/R/4.3.3/bin/Rscript ESMFold_analysis/Scripts/pca/convex_hull_analysis.R \
    ESMFold_analysis/ESMFold_output/pca_results_no_PIGR_single_domain_new_tree_pca_scores.csv \
    ESMFold_analysis/ESMFold_output/convex_hull_no_PIGR_single_domain \
    PC1 PC2 TRUE

echo "Convex hull analysis complete!"

