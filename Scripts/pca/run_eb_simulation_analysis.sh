#!/bin/bash
#SBATCH --job-name=eb_sim_analysis
#SBATCH --output=/home5/ibirchl/Yoder_Lab/CD300_Project/Cetacean/Evolutionary_Diversification/ESMFold_analysis/Scripts/pca/eb_simulation_analysis_%j.out
#SBATCH --error=/home5/ibirchl/Yoder_Lab/CD300_Project/Cetacean/Evolutionary_Diversification/ESMFold_analysis/Scripts/pca/eb_simulation_analysis_%j.err
#SBATCH --partition=standard
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --time=24:00:00

cd /home5/ibirchl/Yoder_Lab/CD300_Project/Cetacean/Evolutionary_Diversification

# Run EB simulation analysis
# This script:
# 1. Fits EB models to PC1 and PC2
# 2. Simulates PC1 and PC2 under EB model 1000x
# 3. Calculates convex hull metrics for each simulation
# 4. Compares observed metrics to simulated distribution
# 5. Calculates p-values to test if observed is outside EB model predictions

/opt/R/4.3.3/bin/Rscript ESMFold_analysis/Scripts/pca/eb_simulation_analysis.R \
    ESMFold_analysis/ESMFold_output/pca_results_no_PIGR_single_domain_new_tree_pca_scores.csv \
    Tree/Tree_for_Analysis/Arteodactyl_ig_hits_no_PIGR_single_domain_iqtree_ultrametric.nwk \
    ESMFold_analysis/ESMFold_output/eb_simulation_no_PIGR_single_domain \
    1000 \
    PC1 \
    PC2

echo "EB simulation analysis complete!"

