#!/bin/bash
#SBATCH --job-name=disparity_analysis
#SBATCH --output=/home5/ibirchl/Yoder_Lab/CD300_Project/Cetacean/Evolutionary_Diversification/ESMFold_analysis/ESMFold_output/disparity_analysis_%j.out
#SBATCH --error=/home5/ibirchl/Yoder_Lab/CD300_Project/Cetacean/Evolutionary_Diversification/ESMFold_analysis/ESMFold_output/disparity_analysis_%j.err
#SBATCH --partition=standard
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --time=6:00:00

# Change to project directory
cd /home5/ibirchl/Yoder_Lab/CD300_Project/Cetacean/Evolutionary_Diversification

# Run disparity analysis with ultrametric tree
python3 ESMFold_analysis/Scripts/pca/disparity_analysis.py \
    ESMFold_analysis/ESMFold_output/pca_results_no_PIGR_single_domain_new_tree_pca_scores.csv \
    Tree/Tree_for_Analysis/Arteodactyl_ig_hits_no_PIGR_single_domain_iqtree_ultrametric.nwk \
    --output-prefix ESMFold_analysis/ESMFold_output/disparity_no_PIGR_single_domain_ultrametric \
    --exclude-humans

echo "Disparity analysis complete!"

