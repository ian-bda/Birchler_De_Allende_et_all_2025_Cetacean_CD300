#!/bin/bash
#SBATCH --job-name=evolutionary_models
#SBATCH --output=/home5/ibirchl/Yoder_Lab/CD300_Project/Cetacean/Evolutionary_Diversification/ESMFold_analysis/ESMFold_output/evolutionary_models_%j.out
#SBATCH --error=/home5/ibirchl/Yoder_Lab/CD300_Project/Cetacean/Evolutionary_Diversification/ESMFold_analysis/ESMFold_output/evolutionary_models_%j.err
#SBATCH --partition=standard
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --time=12:00:00

# Change to project directory
cd /home5/ibirchl/Yoder_Lab/CD300_Project/Cetacean/Evolutionary_Diversification

# Set PC components (can be overridden by environment variables)
export PC_COMPONENTS=${PC_COMPONENTS:-"PC1,PC2"}

# Load R module if needed (adjust based on your cluster)
# module load R/4.3.3

# Use full path to Rscript
/opt/R/4.3.3/bin/Rscript ESMFold_analysis/Scripts/pca/evolutionary_models.R \
    ESMFold_analysis/ESMFold_output/pca_results_no_PIGR_single_domain_new_tree_pca_scores.csv \
    Tree/Tree_for_Analysis/Arteodactyl_ig_hits_no_PIGR_single_domain_iqtree_ultrametric.nwk \
    ESMFold_analysis/ESMFold_output/evolutionary_models_no_PIGR_single_domain_ultrametric \
    ${PC_COMPONENTS}

echo "Evolutionary models analysis complete!"

