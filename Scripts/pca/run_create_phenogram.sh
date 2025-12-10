#!/bin/bash
#SBATCH --job-name=phenogram
#SBATCH --output=/home5/ibirchl/Yoder_Lab/CD300_Project/Cetacean/Evolutionary_Diversification/ESMFold_analysis/Scripts/pca/create_phenogram_%j.out
#SBATCH --error=/home5/ibirchl/Yoder_Lab/CD300_Project/Cetacean/Evolutionary_Diversification/ESMFold_analysis/Scripts/pca/create_phenogram_%j.err
#SBATCH --partition=standard
#SBATCH --mem=8G
#SBATCH --cpus-per-task=2
#SBATCH --time=1:00:00

cd /home5/ibirchl/Yoder_Lab/CD300_Project/Cetacean/Evolutionary_Diversification

# Create phenograms for PC1 and PC2
# Phenograms show trait evolution along the phylogeny with time on x-axis and trait value on y-axis

/opt/R/4.3.3/bin/Rscript ESMFold_analysis/Scripts/pca/create_phenogram.R \
    ESMFold_analysis/ESMFold_output/pca_results_no_PIGR_single_domain_new_tree_pca_scores.csv \
    Tree/Tree_for_Analysis/Arteodactyl_ig_hits_no_PIGR_single_domain_iqtree_ultrametric.nwk \
    ESMFold_analysis/ESMFold_output/phenogram_no_PIGR_single_domain \
    PC1,PC2

echo "Phenogram creation complete!"

