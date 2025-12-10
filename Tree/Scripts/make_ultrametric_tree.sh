#!/bin/bash
#SBATCH --job-name=ultrametric_tree
#SBATCH --output=/home5/ibirchl/Yoder_Lab/CD300_Project/Cetacean/Evolutionary_Diversification/Tree/Scripts/make_ultrametric_tree_%j.out
#SBATCH --error=/home5/ibirchl/Yoder_Lab/CD300_Project/Cetacean/Evolutionary_Diversification/Tree/Scripts/make_ultrametric_tree_%j.err
#SBATCH --partition=standard
#SBATCH --mem=8G
#SBATCH --cpus-per-task=2
#SBATCH --time=2:00:00

# Change to project directory
cd /home5/ibirchl/Yoder_Lab/CD300_Project/Cetacean/Evolutionary_Diversification

# Ensure output directory exists
mkdir -p Tree/Tree_for_Analysis

# Run R script with improved accuracy methods
# Method hierarchy: chronos() -> force.ultrametric(nnls) -> force.ultrametric(extend)
/opt/R/4.3.3/bin/Rscript Tree/scripts/make_ultrametric_tree.R \
    Tree/Tree_for_Analysis/Arteodactyl_ig_hits_no_PIGR_single_domain_iqtree.treefile \
    Tree/Tree_for_Analysis/Arteodactyl_ig_hits_no_PIGR_single_domain_iqtree_ultrametric.nwk

echo "Ultrametric tree creation complete!"

