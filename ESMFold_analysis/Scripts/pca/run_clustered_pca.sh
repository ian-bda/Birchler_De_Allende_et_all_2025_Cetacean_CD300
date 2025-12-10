#!/bin/bash
#SBATCH --job-name=clustered_pca
#SBATCH --output=/home5/ibirchl/Yoder_Lab/CD300_Project/Cetacean/Evolutionary_Diversification/ESMFold_analysis/Scripts/pca/clustered_pca_%j.out
#SBATCH --error=/home5/ibirchl/Yoder_Lab/CD300_Project/Cetacean/Evolutionary_Diversification/ESMFold_analysis/Scripts/pca/clustered_pca_%j.err
#SBATCH --partition=standard
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8
#SBATCH --time=12:00:00

cd /home5/ibirchl/Yoder_Lab/CD300_Project/Cetacean/Evolutionary_Diversification

# Check for Pillow (PIL) - required for matplotlib to save PDFs
echo "Checking for Pillow (PIL)..."
if python3 -c "from PIL import Image" 2>/dev/null; then
    echo "Pillow is available"
else
    echo "WARNING: Pillow (PIL) not available - plots will be skipped"
    echo "To install Pillow, run: python3 -m pip install --user Pillow"
    echo "Or use a Python environment that has Pillow installed"
fi

# Run clustered PCA analysis
# This divides sequences into 3 groups based on PC1/PC2 thresholds and re-runs PCA within each group
# Group 1: PC1 < 5 AND PC2 < 10
# Group 2: PC1 >= 5 AND PC2 < 20
# Group 3: PC2 >= 20 AND PC1 < -10

python3 ESMFold_analysis/Scripts/pca/clustered_pca_analysis.py \
    ESMFold_analysis/ESMFold_output/pca_results_no_PIGR_single_domain_new_tree_pca_scores.csv \
    ESMFold_analysis/ESMFold_output/structures/esmfold_coordinates_no_PIGR_single_domain.npz \
    ESMFold_analysis/ESMFold_output/structures/coordinate_metadata_no_PIGR_single_domain.csv \
    ESMFold_analysis/ESMFold_output/clustered_pca_no_PIGR_single_domain \
    --min-cluster-size 5

echo ""
echo "All outputs saved to: ESMFold_analysis/ESMFold_output/clustered_pca/"

echo "Clustered PCA analysis complete!"

