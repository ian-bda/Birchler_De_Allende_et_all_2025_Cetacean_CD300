#!/bin/bash
#SBATCH --job-name=esmfold_prediction
#SBATCH --output=esmfold_%j.out
#SBATCH --error=esmfold_%j.err
#SBATCH --partition=gpu
#SBATCH --mem=200G
#SBATCH --cpus-per-task=8
#SBATCH --time=48:00:00

# Load required modules (if available)
# module load cuda/11.8

# Use system Python
export PATH="/home5/ibirchl/miniforge3/bin:$PATH"

# Install required packages if not already installed
# pip install torch biotite fair-esm

# Navigate to the working directory
cd /home5/ibirchl/Yoder_Lab/CD300_Project/Cetacean/Evolutionary_Diversification/ESMFold_analysis

# Run ESMFold prediction
python Scripts/esmfold/esmfold_predict.py \
    /home5/ibirchl/Yoder_Lab/CD300_Project/Cetacean/Evolutionary_Diversification/ESMFold_analysis/Arteodactyl_ig_hits.fasta \
    --output-dir ESMFold_output/structures \
    --gpu

echo "ESMFold prediction completed!"