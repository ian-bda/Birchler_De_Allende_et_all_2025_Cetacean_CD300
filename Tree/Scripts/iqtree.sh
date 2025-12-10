#!/bin/bash
#SBATCH --job-name=iqtree
#SBATCH --output=iqtree_%j.out
#SBATCH --error=iqtree_%j.err
#SBATCH --partition=gpu
#SBATCH --mem=200G
#SBATCH --cpus-per-task=8
#SBATCH --time=48:00:00

# Input alignment
ALIGNMENT="/home5/ibirchl/Yoder_Lab/CD300_Project/Cetacean/Evolutionary_Diversification/Tree/Tree_with_PIGRs/Arteodactyl_ig_hits_yes_PIGR_single_domain.aln"

# Number of threads (match SLURM cpus-per-task)
THREADS=8

# Output prefix
PREFIX="Arteodactyl_ig_hits_yes_PIGR_single_domain_iqtree"

# Run IQ-TREE
echo "Running IQ-TREE on $ALIGNMENT using $THREADS threads..."

/home5/ibirchl/Bioinformatics_tools/iqtree-2.2.2.6-Linux/bin/iqtree2 -s "$ALIGNMENT" \
        -nt $THREADS \
        -m MFP \
        -bb 1000 \
        -alrt 1000 \
        -pre "$PREFIX" \
        -quiet

echo "IQ-TREE analysis completed. Output files prefixed with $PREFIX"
