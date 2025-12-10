#!/bin/bash
#SBATCH --job-name=famsa_alignment
#SBATCH --output=famsa_alignment_%j.out
#SBATCH --error=famsa_alignment_%j.err
#SBATCH --partition=gpu
#SBATCH --mem=200G
#SBATCH --cpus-per-task=8
#SBATCH --time=48:00:00

# Set input/output file names
INPUT="/home5/ibirchl/Yoder_Lab/CD300_Project/Cetacean/Evolutionary_Diversification/Data/Arteodactyl_ig_hits_yes_PIGR_single_domain.fasta"
OUTPUT="/home5/ibirchl/Yoder_Lab/CD300_Project/Cetacean/Evolutionary_Diversification/Tree/Tree_with_PIGRs/Arteodactyl_ig_hits_yes_PIGR_single_domain.aln"

# Number of threads
THREADS=8

# Run FAMSA
echo "Aligning $INPUT with $THREADS threads..."
/home5/ibirchl/Bioinformatics_tools/FAMSA/famsa -t $THREADS "$INPUT" "$OUTPUT"

echo "Alignment saved to $OUTPUT"
