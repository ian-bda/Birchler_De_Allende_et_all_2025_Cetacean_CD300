#!/usr/bin/env python3
"""
Disparity Analysis Script

Analyzes evolutionary disparity by comparing:
1. Morphological distances (Euclidean distance in PCA space)
2. Phylogenetic distances (branch length distances in tree)
3. Tests against Brownian Motion expectations

The key question: Are distantly related taxa more morphologically different
than closely related taxa? Does morphological evolution follow Brownian motion?

Usage:
    python disparity_analysis.py <pca_scores_csv> <tree_file> [options]
"""

import argparse
import os
import sys
import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist, squareform
from scipy.stats import pearsonr, spearmanr, linregress
from scipy import stats
from Bio import Phylo
import warnings
warnings.filterwarnings('ignore')


def is_human(seq_id):
    """Check if sequence is from human (Homo sapiens)."""
    if pd.isna(seq_id):
        return False
    seq_id_str = str(seq_id)
    # Check if it starts with "Homo" or contains "Homo_sapiens"
    return seq_id_str.startswith("Homo") or "Homo_sapiens" in seq_id_str


def load_pca_scores(pca_scores_file, exclude_humans=False):
    """Load PCA scores from CSV file."""
    print("\n" + "="*60)
    print("Loading PCA Scores")
    print("="*60)
    
    df = pd.read_csv(pca_scores_file)
    print(f"Loaded {len(df)} sequences")
    
    # Determine ID column
    if 'original_sequence_id' in df.columns:
        id_col = 'original_sequence_id'
    elif 'clean_id' in df.columns:
        id_col = 'clean_id'
    else:
        id_col = df.columns[0]
    
    print(f"Using '{id_col}' as sequence identifier")
    
    # Filter out humans if requested
    if exclude_humans:
        before_count = len(df)
        df = df[~df[id_col].apply(is_human)]
        removed_count = before_count - len(df)
        if removed_count > 0:
            print(f"Excluded {removed_count} human sequences")
        print(f"Using {len(df)} sequences (excluding humans)")
    
    # Get PC columns
    pc_cols = [col for col in df.columns if col.startswith('PC')]
    print(f"Found {len(pc_cols)} principal components")
    
    # Extract PC scores matrix
    pc_matrix = df[pc_cols].values
    sequence_ids = df[id_col].values
    
    return pc_matrix, sequence_ids, df, pc_cols


def load_phylogenetic_tree(tree_file):
    """Load phylogenetic tree from Newick file."""
    print("\n" + "="*60)
    print("Loading Phylogenetic Tree")
    print("="*60)
    
    tree = Phylo.read(tree_file, 'newick')
    
    # Get all tip names
    tip_names = [tip.name for tip in tree.get_terminals()]
    print(f"Tree has {len(tip_names)} tips")
    
    return tree, tip_names


def calculate_morphological_distances(pc_matrix, sequence_ids):
    """
    Calculate pairwise Euclidean distances in PCA space.
    
    This represents morphological/structural differences between sequences.
    """
    print("\n" + "="*60)
    print("Calculating Morphological Distances")
    print("="*60)
    
    # Calculate pairwise Euclidean distances
    # pdist computes distances between all pairs, returns condensed matrix
    morph_distances = pdist(pc_matrix, metric='euclidean')
    
    # Convert to square matrix for easier indexing
    morph_dist_matrix = squareform(morph_distances)
    
    print(f"Calculated {len(morph_distances)} pairwise distances")
    print(f"  Min distance: {morph_distances.min():.4f}")
    print(f"  Max distance: {morph_distances.max():.4f}")
    print(f"  Mean distance: {morph_distances.mean():.4f}")
    
    return morph_distances, morph_dist_matrix, sequence_ids


def calculate_phylogenetic_distances(tree, sequence_ids):
    """
    Calculate pairwise phylogenetic distances from tree.
    
    This represents evolutionary time/divergence between sequences.
    """
    print("\n" + "="*60)
    print("Calculating Phylogenetic Distances")
    print("="*60)
    
    # Get all tip names and create mapping
    tip_names = [tip.name for tip in tree.get_terminals()]
    
    # Clean sequence IDs to match tree tip labels (remove spaces, etc.)
    # Tree tips might have spaces replaced with underscores
    sequence_ids_clean = [s.replace(' ', '_') for s in sequence_ids]
    
    # Find common sequences
    common_seqs = set(sequence_ids_clean) & set(tip_names)
    print(f"Found {len(common_seqs)} common sequences between tree and PCA data")
    
    if len(common_seqs) < 3:
        raise ValueError(f"Too few common sequences ({len(common_seqs)}). Need at least 3.")
    
    # Create distance dictionary
    # Phylo can calculate distances between any two tips
    phylo_distances_list = []
    phylo_pairs = []
    
    # Calculate pairwise distances for common sequences only
    common_list = sorted(list(common_seqs))
    n_common = len(common_list)
    
    print(f"Calculating phylogenetic distances for {n_common} sequences...")
    print(f"  This will compute {n_common * (n_common - 1) // 2} pairwise distances")
    
    for i in range(n_common):
        if (i + 1) % 50 == 0:
            print(f"  Processed {i+1}/{n_common} sequences...")
        for j in range(i + 1, n_common):
            seq1 = common_list[i]
            seq2 = common_list[j]
            
            # Find tips in tree
            tip1 = None
            tip2 = None
            for tip in tree.get_terminals():
                if tip.name == seq1:
                    tip1 = tip
                if tip.name == seq2:
                    tip2 = tip
                    if tip1 is not None:
                        break
            
            if tip1 is not None and tip2 is not None:
                # Calculate distance (sum of branch lengths along path)
                dist = tree.distance(tip1, tip2)
                phylo_distances_list.append(dist)
                phylo_pairs.append((seq1, seq2))
    
    print(f"Calculated {len(phylo_distances_list)} pairwise phylogenetic distances")
    if len(phylo_distances_list) > 0:
        print(f"  Min distance: {min(phylo_distances_list):.6f}")
        print(f"  Max distance: {max(phylo_distances_list):.6f}")
        print(f"  Mean distance: {np.mean(phylo_distances_list):.6f}")
    
    return np.array(phylo_distances_list), phylo_pairs, common_list


def match_distances(morph_dist_matrix, morph_sequence_ids, 
                   phylo_distances, phylo_pairs, common_seqs):
    """
    Match morphological and phylogenetic distances for common sequences.
    
    Returns arrays of matched distances.
    """
    print("\n" + "="*60)
    print("Matching Morphological and Phylogenetic Distances")
    print("="*60)
    
    # Clean sequence IDs
    morph_ids_clean = [s.replace(' ', '_') for s in morph_sequence_ids]
    
    # Create mapping from sequence ID to index in morph matrix
    morph_id_to_idx = {seq_id: idx for idx, seq_id in enumerate(morph_ids_clean)}
    
    # Match distances
    matched_morph = []
    matched_phylo = []
    matched_pairs = []
    
    for (seq1, seq2), phylo_dist in zip(phylo_pairs, phylo_distances):
        if seq1 in morph_id_to_idx and seq2 in morph_id_to_idx:
            idx1 = morph_id_to_idx[seq1]
            idx2 = morph_id_to_idx[seq2]
            
            morph_dist = morph_dist_matrix[idx1, idx2]
            
            matched_morph.append(morph_dist)
            matched_phylo.append(phylo_dist)
            matched_pairs.append((seq1, seq2))
    
    matched_morph = np.array(matched_morph)
    matched_phylo = np.array(matched_phylo)
    
    print(f"Matched {len(matched_morph)} pairs of distances")
    print(f"  Morphological distances: min={matched_morph.min():.4f}, max={matched_morph.max():.4f}, mean={matched_morph.mean():.4f}")
    print(f"  Phylogenetic distances: min={matched_phylo.min():.6f}, max={matched_phylo.max():.6f}, mean={matched_phylo.mean():.6f}")
    
    return matched_morph, matched_phylo, matched_pairs


def test_brownian_motion(morph_dist, phylo_dist):
    """
    Test if morphological evolution follows Brownian motion.
    
    Under Brownian motion:
    - Morphological distance should increase proportionally with phylogenetic distance
    - Linear relationship: morph_dist ~ alpha * phylo_dist
    - Correlation should be positive and significant
    """
    print("\n" + "="*60)
    print("Testing Brownian Motion Model")
    print("="*60)
    
    # Remove zero phylogenetic distances (same species, no divergence)
    valid_mask = phylo_dist > 0
    if np.sum(valid_mask) < len(phylo_dist):
        print(f"  Warning: {len(phylo_dist) - np.sum(valid_mask)} pairs have zero phylogenetic distance")
        print(f"  Using {np.sum(valid_mask)} pairs for analysis")
    
    morph_valid = morph_dist[valid_mask]
    phylo_valid = phylo_dist[valid_mask]
    
    # Linear regression
    slope, intercept, r_value, p_value, std_err = linregress(phylo_valid, morph_valid)
    
    print(f"\nLinear Regression Results:")
    print(f"  Slope (alpha): {slope:.6f} +/- {std_err:.6f}")
    print(f"  Intercept: {intercept:.6f}")
    print(f"  R-squared: {r_value**2:.4f}")
    print(f"  P-value: {p_value:.2e}")
    print(f"  Correlation: {r_value:.4f}")
    
    # Correlation tests
    pearson_r, pearson_p = pearsonr(phylo_valid, morph_valid)
    spearman_r, spearman_p = spearmanr(phylo_valid, morph_valid)
    
    print(f"\nCorrelation Tests:")
    print(f"  Pearson correlation: r = {pearson_r:.4f}, p = {pearson_p:.2e}")
    print(f"  Spearman correlation: rho = {spearman_r:.4f}, p = {spearman_p:.2e}")
    
    # Interpretation
    print(f"\nInterpretation:")
    if p_value < 0.05:
        print(f"  [+] Significant relationship (p < 0.05)")
        if r_value > 0:
            print(f"  [+] Positive correlation: more distantly related = more morphologically different")
        else:
            print(f"  [-] Negative correlation: more distantly related = LESS morphologically different (unexpected!)")
    else:
        print(f"  [-] No significant relationship (p >= 0.05)")
        print(f"    Morphological differences are NOT proportional to phylogenetic distance")
    
    if r_value**2 > 0.5:
        print(f"  [+] Strong relationship (R-squared > 0.5): Brownian motion explains >50% of variance")
    elif r_value**2 > 0.3:
        print(f"  [~] Moderate relationship (R-squared = {r_value**2:.2f}): Brownian motion explains ~{r_value**2*100:.0f}% of variance")
    else:
        print(f"  [-] Weak relationship (R-squared < 0.3): Brownian motion explains <30% of variance")
        print(f"    Other factors (selection, drift, constraint) may be important")
    
    return {
        'slope': slope,
        'intercept': intercept,
        'r_squared': r_value**2,
        'p_value': p_value,
        'pearson_r': pearson_r,
        'pearson_p': pearson_p,
        'spearman_r': spearman_r,
        'spearman_p': spearman_p,
        'n_pairs': len(morph_valid)
    }


def save_results(morph_dist, phylo_dist, results, matched_pairs, output_prefix):
    """Save analysis results to CSV files."""
    print("\n" + "="*60)
    print("Saving Results")
    print("="*60)
    
    # Filter valid pairs
    valid_mask = phylo_dist > 0
    morph_valid = morph_dist[valid_mask]
    phylo_valid = phylo_dist[valid_mask]
    pairs_valid = [matched_pairs[i] for i in range(len(matched_pairs)) if valid_mask[i]]
    
    # Create dataframe
    df = pd.DataFrame({
        'sequence_1': [p[0] for p in pairs_valid],
        'sequence_2': [p[1] for p in pairs_valid],
        'morphological_distance': morph_valid,
        'phylogenetic_distance': phylo_valid
    })
    
    # Save pairwise distances
    output_file = f"{output_prefix}_pairwise_distances.csv"
    df.to_csv(output_file, index=False)
    print(f"Saved: {output_file}")
    
    # Save summary statistics
    summary = pd.DataFrame({
        'statistic': [
            'n_pairs',
            'mean_morphological_distance',
            'mean_phylogenetic_distance',
            'slope',
            'intercept',
            'r_squared',
            'pearson_r',
            'pearson_p',
            'spearman_r',
            'spearman_p'
        ],
        'value': [
            results['n_pairs'],
            morph_valid.mean(),
            phylo_valid.mean(),
            results['slope'],
            results['intercept'],
            results['r_squared'],
            results['pearson_r'],
            results['pearson_p'],
            results['spearman_r'],
            results['spearman_p']
        ]
    })
    
    summary_file = f"{output_prefix}_disparity_summary.csv"
    summary.to_csv(summary_file, index=False)
    print(f"Saved: {summary_file}")


def main():
    parser = argparse.ArgumentParser(
        description="Analyze evolutionary disparity: compare morphological and phylogenetic distances",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic usage
  python disparity_analysis.py pca_scores.csv tree.nwk
  
  # With custom output prefix
  python disparity_analysis.py pca_scores.csv tree.nwk -o disparity_results
  
  # Use full paths
  python disparity_analysis.py /path/to/pca_scores.csv /path/to/tree.nwk
        """
    )
    
    parser.add_argument(
        "pca_scores_csv",
        help="Input CSV file with PCA scores (from structural_pca.py)"
    )
    
    parser.add_argument(
        "tree_file",
        help="Input Newick tree file (from IQ-TREE or similar)"
    )
    
    parser.add_argument(
        "-o", "--output-prefix",
        default="disparity",
        help="Output file prefix (default: disparity)"
    )
    
    parser.add_argument(
        "--pc-components",
        type=int,
        nargs='+',
        help="Which PC components to use for morphological distance (default: all)"
    )
    
    parser.add_argument(
        "--exclude-humans",
        action="store_true",
        help="Exclude human (Homo) sequences from analysis"
    )
    
    args = parser.parse_args()
    
    # Check input files
    if not os.path.exists(args.pca_scores_csv):
        print(f"Error: PCA scores file '{args.pca_scores_csv}' not found")
        sys.exit(1)
    
    if not os.path.exists(args.tree_file):
        print(f"Error: Tree file '{args.tree_file}' not found")
        sys.exit(1)
    
    print("="*60)
    print("Evolutionary Disparity Analysis")
    print("="*60)
    print(f"PCA scores file: {args.pca_scores_csv}")
    print(f"Tree file: {args.tree_file}")
    print(f"Output prefix: {args.output_prefix}")
    if args.exclude_humans:
        print("Excluding human sequences from analysis")
    
    # Load data
    pc_matrix, sequence_ids, pca_df, pc_cols = load_pca_scores(args.pca_scores_csv, exclude_humans=args.exclude_humans)
    
    # Filter PC components if specified
    if args.pc_components:
        pc_indices = [i for i, col in enumerate(pc_cols) if int(col.replace('PC', '')) in args.pc_components]
        if pc_indices:
            pc_matrix = pc_matrix[:, pc_indices]
            print(f"Using PC components: {args.pc_components}")
        else:
            print(f"Warning: No valid PC components found in {args.pc_components}, using all")
    
    tree, tip_names = load_phylogenetic_tree(args.tree_file)
    
    # Filter tree to exclude human tips if requested
    if args.exclude_humans:
        human_tips = [tip for tip in tip_names if is_human(tip)]
        if human_tips:
            print(f"\nRemoving {len(human_tips)} human tips from tree...")
            # Get all non-human tips
            non_human_tips = [tip for tip in tip_names if not is_human(tip)]
            # Prune tree (this is a simplified approach - Bio.Phylo doesn't have easy pruning)
            # We'll filter during distance calculation instead
            print(f"Tree will be filtered during distance calculation")
    
    # Calculate distances
    morph_distances, morph_dist_matrix, morph_ids = calculate_morphological_distances(pc_matrix, sequence_ids)
    phylo_distances, phylo_pairs, common_seqs = calculate_phylogenetic_distances(tree, sequence_ids)
    
    # Match distances
    matched_morph, matched_phylo, matched_pairs = match_distances(
        morph_dist_matrix, sequence_ids, phylo_distances, phylo_pairs, common_seqs
    )
    
    # Test Brownian motion
    results = test_brownian_motion(matched_morph, matched_phylo)
    
    # Save results (CSV only, no plots)
    save_results(matched_morph, matched_phylo, results, matched_pairs, args.output_prefix)
    
    print("\n" + "="*60)
    print("Analysis Complete!")
    print("="*60)
    print(f"\nOutput files:")
    print(f"  - {args.output_prefix}_pairwise_distances.csv")
    print(f"  - {args.output_prefix}_disparity_summary.csv")


if __name__ == "__main__":
    main()

