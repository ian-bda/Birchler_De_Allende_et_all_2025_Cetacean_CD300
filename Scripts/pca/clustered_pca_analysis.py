#!/usr/bin/env python3
"""
Clustered PCA Analysis Script

Divides sequences into 3 groups based on PC1/PC2 thresholds and re-runs PCA within each group.
1. Loads existing PCA scores
2. Removes humans
3. Divides PC1/PC2 space into 3 groups using fixed thresholds:
   - Group 1: PC1 < 5 AND PC2 < 10
   - Group 2: PC1 >= 5 AND PC2 < 20
   - Group 3: PC2 >= 20 AND PC1 < -10
4. Re-runs structural PCA for each group separately
5. Creates PC1 vs PC2 plots colored by species group (matching phylomorphospace_colored.R)
6. Checks if groups represent single-copy clades

Usage:
    python clustered_pca_analysis.py <pca_scores_csv> <npz_file> <metadata_csv> <output_prefix> [--min-cluster-size N]

The script will attempt to create plots during analysis. If matplotlib/PIL is not available,
it will try to create plots from saved CSV files at the end (if matplotlib becomes available).
"""

import argparse
import os
import sys
import numpy as np
import pandas as pd
from scipy.spatial.transform import Rotation
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

# Import matplotlib BEFORE importing from structural_pca.py
# This ensures matplotlib is available when structural_pca.py imports it
try:
    import matplotlib
    matplotlib.use('Agg')  # Use non-interactive backend
    import matplotlib.pyplot as plt
    # Test if matplotlib can actually create a figure (checks PIL dependency)
    try:
        fig = plt.figure()
        plt.close(fig)
        HAS_MATPLOTLIB = True
    except Exception as e:
        HAS_MATPLOTLIB = False
        print(f"Warning: matplotlib available but cannot create plots (likely missing PIL/Pillow): {e}")
        print("  Plots will be skipped.")
except (ImportError, ModuleNotFoundError) as e:
    HAS_MATPLOTLIB = False
    print(f"Warning: matplotlib not available: {e}")
    print("  Plots will be skipped.")
    # Create dummy plt object
    class DummyPlot:
        def __getattr__(self, name):
            def dummy(*args, **kwargs):
                pass
            return dummy
    plt = DummyPlot()

try:
    import seaborn as sns
    HAS_SEABORN = True
except ImportError:
    HAS_SEABORN = False
    # Don't print warning for seaborn - it's optional

# Import functions from structural_pca.py (now matplotlib is already imported/handled)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from structural_pca import (
    kabsch_alignment,
    align_to_reference,
    prepare_coordinate_matrix,
    standardize_coordinates,
    perform_pca,
    save_results,
    plot_results
)


def is_human(seq_id):
    """Check if sequence ID is human."""
    seq_id_str = str(seq_id)
    return seq_id_str.startswith("Homo") or "Homo_sapiens" in seq_id_str


def classify_species_group(seq_id):
    """
    Classify species group from sequence ID (same as phylomorphospace_colored.R).
    Returns: 'Cetacean', 'Terrestrial_Artiodactyl', or 'Other'
    """
    seq_id_str = str(seq_id)
    
    # Extract genus (first part before underscore or pipe)
    if "_" in seq_id_str:
        genus = seq_id_str.split("_")[0]
    elif "|" in seq_id_str:
        genus = seq_id_str.split("|")[0]
    else:
        genus = seq_id_str
    
    # Cetacean genera
    cetacean_genera = [
        "Balaenoptera", "Delphinapterus", "Delphinus", "Eschrichtius", 
        "Eubalaena", "Globicephala", "Kogia", "Lagenorhynchus", 
        "Lipotes", "Monodon", "Neophocaena", "Orcinus", "Phocoena", 
        "Physeter", "Pseudorca", "Tursiops", "Sousa", "Mesoplodon"
    ]
    
    # Terrestrial Artiodactyl genera
    artiodactyl_genera = [
        "Bison", "Bos", "Bubalus", "Budorcas", "Camelus", "Capra", 
        "Capricornis", "Cervus", "Dama", "Hippopotamus", "Moschus", 
        "Muntiacus", "Odocoileus", "Oryx", "Ovibos", "Ovis", 
        "Pantholops", "Phacochoerus", "Rangifer", "Sus", "Vicugna", "Homo"
    ]
    
    if genus in cetacean_genera:
        return "Cetacean"
    elif genus in artiodactyl_genera:
        return "Terrestrial_Artiodactyl"
    else:
        return "Other"


def get_species_group_color(group):
    """
    Get color for species group (same as phylomorphospace_colored.R).
    Returns color name or hex code.
    """
    color_map = {
        "Cetacean": "steelblue",
        "Terrestrial_Artiodactyl": "darkorange",
        "Other": "gray50"
    }
    return color_map.get(group, "gray50")


def create_plot_from_csv(pca_scores_file, output_file, cluster_id):
    """
    Create PC1 vs PC2 scatter plot for a cluster from saved CSV file.
    This is a fallback function if plots couldn't be created during analysis.
    """
    if not HAS_MATPLOTLIB:
        return False
    
    try:
        # Read PCA scores
        df = pd.read_csv(pca_scores_file)
        
        # Check for PC1 and PC2 columns
        if 'PC1' not in df.columns or 'PC2' not in df.columns:
            print(f"  Error: PC1 or PC2 not found in {pca_scores_file}")
            return False
        
        # Determine sequence ID column
        id_col = None
        for col in ['clean_id', 'original_sequence_id', 'sequence_id']:
            if col in df.columns:
                id_col = col
                break
        
        if id_col is None:
            # Use first column as ID
            id_col = df.columns[0]
        
        # Classify species groups for coloring (same as phylomorphospace_colored.R)
        df['species_group'] = df[id_col].apply(classify_species_group)
        
        # Calculate explained variance if available
        variance_file = pca_scores_file.replace('_pca_scores.csv', '_explained_variance.csv')
        pc1_var = 0
        pc2_var = 0
        
        if os.path.exists(variance_file):
            var_df = pd.read_csv(variance_file)
            if len(var_df) > 0:
                pc1_var = var_df.iloc[0]['Explained_Variance_Ratio'] * 100
            if len(var_df) > 1:
                pc2_var = var_df.iloc[1]['Explained_Variance_Ratio'] * 100
        
        # Color map (same as phylomorphospace_colored.R)
        group_colors_map = {
            "Cetacean": "steelblue",
            "Terrestrial_Artiodactyl": "darkorange",
            "Other": "gray50"
        }
        
        # Create plot
        fig = plt.figure(figsize=(10, 8))
        
        # Plot each group separately with appropriate colors
        unique_groups = df['species_group'].unique()
        for group in unique_groups:
            group_mask = df['species_group'] == group
            color = group_colors_map.get(group, "gray50")
            # Use alpha=0.7 to match phylomorphospace_colored.R
            plt.scatter(df.loc[group_mask, 'PC1'], df.loc[group_mask, 'PC2'], 
                       c=color, alpha=0.7, s=50, label=group)
        
        plt.xlabel(f'PC1 ({pc1_var:.1f}% variance)', fontsize=12)
        plt.ylabel(f'PC2 ({pc2_var:.1f}% variance)', fontsize=12)
        plt.title(f'Group {cluster_id + 1}: PC1 vs PC2 Structural Morphospace', fontsize=14)
        plt.legend(loc='best', fontsize=9)
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close(fig)
        
        return True
    except Exception as e:
        print(f"  Error creating plot from CSV: {e}")
        if HAS_MATPLOTLIB:
            try:
                plt.close('all')
            except:
                pass
        return False


def assign_clusters(pc1, pc2):
    """
    Assign sequences to 3 clusters based on PC1 and PC2.
    
    Group 1: PC1 < 5 AND PC2 < 10 (left of PC1=5, below PC2=10)
    Group 2: PC1 >= 5 AND PC2 < 20 (right of PC1=5, below PC2=20)
    Group 3: PC2 >= 20 AND PC1 < -10 (above PC2=20, below PC1=-10)
    
    Parameters:
    -----------
    pc1 : array
        PC1 values
    pc2 : array
        PC2 values
    
    Returns:
    --------
    clusters : array, shape (n_samples,)
        Cluster assignments (0, 1, 2, or -1 for unassigned)
    """
    clusters = np.full(len(pc1), -1, dtype=int)  # -1 means unassigned
    
    # Group 1: PC1 < 5 AND PC2 < 10
    clusters[(pc1 < 5.0) & (pc2 < 10.0)] = 0
    
    # Group 2: PC1 >= 5 AND PC2 < 20
    clusters[(pc1 >= 5.0) & (pc2 < 20.0)] = 1
    
    # Group 3: PC2 >= 20 AND PC1 < -10
    clusters[(pc2 >= 20.0) & (pc1 < -10.0)] = 2
    
    return clusters


def check_single_copy_clade(sequence_ids, tree_file=None):
    """
    Check if a cluster represents a single-copy clade.
    This is a simplified check - looks for taxonomic patterns.
    """
    # Extract genera from sequence IDs
    genera = []
    for seq_id in sequence_ids:
        seq_str = str(seq_id)
        if "_" in seq_str:
            genus = seq_str.split("_")[0]
        elif "|" in seq_str:
            genus = seq_str.split("|")[0]
        else:
            genus = seq_str
        genera.append(genus)
    
    unique_genera = set(genera)
    n_genera = len(unique_genera)
    n_sequences = len(sequence_ids)
    
    # If all sequences are from the same genus, likely single-copy
    is_single_copy = (n_genera == 1) or (n_genera / n_sequences < 0.3)
    
    return {
        'is_single_copy': is_single_copy,
        'n_genera': n_genera,
        'n_sequences': n_sequences,
        'unique_genera': list(unique_genera),
        'genus_diversity': n_genera / n_sequences if n_sequences > 0 else 0
    }


def main():
    parser = argparse.ArgumentParser(
        description="Perform clustered PCA analysis",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    parser.add_argument(
        "pca_scores_csv",
        help="Input CSV file with existing PCA scores"
    )
    
    parser.add_argument(
        "npz_file",
        help="Input NPZ file containing coordinate arrays"
    )
    
    parser.add_argument(
        "metadata_csv",
        help="Input CSV file with coordinate metadata"
    )
    
    parser.add_argument(
        "output_prefix",
        help="Output file prefix"
    )
    
    parser.add_argument(
        "--min-cluster-size",
        type=int,
        default=5,
        help="Minimum number of sequences per cluster (default: 5)"
    )
    
    args = parser.parse_args()
    
    # Check input files
    if not os.path.exists(args.pca_scores_csv):
        print(f"Error: PCA scores file '{args.pca_scores_csv}' not found")
        sys.exit(1)
    
    if not os.path.exists(args.npz_file):
        print(f"Error: NPZ file '{args.npz_file}' not found")
        sys.exit(1)
    
    if not os.path.exists(args.metadata_csv):
        print(f"Error: Metadata CSV file '{args.metadata_csv}' not found")
        sys.exit(1)
    
    # Create output directory (clustered_pca subdirectory)
    output_dir = os.path.join(os.path.dirname(args.output_prefix), "clustered_pca")
    os.makedirs(output_dir, exist_ok=True)
    
    # Update output prefix to be in clustered_pca directory
    base_name = os.path.basename(args.output_prefix)
    args.output_prefix = os.path.join(output_dir, base_name)
    
    print("="*60)
    print("Clustered PCA Analysis")
    print("="*60)
    print(f"PCA scores file: {args.pca_scores_csv}")
    print(f"NPZ file: {args.npz_file}")
    print(f"Metadata file: {args.metadata_csv}")
    print(f"Output directory: {output_dir}")
    print(f"Output prefix: {args.output_prefix}")
    print("="*60)
    
    # 1) Load existing PCA scores
    print("\n1. Loading existing PCA scores...")
    pca_scores_df = pd.read_csv(args.pca_scores_csv)
    
    # Determine ID column - use clean_id as primary (it's the first column in the CSV)
    # This ensures we use the exact same sequences as in the CSV
    if "clean_id" in pca_scores_df.columns:
        id_col = "clean_id"
    elif "original_sequence_id" in pca_scores_df.columns:
        id_col = "original_sequence_id"
    else:
        id_col = pca_scores_df.columns[0]
    
    print(f"Using ID column: {id_col}")
    print(f"Total sequences in CSV: {len(pca_scores_df)}")
    
    # Remove humans
    print("\n2. Removing humans...")
    human_mask = pca_scores_df[id_col].apply(is_human)
    n_humans = human_mask.sum()
    print(f"  Found {n_humans} human sequences")
    
    pca_scores_df = pca_scores_df[~human_mask].copy()
    print(f"  Remaining sequences (after removing humans): {len(pca_scores_df)}")
    
    # Store the exact sequence IDs we'll use (this is our master list)
    exact_sequence_ids = set(pca_scores_df[id_col].values)
    print(f"  Exact sequence IDs to use: {len(exact_sequence_ids)}")
    
    # Also store mapping from id_col to other possible ID formats for matching
    id_mapping = {}
    if "clean_id" in pca_scores_df.columns and id_col != "clean_id":
        for idx, row in pca_scores_df.iterrows():
            main_id = row[id_col]
            clean_id = row.get("clean_id", "")
            if pd.notna(clean_id):
                id_mapping[str(clean_id)] = main_id
    if "original_sequence_id" in pca_scores_df.columns and id_col != "original_sequence_id":
        for idx, row in pca_scores_df.iterrows():
            main_id = row[id_col]
            orig_id = row.get("original_sequence_id", "")
            if pd.notna(orig_id):
                id_mapping[str(orig_id)] = main_id
    
    # Check for PC1 and PC2
    if "PC1" not in pca_scores_df.columns or "PC2" not in pca_scores_df.columns:
        print("Error: PC1 and PC2 not found in PCA scores file")
        sys.exit(1)
    
    # 3) Assign clusters
    print("\n3. Assigning sequences to 3 clusters based on PC1/PC2 thresholds...")
    pc1 = pca_scores_df["PC1"].values
    pc2 = pca_scores_df["PC2"].values
    
    clusters = assign_clusters(pc1, pc2)
    pca_scores_df["cluster"] = clusters
    
    print("\n  Cluster definitions:")
    print(f"    Group 1: PC1 < 5.0 AND PC2 < 10.0 (left of PC1=5, below PC2=10)")
    print(f"    Group 2: PC1 >= 5.0 AND PC2 < 20.0 (right of PC1=5, below PC2=20)")
    print(f"    Group 3: PC2 >= 20.0 AND PC1 < -10.0 (above PC2=20, below PC1=-10)")
    print("\n  Cluster sizes:")
    for cluster_id in range(3):
        n_in_cluster = (clusters == cluster_id).sum()
        print(f"    Group {cluster_id + 1}: {n_in_cluster} sequences")
    
    # Check for unassigned sequences
    n_unassigned = (clusters == -1).sum()
    if n_unassigned > 0:
        print(f"    Unassigned: {n_unassigned} sequences (do not match any group criteria)")
    
    # 4) Load coordinates and metadata
    print("\n4. Loading coordinate data...")
    coords_data = np.load(args.npz_file)
    coords_dict = {key: coords_data[key] for key in coords_data.files}
    print(f"  Loaded {len(coords_dict)} coordinate sets")
    
    metadata_df = pd.read_csv(args.metadata_csv)
    print(f"  Loaded metadata for {len(metadata_df)} structures")
    
    # Determine metadata ID column
    if "clean_id" in metadata_df.columns:
        metadata_id_col = "clean_id"
    elif "original_sequence_id" in metadata_df.columns:
        metadata_id_col = "original_sequence_id"
    else:
        metadata_id_col = metadata_df.columns[0]
    
    # 5) Process each cluster
    print("\n" + "="*60)
    print("Processing Each Cluster")
    print("="*60)
    
    cluster_results = {}
    
    for cluster_id in range(3):
        print(f"\n{'='*60}")
        print(f"Group {cluster_id + 1} (Cluster {cluster_id})")
        print(f"{'='*60}")
        
        # Get sequences in this cluster
        cluster_mask = clusters == cluster_id
        cluster_sequences = pca_scores_df[cluster_mask][id_col].values
        
        if len(cluster_sequences) < args.min_cluster_size:
            print(f"  Skipping Group {cluster_id + 1}: too few sequences ({len(cluster_sequences)} < {args.min_cluster_size})")
            continue
        
        print(f"  Group {cluster_id + 1} (Cluster {cluster_id}):")
        
        print(f"  Sequences in cluster: {len(cluster_sequences)}")
        
        # Check if single-copy clade
        single_copy_info = check_single_copy_clade(cluster_sequences)
        print(f"  Single-copy check:")
        print(f"    Number of genera: {single_copy_info['n_genera']}")
        print(f"    Number of sequences: {single_copy_info['n_sequences']}")
        print(f"    Genus diversity: {single_copy_info['genus_diversity']:.3f}")
        print(f"    Likely single-copy: {single_copy_info['is_single_copy']}")
        if single_copy_info['unique_genera']:
            print(f"    Genera: {', '.join(list(single_copy_info['unique_genera'])[:10])}")
        
        # Match cluster sequences to coordinates
        # IMPORTANT: Only use sequences that are in the exact_sequence_ids set
        cluster_coords = {}
        matched_sequences = []
        
        # Create mapping: clean_id -> coordinates from metadata
        clean_id_to_coords = {}
        if "clean_id" in metadata_df.columns:
            for idx, row in metadata_df.iterrows():
                clean_id = str(row["clean_id"])
                if clean_id in coords_dict:
                    clean_id_to_coords[clean_id] = coords_dict[clean_id]
        
        # Also map original_sequence_id if available
        original_id_to_coords = {}
        if "original_sequence_id" in metadata_df.columns:
            for idx, row in metadata_df.iterrows():
                orig_id = str(row["original_sequence_id"])
                clean_id = str(row.get("clean_id", ""))
                # Try to find coordinates for this sequence
                if clean_id in coords_dict:
                    original_id_to_coords[orig_id] = coords_dict[clean_id]
                elif orig_id in coords_dict:
                    original_id_to_coords[orig_id] = coords_dict[orig_id]
        
        # Also try direct matching with coordinate keys
        for coord_key in coords_dict.keys():
            coord_key_str = str(coord_key)
            if coord_key_str not in clean_id_to_coords:
                clean_id_to_coords[coord_key_str] = coords_dict[coord_key]
        
        for seq_id in cluster_sequences:
            # Only process if this sequence is in our exact list (or maps to one)
            seq_str = str(seq_id)
            main_seq_id = seq_id
            if seq_id not in exact_sequence_ids:
                # Check if this ID maps to one in our exact list
                if seq_str in id_mapping:
                    main_seq_id = id_mapping[seq_str]
                else:
                    continue
                
            seq_str = str(main_seq_id)
            matched = False
            
            # Try multiple matching strategies to find coordinates
            # 1. Direct match with coordinate keys (using clean_id)
            if seq_str in coords_dict:
                cluster_coords[main_seq_id] = coords_dict[seq_str]
                matched_sequences.append(main_seq_id)
                matched = True
            # 2. Match via clean_id mapping
            elif seq_str in clean_id_to_coords:
                cluster_coords[main_seq_id] = clean_id_to_coords[seq_str]
                matched_sequences.append(main_seq_id)
                matched = True
            # 3. Match via original_sequence_id mapping
            elif seq_str in original_id_to_coords:
                cluster_coords[main_seq_id] = original_id_to_coords[seq_str]
                matched_sequences.append(main_seq_id)
                matched = True
            # 4. Try matching via metadata
            else:
                # Try matching by clean_id in metadata
                if "clean_id" in metadata_df.columns:
                    metadata_match = metadata_df[metadata_df["clean_id"] == main_seq_id]
                    if len(metadata_match) > 0:
                        clean_id = str(metadata_match["clean_id"].iloc[0])
                        if clean_id in coords_dict:
                            cluster_coords[main_seq_id] = coords_dict[clean_id]
                            matched_sequences.append(main_seq_id)
                            matched = True
                
                # Try matching by original_sequence_id in metadata
                if not matched and "original_sequence_id" in metadata_df.columns:
                    metadata_match = metadata_df[metadata_df["original_sequence_id"] == main_seq_id]
                    if len(metadata_match) > 0:
                        # Get clean_id from metadata to find coordinates
                        if "clean_id" in metadata_match.columns:
                            clean_id = str(metadata_match["clean_id"].iloc[0])
                            if clean_id in coords_dict:
                                cluster_coords[main_seq_id] = coords_dict[clean_id]
                                matched_sequences.append(main_seq_id)
                                matched = True
                
                # Try fuzzy matching as last resort
                if not matched:
                    for coord_key in coords_dict.keys():
                        coord_key_str = str(coord_key)
                        # Check if seq_id is contained in coord_key or vice versa
                        if seq_str in coord_key_str or coord_key_str in seq_str:
                            cluster_coords[main_seq_id] = coords_dict[coord_key]
                            matched_sequences.append(main_seq_id)
                            matched = True
                            break
        
        if len(cluster_coords) < args.min_cluster_size:
            print(f"  Skipping Group {cluster_id + 1}: too few matched coordinates ({len(cluster_coords)} < {args.min_cluster_size})")
            continue
        
        print(f"  Matched coordinates: {len(cluster_coords)}")
        
        # Re-run PCA for this cluster
        print(f"  Re-running structural PCA for Group {cluster_id + 1}...")
        
        try:
            # Align structures
            aligned_coords, reference_key = align_to_reference(cluster_coords)
            
            # Prepare coordinate matrix
            X, structure_ids, actual_lengths = prepare_coordinate_matrix(aligned_coords)
            
            # Standardize
            X_scaled, scaler = standardize_coordinates(X)
            
            # Perform PCA
            pca, X_pca = perform_pca(X_scaled, n_components=None)
            
            # Save results
            cluster_output_prefix = f"{args.output_prefix}_cluster_{cluster_id}"
            scores_file, variance_file, loadings_file = save_results(
                X_pca, structure_ids, pca, metadata_df, cluster_output_prefix
            )
            
            # Create PC1 vs PC2 scatter plot for this cluster
            print(f"  Creating PC1 vs PC2 plot for Group {cluster_id + 1}...")
            plot_created = False
            if HAS_MATPLOTLIB:
                try:
                    # Only create PC1 vs PC2 scatter plot
                    if X_pca.shape[1] >= 2:
                        fig = plt.figure(figsize=(10, 8))
                        
                        # Classify species groups for coloring (same as phylomorphospace_colored.R)
                        structure_groups = [classify_species_group(seq_id) for seq_id in structure_ids]
                        
                        # Color points by species group
                        unique_groups = list(set(structure_groups))
                        group_colors_map = {
                            "Cetacean": "steelblue",
                            "Terrestrial_Artiodactyl": "darkorange",
                            "Other": "gray50"
                        }
                        
                        # Plot each group separately with appropriate colors
                        for group in unique_groups:
                            group_mask = [g == group for g in structure_groups]
                            group_pc1 = X_pca[group_mask, 0]
                            group_pc2 = X_pca[group_mask, 1]
                            color = group_colors_map.get(group, "gray50")
                            # Use alpha=0.7 to match phylomorphospace_colored.R
                            plt.scatter(group_pc1, group_pc2, c=color, alpha=0.7, s=50, label=group)
                        
                        plt.xlabel(f'PC1 ({pca.explained_variance_ratio_[0]*100:.1f}% variance)', fontsize=12)
                        plt.ylabel(f'PC2 ({pca.explained_variance_ratio_[1]*100:.1f}% variance)', fontsize=12)
                        plt.title(f'Group {cluster_id + 1}: PC1 vs PC2 Structural Morphospace', fontsize=14)
                        plt.legend(loc='best', fontsize=9)
                        plt.grid(True, alpha=0.3)
                        plt.tight_layout()
                        plot_file = f"{cluster_output_prefix}_PC1_PC2_scatter.pdf"
                        plt.savefig(plot_file, dpi=300, bbox_inches='tight')
                        plt.close(fig)
                        print(f"  [OK] PC1 vs PC2 plot saved: {plot_file}")
                        plot_created = True
                    else:
                        print(f"  [WARNING] Not enough components for PC1 vs PC2 plot (only {X_pca.shape[1]} components)")
                except Exception as e:
                    print(f"  [ERROR] Failed to create plot for Group {cluster_id + 1}: {e}")
                    print("  Will try to create plot from CSV files at the end...")
                    if HAS_MATPLOTLIB:
                        try:
                            plt.close('all')
                        except:
                            pass
            else:
                print("  Skipping plot (matplotlib not available)")
            
            # Store results for this cluster
            cluster_results[cluster_id] = {
                'n_sequences': len(cluster_sequences),
                'n_matched': len(cluster_coords),
                'single_copy_info': single_copy_info,
                'pca': pca,
                'scores_file': scores_file,
                'variance_file': variance_file,
                'plot_created': plot_created,
                'plot_file': f"{cluster_output_prefix}_PC1_PC2_scatter.pdf"
            }
            
            print(f"  [OK] Group {cluster_id + 1} PCA complete")
            print(f"    Explained variance (PC1): {pca.explained_variance_ratio_[0]:.4f}")
            print(f"    Explained variance (PC2): {pca.explained_variance_ratio_[1]:.4f}")
            
        except Exception as e:
            print(f"  [ERROR] Error processing Group {cluster_id + 1}: {e}")
            import traceback
            traceback.print_exc()
            continue
    
    # 6) Create summary visualization
    print("\n" + "="*60)
    print("Creating Summary Visualization")
    print("="*60)
    
    if not HAS_MATPLOTLIB:
        print("  Skipping visualization (matplotlib not available)")
    else:
        # Plot original PC1/PC2 with cluster assignments
        fig, axes = plt.subplots(1, 2, figsize=(16, 8))
        
        # Left plot: Original PCA with cluster colors
        ax1 = axes[0]
        colors = ['red', 'blue', 'green']
        for cluster_id in range(3):
            cluster_mask = clusters == cluster_id
            if cluster_mask.sum() > 0:
                ax1.scatter(
                    pc1[cluster_mask],
                    pc2[cluster_mask],
                    c=colors[cluster_id],
                    label=f'Group {cluster_id + 1} (n={cluster_mask.sum()})',
                    alpha=0.6,
                    s=50
                )
        
        # Plot unassigned points in gray
        unassigned_mask = clusters == -1
        if unassigned_mask.sum() > 0:
            ax1.scatter(
                pc1[unassigned_mask],
                pc2[unassigned_mask],
                c='gray',
                label=f'Unassigned (n={unassigned_mask.sum()})',
                alpha=0.3,
                s=30
            )
        
        # Add threshold lines
        ax1.axvline(5.0, color='black', linestyle='--', alpha=0.5, label='PC1 = 5')
        ax1.axvline(-10.0, color='black', linestyle=':', alpha=0.5, label='PC1 = -10')
        ax1.axhline(10.0, color='black', linestyle='--', alpha=0.5, label='PC2 = 10')
        ax1.axhline(20.0, color='black', linestyle='--', alpha=0.5, label='PC2 = 20')
        
        ax1.set_xlabel('PC1')
        ax1.set_ylabel('PC2')
        ax1.set_title('Original PCA with 3-Group Assignment')
        ax1.legend()
        ax1.grid(True, alpha=0.3)
        
        # Right plot: Group sizes
        ax2 = axes[1]
        group_sizes = [(clusters == i).sum() for i in range(3)]
        group_labels = [f'Group {i+1}\n(n={size})' for i, size in enumerate(group_sizes)]
        ax2.bar(range(3), group_sizes, color=colors, alpha=0.7)
        ax2.set_xlabel('Group')
        ax2.set_ylabel('Number of Sequences')
        ax2.set_title('Group Sizes')
        ax2.set_xticks(range(3))
        ax2.set_xticklabels([f'Group {i+1}' for i in range(3)])
        ax2.grid(True, alpha=0.3, axis='y')
        
        plt.tight_layout()
        summary_plot_file = f"{args.output_prefix}_cluster_summary.pdf"
        plt.savefig(summary_plot_file, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"Saved: {summary_plot_file}")
    
    # 7) Try to create plots from CSV files if they weren't created during analysis
    print("\n" + "="*60)
    print("Creating Plots from CSV Files (if needed)")
    print("="*60)
    
    plots_created_from_csv = 0
    for cluster_id, results in cluster_results.items():
        if not results.get('plot_created', False):
            scores_file = results.get('scores_file')
            plot_file = results.get('plot_file')
            if scores_file and plot_file and os.path.exists(scores_file):
                print(f"\nCreating plot for Group {cluster_id + 1} from CSV file...")
                if create_plot_from_csv(scores_file, plot_file, cluster_id):
                    print(f"  [OK] Plot created: {plot_file}")
                    plots_created_from_csv += 1
                    results['plot_created'] = True
                else:
                    print(f"  [WARNING] Failed to create plot for Group {cluster_id + 1}")
    
    if plots_created_from_csv > 0:
        print(f"\nCreated {plots_created_from_csv} plots from CSV files")
    
    # 8) Save cluster assignments
    print("\n" + "="*60)
    print("Saving Results")
    print("="*60)
    
    # Save cluster assignments
    cluster_assignments_file = f"{args.output_prefix}_cluster_assignments.csv"
    pca_scores_df.to_csv(cluster_assignments_file, index=False)
    print(f"Saved: {cluster_assignments_file}")
    
    # Save summary
    summary_data = []
    for cluster_id, results in cluster_results.items():
        summary_data.append({
            'group': cluster_id + 1,
            'cluster_id': cluster_id,
            'n_sequences': results['n_sequences'],
            'n_matched_coordinates': results['n_matched'],
            'is_single_copy': results['single_copy_info']['is_single_copy'],
            'n_genera': results['single_copy_info']['n_genera'],
            'genus_diversity': results['single_copy_info']['genus_diversity'],
            'pc1_variance_explained': results['pca'].explained_variance_ratio_[0] if len(results['pca'].explained_variance_ratio_) > 0 else 0,
            'pc2_variance_explained': results['pca'].explained_variance_ratio_[1] if len(results['pca'].explained_variance_ratio_) > 1 else 0
        })
    
    summary_df = pd.DataFrame(summary_data)
    summary_file = f"{args.output_prefix}_cluster_summary.csv"
    summary_df.to_csv(summary_file, index=False)
    print(f"Saved: {summary_file}")
    
    print("\n" + "="*60)
    print("Analysis Complete!")
    print("="*60)
    print(f"\nSummary:")
    print(f"  Total sequences (excluding humans): {len(pca_scores_df)}")
    print(f"  Groups processed: {len(cluster_results)}")
    for cluster_id, results in cluster_results.items():
        print(f"    Group {cluster_id + 1}: {results['n_sequences']} sequences, "
              f"single-copy: {results['single_copy_info']['is_single_copy']}")
    
    n_unassigned = (clusters == -1).sum()
    if n_unassigned > 0:
        print(f"    Unassigned: {n_unassigned} sequences")


if __name__ == "__main__":
    main()

