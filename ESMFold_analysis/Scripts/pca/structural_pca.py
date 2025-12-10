#!/usr/bin/env python3
"""
Structural PCA Analysis Script

Performs structural alignment and PCA on protein coordinate data:
1. Loads CA coordinates from NPZ file
2. Structural alignment (superposition) using Kabsch algorithm
3. Prepare coordinate data matrix
4. Standardization (z-score normalization)
5. Perform PCA

Usage:
    python structural_pca.py <npz_file> <metadata_csv> [options]
"""

import argparse
import os
import sys
import numpy as np
import pandas as pd
from scipy.spatial.transform import Rotation
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

# Try to import matplotlib, handle gracefully if PIL is missing
try:
    import matplotlib
    matplotlib.use('Agg')  # Use non-interactive backend
    import matplotlib.pyplot as plt
    import seaborn as sns
    HAS_MATPLOTLIB = True
except (ImportError, ModuleNotFoundError) as e:
    HAS_MATPLOTLIB = False
    plt = None
    sns = None
    # Create dummy objects to avoid NameError
    class DummyPlot:
        def __getattr__(self, name):
            def dummy(*args, **kwargs):
                pass
            return dummy
    plt = DummyPlot()
    sns = DummyPlot()


def kabsch_alignment(P, Q):
    """
    Align structure P to structure Q using Kabsch algorithm.
    
    Parameters:
    -----------
    P : ndarray, shape (n_atoms, 3)
        Source coordinates to be aligned
    Q : ndarray, shape (n_atoms, 3)
        Target coordinates (reference)
    
    Returns:
    --------
    P_aligned : ndarray, shape (n_atoms, 3)
        Aligned coordinates of P
    R : ndarray, shape (3, 3)
        Rotation matrix
    """
    # Center both structures
    P_centered = P - P.mean(axis=0)
    Q_centered = Q - Q.mean(axis=0)
    
    # Compute covariance matrix
    H = P_centered.T @ Q_centered
    
    # SVD
    U, S, Vt = np.linalg.svd(H)
    
    # Compute rotation matrix
    R = Vt.T @ U.T
    
    # Ensure proper rotation (det(R) = 1)
    if np.linalg.det(R) < 0:
        Vt[-1, :] *= -1
        R = Vt.T @ U.T
    
    # Apply rotation and translation
    P_aligned = (P_centered @ R.T) + Q.mean(axis=0)
    
    return P_aligned, R


def align_to_reference(coords_dict, reference_key=None):
    """
    Align all structures to a reference structure.
    
    Parameters:
    -----------
    coords_dict : dict
        Dictionary of {key: coords_array} where coords_array is (n_residues, 3)
    reference_key : str, optional
        Key of reference structure. If None, uses the longest structure.
    
    Returns:
    --------
    aligned_coords : dict
        Dictionary of aligned coordinates
    reference_key : str
        Key of the reference structure used
    """
    print("\n" + "="*60)
    print("Step 3: Structural Alignment (Superposition)")
    print("="*60)
    
    # Find reference structure (longest if not specified)
    if reference_key is None:
        reference_key = max(coords_dict.keys(), key=lambda k: len(coords_dict[k]))
        print(f"Using longest structure as reference: {reference_key}")
        print(f"  Reference length: {len(coords_dict[reference_key])} residues")
    else:
        if reference_key not in coords_dict:
            raise ValueError(f"Reference key '{reference_key}' not found in coordinates")
        print(f"Using specified reference: {reference_key}")
    
    reference = coords_dict[reference_key]
    aligned_coords = {reference_key: reference.copy()}  # Reference stays as-is
    
    # Align each structure to reference
    # Note: For structures of different lengths, we need to handle this carefully
    # For now, we'll align structures of the same length, or use a subset
    ref_len = len(reference)
    
    print(f"\nAligning structures to reference (length: {ref_len} residues)...")
    aligned_count = 0
    skipped_count = 0
    
    for key, coords in coords_dict.items():
        if key == reference_key:
            continue
        
        # For structures of different lengths, we'll use the minimum length
        # This is a simplified approach - for production, consider sequence alignment first
        min_len = min(ref_len, len(coords))
        
        if min_len < 10:  # Skip very short structures
            print(f"  Skipping {key}: too short ({len(coords)} residues)")
            skipped_count += 1
            continue
        
        # Align first min_len residues
        P_subset = coords[:min_len]
        Q_subset = reference[:min_len]
        
        try:
            P_aligned, _ = kabsch_alignment(P_subset, Q_subset)
            
            # If original was longer, pad with aligned coordinates
            if len(coords) > min_len:
                # For padding, we can extend the last aligned coordinate
                # or use the original unaligned coordinates for the rest
                # Simple approach: keep original for extra residues
                aligned_full = np.vstack([P_aligned, coords[min_len:]])
            else:
                aligned_full = P_aligned
            
            aligned_coords[key] = aligned_full
            aligned_count += 1
            
            if aligned_count % 50 == 0:
                print(f"  Aligned {aligned_count} structures...")
                
        except Exception as e:
            print(f"  Warning: Failed to align {key}: {e}")
            skipped_count += 1
            continue
    
    print(f"\nAlignment complete:")
    print(f"  Successfully aligned: {aligned_count} structures")
    print(f"  Skipped: {skipped_count} structures")
    print(f"  Reference: {reference_key}")
    
    return aligned_coords, reference_key


def prepare_coordinate_matrix(aligned_coords, min_length=None, max_length=None):
    """
    Prepare coordinate data matrix for PCA.
    Each structure becomes a feature vector of concatenated x, y, z coordinates.
    
    Parameters:
    -----------
    aligned_coords : dict
        Dictionary of aligned coordinates
    min_length : int, optional
        Minimum number of residues (filter shorter structures)
    max_length : int, optional
        Maximum number of residues (pad or truncate to this length)
    
    Returns:
    --------
    X : ndarray, shape (n_structures, 3*N)
        Coordinate matrix where each row is a flattened structure
    structure_ids : list
        List of structure IDs in same order as X
    actual_lengths : dict
        Dictionary mapping structure IDs to actual lengths
    """
    print("\n" + "="*60)
    print("Step 4: Coordinate Data Preparation")
    print("="*60)
    
    # Filter structures by length if specified
    filtered_coords = {}
    if min_length is not None:
        filtered_coords = {k: v for k, v in aligned_coords.items() 
                          if len(v) >= min_length}
        print(f"Filtered to structures with >= {min_length} residues: {len(filtered_coords)} structures")
    else:
        filtered_coords = aligned_coords
    
    # Determine target length
    if max_length is None:
        # Use the most common length, or median length
        lengths = [len(coords) for coords in filtered_coords.values()]
        target_length = int(np.median(lengths))
        print(f"Using median length as target: {target_length} residues")
    else:
        target_length = max_length
        print(f"Using specified target length: {target_length} residues")
    
    # Prepare matrix
    structure_ids = []
    actual_lengths = {}
    coordinate_vectors = []
    
    for key, coords in filtered_coords.items():
        n_residues = len(coords)
        actual_lengths[key] = n_residues
        
        # Pad or truncate to target length
        if n_residues < target_length:
            # Pad with last coordinate (or could use interpolation)
            padding = np.tile(coords[-1:], (target_length - n_residues, 1))
            coords_padded = np.vstack([coords, padding])
        elif n_residues > target_length:
            # Truncate to target length
            coords_padded = coords[:target_length]
        else:
            coords_padded = coords
        
        # Flatten to vector: [x1, y1, z1, x2, y2, z2, ...]
        coords_vector = coords_padded.flatten()
        coordinate_vectors.append(coords_vector)
        structure_ids.append(key)
    
    X = np.array(coordinate_vectors)
    
    print(f"\nCoordinate matrix prepared:")
    print(f"  Number of structures: {X.shape[0]}")
    print(f"  Feature dimension: {X.shape[1]} (3 × {target_length} residues)")
    print(f"  Matrix shape: {X.shape}")
    
    return X, structure_ids, actual_lengths


def standardize_coordinates(X):
    """
    Standardize coordinate data using z-score normalization.
    
    Parameters:
    -----------
    X : ndarray, shape (n_samples, n_features)
        Coordinate matrix
    
    Returns:
    --------
    X_scaled : ndarray, shape (n_samples, n_features)
        Standardized coordinate matrix
    scaler : StandardScaler
        Fitted scaler object
    """
    print("\n" + "="*60)
    print("Step 5: Standardization (Z-score Normalization)")
    print("="*60)
    
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)
    
    print(f"\nStandardization complete:")
    print(f"  Mean after scaling: {X_scaled.mean(axis=0)[:5]}... (should be ~0)")
    print(f"  Std after scaling: {X_scaled.std(axis=0)[:5]}... (should be ~1)")
    print(f"  Overall mean: {X_scaled.mean():.6f}")
    print(f"  Overall std: {X_scaled.std():.6f}")
    
    return X_scaled, scaler


def perform_pca(X_scaled, n_components=None):
    """
    Perform Principal Component Analysis on standardized coordinates.
    
    Parameters:
    -----------
    X_scaled : ndarray, shape (n_samples, n_features)
        Standardized coordinate matrix
    n_components : int, optional
        Number of components to compute. If None, computes all components.
    
    Returns:
    --------
    pca : PCA
        Fitted PCA object
    X_pca : ndarray, shape (n_samples, n_components)
        Principal component scores
    """
    print("\n" + "="*60)
    print("Step 6: Principal Component Analysis")
    print("="*60)
    
    if n_components is None:
        n_components = min(X_scaled.shape[0], X_scaled.shape[1])
    
    pca = PCA(n_components=n_components)
    X_pca = pca.fit_transform(X_scaled)
    
    # Calculate explained variance
    explained_variance = pca.explained_variance_ratio_
    cumulative_variance = np.cumsum(explained_variance)
    
    print(f"\nPCA complete:")
    print(f"  Number of components: {len(explained_variance)}")
    print(f"  Number of samples: {X_pca.shape[0]}")
    print(f"\nExplained variance (top 10 components):")
    for i in range(min(10, len(explained_variance))):
        print(f"  PC{i+1}: {explained_variance[i]:.4f} ({explained_variance[i]*100:.2f}%) "
              f"[Cumulative: {cumulative_variance[i]:.4f} ({cumulative_variance[i]*100:.2f}%)]")
    
    return pca, X_pca


def save_results(X_pca, structure_ids, pca, metadata_df, output_prefix):
    """Save PCA results to files."""
    print("\n" + "="*60)
    print("Saving Results")
    print("="*60)
    
    # Create output directory if needed
    output_dir = os.path.dirname(output_prefix) if os.path.dirname(output_prefix) else "."
    os.makedirs(output_dir, exist_ok=True)
    
    # Save PC scores
    scores_df = pd.DataFrame(
        X_pca,
        index=structure_ids,
        columns=[f"PC{i+1}" for i in range(X_pca.shape[1])]
    )
    
    # Merge with metadata if available
    if metadata_df is not None:
        # Try to merge on clean_id
        scores_df = scores_df.reset_index()
        scores_df.columns = ['clean_id'] + list(scores_df.columns[1:])
        scores_df = scores_df.merge(metadata_df, on='clean_id', how='left')
    
    scores_file = f"{output_prefix}_pca_scores.csv"
    scores_df.to_csv(scores_file, index=False)
    print(f"Saved PC scores: {scores_file}")
    
    # Save explained variance
    variance_df = pd.DataFrame({
        'PC': [f"PC{i+1}" for i in range(len(pca.explained_variance_ratio_))],
        'Explained_Variance_Ratio': pca.explained_variance_ratio_,
        'Cumulative_Variance_Ratio': np.cumsum(pca.explained_variance_ratio_),
        'Explained_Variance': pca.explained_variance_
    })
    variance_file = f"{output_prefix}_explained_variance.csv"
    variance_df.to_csv(variance_file, index=False)
    print(f"Saved explained variance: {variance_file}")
    
    # Save loadings (first 10 PCs)
    n_components_to_save = min(10, X_pca.shape[1])
    loadings_df = pd.DataFrame(
        pca.components_[:n_components_to_save].T,
        columns=[f"PC{i+1}" for i in range(n_components_to_save)]
    )
    loadings_file = f"{output_prefix}_pca_loadings.csv"
    loadings_df.to_csv(loadings_file, index=False)
    print(f"Saved PC loadings: {loadings_file}")
    
    return scores_file, variance_file, loadings_file


def plot_results(X_pca, structure_ids, pca, output_prefix):
    """Create visualization plots."""
    if not HAS_MATPLOTLIB:
        print("\n" + "="*60)
        print("Creating Visualizations")
        print("="*60)
        print("Skipping plots (matplotlib not available)")
        return
    
    print("\n" + "="*60)
    print("Creating Visualizations")
    print("="*60)
    
    # Set style
    if sns is not None:
        sns.set_style("whitegrid")
    plt.rcParams['figure.figsize'] = (10, 8)
    
    # 1. Scree plot
    plt.figure()
    n_components_to_plot = min(20, len(pca.explained_variance_ratio_))
    plt.plot(range(1, n_components_to_plot + 1), 
             pca.explained_variance_ratio_[:n_components_to_plot], 
             'bo-', linewidth=2, markersize=8)
    plt.xlabel('Principal Component', fontsize=12)
    plt.ylabel('Explained Variance Ratio', fontsize=12)
    plt.title('Scree Plot: Explained Variance by Principal Component', fontsize=14)
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(f"{output_prefix}_scree_plot.pdf", dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved scree plot: {output_prefix}_scree_plot.pdf")
    
    # 2. Cumulative variance plot
    plt.figure()
    cumulative_variance = np.cumsum(pca.explained_variance_ratio_)
    n_components_to_plot = min(20, len(cumulative_variance))
    plt.plot(range(1, n_components_to_plot + 1), 
             cumulative_variance[:n_components_to_plot], 
             'ro-', linewidth=2, markersize=8)
    plt.xlabel('Number of Principal Components', fontsize=12)
    plt.ylabel('Cumulative Explained Variance Ratio', fontsize=12)
    plt.title('Cumulative Explained Variance', fontsize=14)
    plt.grid(True, alpha=0.3)
    plt.axhline(y=0.8, color='g', linestyle='--', label='80% variance')
    plt.axhline(y=0.9, color='orange', linestyle='--', label='90% variance')
    plt.legend()
    plt.tight_layout()
    plt.savefig(f"{output_prefix}_cumulative_variance.pdf", dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved cumulative variance plot: {output_prefix}_cumulative_variance.pdf")
    
    # 3. PC1 vs PC2 scatter plot
    if X_pca.shape[1] >= 2:
        plt.figure()
        plt.scatter(X_pca[:, 0], X_pca[:, 1], alpha=0.6, s=50)
        plt.xlabel(f'PC1 ({pca.explained_variance_ratio_[0]*100:.1f}% variance)', fontsize=12)
        plt.ylabel(f'PC2 ({pca.explained_variance_ratio_[1]*100:.1f}% variance)', fontsize=12)
        plt.title('PC1 vs PC2: Structural Morphospace', fontsize=14)
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.savefig(f"{output_prefix}_PC1_PC2_scatter.pdf", dpi=300, bbox_inches='tight')
        plt.close()
        print(f"Saved PC1 vs PC2 scatter: {output_prefix}_PC1_PC2_scatter.pdf")
    
    # 4. PC2 vs PC3 scatter plot
    if X_pca.shape[1] >= 3:
        plt.figure()
        plt.scatter(X_pca[:, 1], X_pca[:, 2], alpha=0.6, s=50, color='orange')
        plt.xlabel(f'PC2 ({pca.explained_variance_ratio_[1]*100:.1f}% variance)', fontsize=12)
        plt.ylabel(f'PC3 ({pca.explained_variance_ratio_[2]*100:.1f}% variance)', fontsize=12)
        plt.title('PC2 vs PC3: Structural Morphospace', fontsize=14)
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.savefig(f"{output_prefix}_PC2_PC3_scatter.pdf", dpi=300, bbox_inches='tight')
        plt.close()
        print(f"Saved PC2 vs PC3 scatter: {output_prefix}_PC2_PC3_scatter.pdf")


def main():
    parser = argparse.ArgumentParser(
        description="Perform structural alignment and PCA on protein coordinates",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    # Basic usage
    python structural_pca.py coordinates.npz metadata.csv
    
    # Specify output prefix
    python structural_pca.py coordinates.npz metadata.csv -o results/pca_analysis
    
    # Filter by length and specify reference
    python structural_pca.py coordinates.npz metadata.csv --min-length 80 --reference-key key123
        """
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
        "-o", "--output-prefix",
        default="pca_results",
        help="Output file prefix (default: pca_results)"
    )
    
    parser.add_argument(
        "--min-length",
        type=int,
        help="Minimum number of residues (filter shorter structures)"
    )
    
    parser.add_argument(
        "--max-length",
        type=int,
        help="Maximum number of residues (pad/truncate to this length)"
    )
    
    parser.add_argument(
        "--reference-key",
        help="Key of reference structure for alignment (default: longest structure)"
    )
    
    parser.add_argument(
        "--n-components",
        type=int,
        help="Number of principal components to compute (default: all)"
    )
    
    parser.add_argument(
        "--no-plots",
        action="store_true",
        help="Skip generating plots"
    )
    
    args = parser.parse_args()
    
    # Check input files
    if not os.path.exists(args.npz_file):
        print(f"Error: NPZ file '{args.npz_file}' not found")
        sys.exit(1)
    
    if not os.path.exists(args.metadata_csv):
        print(f"Error: Metadata CSV file '{args.metadata_csv}' not found")
        sys.exit(1)
    
    print("="*60)
    print("Structural PCA Analysis")
    print("="*60)
    print(f"Input NPZ file: {args.npz_file}")
    print(f"Input metadata: {args.metadata_csv}")
    print(f"Output prefix: {args.output_prefix}")
    
    # Load coordinates
    print("\n" + "="*60)
    print("Loading Coordinates")
    print("="*60)
    coords_data = np.load(args.npz_file)
    coords_dict = {key: coords_data[key] for key in coords_data.files}
    print(f"Loaded {len(coords_dict)} coordinate sets")
    
    # Load metadata
    metadata_df = pd.read_csv(args.metadata_csv)
    print(f"Loaded metadata for {len(metadata_df)} structures")
    
    # Step 3: Structural alignment
    aligned_coords, reference_key = align_to_reference(
        coords_dict, 
        reference_key=args.reference_key
    )
    
    # Step 4: Prepare coordinate matrix
    X, structure_ids, actual_lengths = prepare_coordinate_matrix(
        aligned_coords,
        min_length=args.min_length,
        max_length=args.max_length
    )
    
    # Step 5: Standardize
    X_scaled, scaler = standardize_coordinates(X)
    
    # Step 6: Perform PCA
    pca, X_pca = perform_pca(X_scaled, n_components=args.n_components)
    
    # Save results
    save_results(X_pca, structure_ids, pca, metadata_df, args.output_prefix)
    
    # Create plots
    if not args.no_plots:
        plot_results(X_pca, structure_ids, pca, args.output_prefix)
    
    print("\n" + "="*60)
    print("Analysis Complete!")
    print("="*60)
    print(f"\nOutput files:")
    print(f"  - {args.output_prefix}_pca_scores.csv")
    print(f"  - {args.output_prefix}_explained_variance.csv")
    print(f"  - {args.output_prefix}_pca_loadings.csv")
    if not args.no_plots:
        print(f"  - {args.output_prefix}_scree_plot.pdf")
        print(f"  - {args.output_prefix}_cumulative_variance.pdf")
        print(f"  - {args.output_prefix}_PC1_PC2_scatter.pdf")
        print(f"  - {args.output_prefix}_PC2_PC3_scatter.pdf")


if __name__ == "__main__":
    main()

