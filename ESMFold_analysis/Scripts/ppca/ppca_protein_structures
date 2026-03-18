#!/usr/bin/env python3
"""
Apply Phylogenetic Principal Components Analysis (pPCA) to ESMFold protein structure coordinates.

This script follows the methodology from Polly et al. (2013) to apply pPCA to 3D protein
structure coordinates, accounting for phylogenetic relationships among species.

Pipeline Overview:
==================

1. Load ESMFold coordinates and phylogenetic tree
   - Supports NPZ files with coordinate arrays

2. Correspondence identification (Selection Guide):
   
   DEFAULT METHOD: MSA-based Alignment with Complete Coverage (when --alignment provided):
   - Uses multiple sequence alignment to identify evolutionarily homologous positions.
   - Requires 100%% coverage (all structures must have coordinates at each position).
   - Zero imputation bias - uses only positions where all structures have coordinates.
   - RECOMMENDED for most analyses when a sequence alignment is available.
   
   OPTIONAL METHOD 1: MSA with Imputation (--alignment --allow-imputation):
   - Uses MSA but allows incomplete coverage with mean distance imputation.
   - Automatically optimizes coverage threshold to maximize positions while minimizing imputation.
   - Uses mean distance imputation (Polly et al. 2013) for positions with incomplete coverage.
   - Preserves evolutionary signal and maximizes data retention.
   - Use when you need more positions and can accept some imputation.
   
   OPTIONAL METHOD 2: FoldMason Structural Alignment (--use-structural-alignment-only or no --alignment):
   - Uses FoldMason for high-speed progressive multiple structure alignment (MSTA).
   - Method: Foldseek-based structural alphabet (3Di) -> progressive merge.
   - Fully automated, ultra-fast structural core identification.
   - Use when sequence alignment is difficult or poor (highly divergent sequences).
   
   FUNDAMENTAL CONSTRAINT: PCA requires all structures to have the SAME number of features.
   This means we must analyze a consistent set of homologous positions across all structures.
   
HANDLING MISSING DATA IN MSA-BASED ALIGNMENT (DISTANCE-BASED IMPUTATION):
   - MSA-based alignment (DEFAULT METHOD) requires 100%% coverage (no imputation).
   - OPTIONAL: Use --allow-imputation to enable automatic coverage optimization to maximize positions
     while minimizing imputation via distance-based imputation.
   - CRITICAL: Rather than dropping structures with gaps, this script uses
     MEAN DISTANCE IMPUTATION (Polly et al. 2013) to ensure 100% specimen retention.
   
   Automatic Coverage Optimization (DEFAULT when --alignment provided):
      - Tests multiple coverage thresholds (95%, 90%, 85%, 80%, 75%, 70%, 65%, 60%)
      - For each threshold, calculates number of positions and predicted imputation rate
      - Selects optimal threshold that maximizes positions while minimizing imputation
      - Goal: Find best balance between data retention and imputation quality
   
   Gap Detection:
      - When MSA identifies homologous positions with incomplete coverage (e.g., 85% of structures have a residue),
        structures with gaps at that position are flagged as missing.
      - Missing positions are marked with placeholder coordinates (0,0,0) and a missing flag.
      - Only positions meeting the optimized coverage threshold are included.
   
   Distance Feature Calculation:
      - For each structure, pairwise Cα-Cα distances are calculated between all residue pairs.
      - If EITHER residue in a pair is missing (gap), that distance is marked as -1.0 (placeholder).
      - Only distances where BOTH residues exist are calculated from actual coordinates.
   
   Mean Distance Imputation (Polly et al. 2013 method):
      - For each distance pair (e.g., distance between residue i and residue j):
        a) Find all structures where BOTH residues i and j exist (valid distances).
        b) Calculate the mean distance across these valid structures.
        c) Replace all -1.0 placeholders (from structures missing one or both residues) with this mean.
      - This ensures every structure has a complete distance feature vector for PCA.
   
   Coordinate Imputation (for GPA alignment only):
      - For missing residues, coordinates are imputed as the mean position of that residue
        across all structures that have it.
      - This allows Generalized Procrustes Analysis to align all structures.
      - NOTE: PCA uses the imputed distance matrix, NOT the imputed coordinates.
   
   Imputation Rate Calculation:
      - Imputation Rate = (Number of imputed distance values / Total distance values) × 100%
      - This is printed in the log for quality assessment.
      - Both predicted (during optimization) and actual (after distance calculation) rates are reported.
   
   QUALITY GUIDELINES (when using --allow-imputation):
   - < 5% imputation rate: Excellent, results are highly reliable with negligible bias.
   - 5-7% imputation rate: Good, acceptable for most analyses.
   - 7-10% imputation rate: Moderate, acceptable but consider using default (complete coverage) for stricter control.
   - > 10% imputation rate: High bias risk, strongly recommend using default (complete coverage) or --use-structural-alignment-only.
   
   Enabling Imputation:
   - DEFAULT: Complete coverage required (100% coverage, no imputation).
   - Use --allow-imputation to enable imputation for positions with incomplete coverage.
   - This may result in more positions but introduces imputation bias.
   - Recommended when you need more positions and can accept some imputation.

3. Generalized Procrustes Analysis (GPA):
   - Optimal alignment of structures with consistent features (from Step 2).
   - Removes translation, rotation, and scale (Full Procrustes).
   - Iterative alignment to consensus mean structure.

4. Feature extraction (default: distance-based, RECOMMENDED):
   - Computes all pairwise Cα-Cα distances between aligned residues.
   - Uses the Imputed Distance Matrix for PCA.
   - Captures local geometry and is robust to global alignment issues.

5. Build phylogenetic covariance matrix from tree:
   - Shared branch lengths between all pairs of species (MRCA to root).

6. Perform standard PCA and pPCA:
   - pPCA removes phylogenetic signal during dimensionality reduction.
   - PC scores indicate how each protein varies along non-phylogenetic distance axes.

Key Features:
=============

- **100% Structure Retention**: Uses Mean Distance Imputation to keep every structure in your dataset.
- **High-Speed Alignment**: Uses FoldMason (Foldseek-based) for ultra-fast structural comparison.
- **Hardware Optimized**: Specifically configured for AVX2-enabled cluster nodes.
- **Distance-based PCA**: More appropriate for protein folds than raw coordinate-based PCA.


References:
===========
Polly, P. D., et al. (2013). "Phylogenetic principal components analysis and geometric morphometrics."
Hystrix, 24(1), 33-41.

Revell, L. J. (2009). "Size-correction and principal components for interspecific comparative studies."
Evolution, 63(12), 3258-3268.
"""


import numpy as np
import scipy.linalg
from scipy.spatial.distance import pdist, squareform
from scipy.optimize import minimize, linear_sum_assignment
from scipy.linalg import svd
import pickle
from pathlib import Path
import warnings
import re
import subprocess
import os
warnings.filterwarnings('ignore')

try:
    import pandas as pd
    HAS_PANDAS = True
except ImportError:
    HAS_PANDAS = False
    print("Warning: pandas not available. CSV export will be skipped.")

try:
    import matplotlib
    matplotlib.use('Agg')  # Non-interactive backend
    import matplotlib.pyplot as plt
    import seaborn as sns
    HAS_MATPLOTLIB = True
    sns.set_style("whitegrid")
except ImportError:
    HAS_MATPLOTLIB = False
    print("Warning: matplotlib/seaborn not available. Figure export will be skipped.")

try:
    from Bio import Phylo
    from Bio.Phylo.TreeConstruction import DistanceMatrix
    from Bio.PDB import Superimposer
    from Bio.PDB.vectors import calc_angle, calc_dihedral
    HAS_BIOPYTHON = True
except ImportError:
    print("Warning: BioPython not available. Will use alternative tree parsing.")
    Phylo = None
    HAS_BIOPYTHON = False

try:
    from scipy.spatial.distance import cdist
    HAS_SCIPY = True
except ImportError:
    HAS_SCIPY = False

def procrustes_superimpose(X, Y, scale=True):
    """
    Procrustes superimposition: align Y to X.
    
    Parameters:
    -----------
    X : array, shape (n_points, 3)
        Reference structure
    Y : array, shape (n_points, 3)
        Structure to align
    scale : bool, default=True
        If True, scale structures to unit size (Full Procrustes).
        If False, preserve size differences (Partial Procrustes).
    
    Returns:
    --------
    Y_aligned : array, shape (n_points, 3)
        Aligned version of Y
    R : array, shape (3, 3)
        Rotation matrix
    t : array, shape (3,)
        Translation vector
    s : float
        Scale factor (1.0 if scale=False)
    """
    # Center both structures
    X_centered = X - X.mean(axis=0)
    Y_centered = Y - Y.mean(axis=0)
    
    # Calculate centroid sizes (RMSD from centroid)
    X_size = np.sqrt(np.sum(X_centered ** 2))
    Y_size = np.sqrt(np.sum(Y_centered ** 2))
    
    # Scale to unit size if requested (Full Procrustes)
    if scale:
        if X_size > 1e-10:
            X_scaled = X_centered / X_size
        else:
            X_scaled = X_centered
            X_size = 1.0
        if Y_size > 1e-10:
            Y_scaled = Y_centered / Y_size
            scale_factor = 1.0 / Y_size  # Store original scale for return value
        else:
            Y_scaled = Y_centered
            scale_factor = 1.0
            Y_size = 1.0
    else:
        X_scaled = X_centered
        Y_scaled = Y_centered
        scale_factor = 1.0
    
    # Compute cross-covariance matrix
    H = Y_scaled.T @ X_scaled
    
    # SVD to find optimal rotation
    U, s, Vt = svd(H, full_matrices=False)
    R = U @ Vt
    
    # Ensure proper rotation (det(R) = 1)
    if np.linalg.det(R) < 0:
        U[:, -1] *= -1
        R = U @ Vt
    
    # Apply rotation and translation
    if scale:
        # Both structures are at unit size, so result is also at unit size
        # Then translate to X's centroid
        Y_aligned = (Y_scaled @ R.T) + X.mean(axis=0)
    else:
        # Just rotate and translate (preserve original size)
        Y_aligned = (Y_scaled @ R.T) + X.mean(axis=0)
    
    return Y_aligned, R, X.mean(axis=0) - Y.mean(axis=0) @ R.T, scale_factor

def generalized_procrustes_analysis(structures, scale=True):
    """
    Generalized Procrustes Analysis: align all structures to a common reference.
    
    Parameters:
    -----------
    structures : list of arrays, each shape (n_points, 3)
        List of structures to align
    scale : bool, default=True
        If True, scale structures to unit size (Full Procrustes, removes size differences).
        If False, preserve size differences (Partial Procrustes).
    
    Returns:
    --------
    aligned : list of arrays
        Aligned structures
    mean_structure : array, shape (n_points, 3)
        Mean (consensus) structure
    """
    n_structures = len(structures)
    aligned = [s.copy() for s in structures]
    
    # Initialize with first structure as reference
    mean_structure = structures[0].copy()
    
    # Iterative alignment
    for iteration in range(10):  # Usually converges quickly
        old_mean = mean_structure.copy()
        
        # Align each structure to current mean
        for i in range(n_structures):
            aligned[i], _, _, _ = procrustes_superimpose(mean_structure, aligned[i], scale=scale)
        
        # Update mean
        mean_structure = np.mean(aligned, axis=0)
        
        # Check convergence
        if np.linalg.norm(mean_structure - old_mean) < 1e-6:
            break
    
    return aligned, mean_structure

def parse_sequence_alignment(alignment_file):
    """
    Parse a sequence alignment file (FASTA format with gaps).
    
    Parameters:
    -----------
    alignment_file : str or Path
        Path to alignment file (.aln format, FASTA-like)
    
    Returns:
    --------
    alignment_dict : dict
        Dictionary mapping sequence IDs to aligned sequences (as strings with gaps)
    sequence_lengths : dict
        Dictionary mapping sequence IDs to ungapped sequence lengths
    """
    alignment_dict = {}
    sequence_lengths = {}
    current_seq_id = None
    current_seq = []
    
    with open(alignment_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                # Save previous sequence if exists
                if current_seq_id is not None:
                    aligned_seq = ''.join(current_seq)
                    alignment_dict[current_seq_id] = aligned_seq
                    # Calculate ungapped length
                    sequence_lengths[current_seq_id] = len(aligned_seq.replace('-', ''))
                
                # Start new sequence
                current_seq_id = line[1:].strip()  # Remove '>'
                current_seq = []
            elif line and current_seq_id is not None:
                # Append sequence line (may contain gaps)
                current_seq.append(line)
        
        # Don't forget the last sequence
        if current_seq_id is not None:
            aligned_seq = ''.join(current_seq)
            alignment_dict[current_seq_id] = aligned_seq
            sequence_lengths[current_seq_id] = len(aligned_seq.replace('-', ''))
    
    return alignment_dict, sequence_lengths

def normalize_seq_id_for_matching(seq_id):
    """
    Normalize sequence ID for matching between alignment and structure files.
    
    This handles variations in ID formatting (spaces, dots, underscores, pipes, etc.)
    to enable matching between alignment sequence IDs and structure file keys.
    
    Parameters:
    -----------
    seq_id : str
        Original sequence ID
    
    Returns:
    --------
    normalized : str
        Normalized ID for matching
    """
    import re
    # Remove leading/trailing whitespace
    normalized = str(seq_id).strip()
    
    # Remove pipes (alignment files use | as separators, NPZ files don't)
    normalized = normalized.replace('|', '')
    
    # Replace spaces with underscores
    normalized = normalized.replace(' ', '_')
    
    # Handle XP accession numbers: XP_042112342.1 -> XP_0421123421 (remove dot, concatenate version)
    # Pattern: XP_ followed by digits, then .version_number
    normalized = re.sub(r'(XP_\d+)\.(\d+)', r'\1\2', normalized)
    
    # Handle other accession patterns: remove dots before version numbers
    # Pattern: underscore/dot followed by single digit version at end of accession
    normalized = re.sub(r'(\d+)\.(\d+)(?=[^0-9]|$)', r'\1\2', normalized)
    
    # Handle PREDICTED: -> PREDICTED_ (alignment has colon, NPZ has underscore)
    normalized = normalized.replace('PREDICTED:', 'PREDICTED_')
    normalized = re.sub(r'PREDICTED[_:]+', 'PREDICTED_', normalized)
    
    # Normalize other common prefixes
    normalized = re.sub(r'LOW_QUALITY_PROTEIN[_:]+', 'LOW_QUALITY_PROTEIN_', normalized)
    normalized = re.sub(r'TPA[_:]+', 'TPA_', normalized)
    
    # Remove commas and colons (except those already handled)
    normalized = normalized.replace(',', '').replace(':', '_')
    
    # Remove any remaining dots (convert to nothing, not underscore, to match NPZ format)
    normalized = normalized.replace('.', '')
    
    # Clean up multiple underscores
    normalized = re.sub(r'_{2,}', '_', normalized)
    
    # Remove trailing underscores
    normalized = normalized.rstrip('_')
    
    return normalized

def optimize_msa_coverage_threshold(alignment_dict, structure_names, structures, msa_threshold=0.8, 
                                     candidate_columns=None, msa_to_structure_maps_prelim=None, 
                                     matched_sequences=None):
    """
    Optimize coverage threshold to find best balance between number of positions and imputation rate.
    
    Tests multiple coverage thresholds and selects the one that minimizes imputation while maximizing positions.
    
    Returns:
    --------
    optimal_threshold : float
        Best coverage threshold found
    optimal_positions : list
        MSA column indices at optimal threshold
    predicted_imputation_rate : float
        Predicted imputation rate at optimal threshold
    """
    if candidate_columns is None or msa_to_structure_maps_prelim is None or matched_sequences is None:
        # Need to compute these first - this should not happen in normal flow
        return None, None, None
    
    n_matched = len(matched_sequences)
    n_columns = len(list(matched_sequences.values())[0])
    
    # Test coverage thresholds from high to low
    test_thresholds = [0.95, 0.90, 0.85, 0.80, 0.75, 0.70, 0.65, 0.60]
    results = []
    
    print(f"\n   STEP 3b: Optimizing Coverage Threshold (Finding Best Imputation/Position Trade-off)")
    print(f"   {'='*70}")
    print(f"   Testing multiple coverage thresholds to find optimal balance...")
    print(f"   Goal: Minimize imputation rate while maximizing number of positions")
    
    for threshold in test_thresholds:
        min_required = int(np.ceil(threshold * n_matched))
        
        # Find positions meeting this threshold
        positions_at_threshold = []
        for col_idx in candidate_columns:
            n_with_coords = 0
            for struct_idx in matched_sequences.keys():
                msa_to_struct = msa_to_structure_maps_prelim[struct_idx]
                if col_idx in msa_to_struct:
                    struct_res_idx = msa_to_struct[col_idx]
                    if struct_res_idx < structures[struct_idx].shape[0]:
                        n_with_coords += 1
            
            if n_with_coords >= min_required:
                positions_at_threshold.append(col_idx)
        
        if len(positions_at_threshold) == 0:
            continue
        
        # Calculate predicted imputation rate
        # For each position, count missing residues across all structures
        total_residue_slots = n_matched * len(positions_at_threshold)
        missing_residues = 0
        
        for col_idx in positions_at_threshold:
            for struct_idx in matched_sequences.keys():
                msa_to_struct = msa_to_structure_maps_prelim[struct_idx]
                if col_idx not in msa_to_struct:
                    missing_residues += 1
                else:
                    struct_res_idx = msa_to_struct[col_idx]
                    if struct_res_idx >= structures[struct_idx].shape[0]:
                        missing_residues += 1
        
        # Imputation rate: missing distance pairs / total distance pairs
        # For distance pairs: if position i or j is missing, that pair is imputed
        # Rough estimate: ~2 * missing_residue_rate (since each missing residue affects multiple pairs)
        missing_rate = (missing_residues / total_residue_slots) * 100 if total_residue_slots > 0 else 100
        # More accurate: for n positions, we have n*(n-1)/2 pairs
        # If p% of residues are missing, roughly 2p% of pairs are affected (simplified)
        predicted_imputation_rate = min(100, 2 * missing_rate)  # Conservative estimate
        
        results.append({
            'threshold': threshold,
            'n_positions': len(positions_at_threshold),
            'predicted_imputation_rate': predicted_imputation_rate,
            'missing_rate': missing_rate,
            'positions': positions_at_threshold
        })
        
        print(f"   - Coverage {threshold*100:.0f}%: {len(positions_at_threshold)} positions, "
              f"~{predicted_imputation_rate:.1f}% predicted imputation")
    
    if len(results) == 0:
        return None, None, None
    
    # Find optimal threshold using scoring function
    # Score = (n_positions / max_positions) - (imputation_rate / 10)
    # This balances maximizing positions while minimizing imputation
    max_positions = max(r['n_positions'] for r in results)
    
    best_score = -np.inf
    optimal_result = None
    
    for r in results:
        # Normalize: positions (higher is better), imputation (lower is better)
        position_score = r['n_positions'] / max_positions if max_positions > 0 else 0
        imputation_penalty = r['predicted_imputation_rate'] / 10.0  # Penalty increases with imputation
        
        # Combined score: favor more positions, but penalize high imputation
        # Weight imputation more heavily (imputation >10% is bad)
        score = position_score - (imputation_penalty * 0.5)
        
        # Prefer thresholds with imputation <10%
        if r['predicted_imputation_rate'] <= 10.0:
            score += 0.2  # Bonus for acceptable imputation
        
        if score > best_score:
            best_score = score
            optimal_result = r
    
    if optimal_result is None:
        optimal_result = results[0]  # Fallback to first result
    
    print(f"\n   OPTIMAL THRESHOLD SELECTED:")
    print(f"   - Coverage: {optimal_result['threshold']*100:.0f}%")
    print(f"   - Positions: {optimal_result['n_positions']}")
    print(f"   - Predicted imputation rate: ~{optimal_result['predicted_imputation_rate']:.1f}%")
    
    return optimal_result['threshold'], optimal_result['positions'], optimal_result['predicted_imputation_rate']


def extract_homologous_positions_from_msa(alignment_dict, structure_names, structures, msa_threshold=0.8, min_coverage=0.7, allow_imputation=False, optimize_coverage=True):
    """
    Extract homologous positions from structures based on multiple sequence alignment.
    
    This function:
    1. Identifies columns in the MSA where sequences have non-gap residues
    2. Maps MSA positions to structure residue indices (accounting for gaps)
    3. Extracts corresponding residues from each structure
    
    Parameters:
    -----------
    alignment_dict : dict
        Dictionary mapping sequence IDs to aligned sequences (with gaps)
    structure_names : list
        List of structure names/IDs (keys from NPZ file)
    structures : list
        List of structure coordinate arrays
    msa_threshold : float, default=0.8
        Minimum proportion of sequences that must have non-gap residues at a column
        to be considered as candidate. 0.8 = 80% (default, good for most proteins),
        0.5-0.6 for highly variable domains (e.g., LRR).
    min_coverage : float, default=0.7
        Minimum proportion of structures that must have residues at a column to be
        included in final alignment. 0.7 = 70% (allows some missing structures),
        1.0 = all structures must have residues (strictest).
    
    Returns:
    --------
    homologous_structures : list of arrays
        Structures with only homologous positions (same shape)
    msa_to_structure_maps : list of arrays
        Maps MSA column indices to structure residue indices for each structure
    common_columns : array
        MSA column indices that meet the threshold requirement
    """
    # Normalize structure names for matching
    normalized_structure_names = {normalize_seq_id_for_matching(name): i 
                                  for i, name in enumerate(structure_names)}
    
    # Find sequences in alignment that match structures
    matched_sequences = {}
    matched_indices = {}
    
    for align_seq_id, aligned_seq in alignment_dict.items():
        normalized_align_id = normalize_seq_id_for_matching(align_seq_id)
        
        # Try to match with structure names
        if normalized_align_id in normalized_structure_names:
            struct_idx = normalized_structure_names[normalized_align_id]
            matched_sequences[struct_idx] = aligned_seq
            matched_indices[struct_idx] = align_seq_id
        else:
            # Try partial matching (in case alignment ID is subset of structure ID or vice versa)
            for norm_struct_id, struct_idx in normalized_structure_names.items():
                if normalized_align_id in norm_struct_id or norm_struct_id in normalized_align_id:
                    matched_sequences[struct_idx] = aligned_seq
                    matched_indices[struct_idx] = align_seq_id
                    break
    
    if len(matched_sequences) == 0:
        raise ValueError("No sequences in alignment matched structure names. "
                        "Check sequence ID formatting.")
    
    if len(matched_sequences) < len(structures):
        print(f"   Warning: Only {len(matched_sequences)}/{len(structures)} structures matched to alignment")
        print(f"   Matched: {list(matched_indices.values())[:5]}...")
    
    # Find MSA columns where sequences meet the threshold requirements
    # Strategy: Use msa_threshold to find candidate columns, then apply min_coverage
    # to allow some structures to have gaps (handles missing data)
    n_columns = len(list(matched_sequences.values())[0])
    n_matched = len(matched_sequences)
    min_required_candidate = int(np.ceil(msa_threshold * n_matched))
    min_required_final = int(np.ceil(min_coverage * n_matched))
    
    # First pass: find candidate columns meeting msa_threshold
    candidate_columns = []
    for col_idx in range(n_columns):
        n_with_residue = sum(
            1 for struct_idx in matched_sequences.keys()
            if matched_sequences[struct_idx][col_idx] != '-'
        )
        if n_with_residue >= min_required_candidate:
            candidate_columns.append(col_idx)
    
    # Build MSA-to-structure mappings first to check actual structural coverage
    # Structures may have fewer residues than alignment expects
    msa_to_structure_maps_prelim = {}
    for struct_idx in sorted(matched_sequences.keys()):
        aligned_seq = matched_sequences[struct_idx]
        msa_to_struct = {}
        struct_res_idx = 0
        for msa_col_idx in range(len(aligned_seq)):
            if aligned_seq[msa_col_idx] != '-':
                msa_to_struct[msa_col_idx] = struct_res_idx
                struct_res_idx += 1
        msa_to_structure_maps_prelim[struct_idx] = msa_to_struct
    
    # Second pass: from candidates, keep columns where structures ACTUALLY have coordinates
    # Only count positions where structure has residue
    common_columns = []
    for col_idx in candidate_columns:
        n_with_actual_residue = 0
        for struct_idx in matched_sequences.keys():
            # Check if alignment has residue AND structure has coordinate at that position
            if matched_sequences[struct_idx][col_idx] != '-':
                msa_to_struct = msa_to_structure_maps_prelim[struct_idx]
                if col_idx in msa_to_struct:
                    struct_res_idx = msa_to_struct[col_idx]
                    # Critical: Check if structure actually has this residue
                    if struct_res_idx < structures[struct_idx].shape[0]:
                        n_with_actual_residue += 1
        
        # Only include if enough structures have actual coordinates (not just alignment residues)
        if n_with_actual_residue >= min_required_final:
            common_columns.append(col_idx)
    
    if len(common_columns) == 0:
        raise ValueError(
            f"No columns found meeting coverage requirement (>={min_coverage*100:.0f}% structures with residues, "
            f"from {len(candidate_columns)} candidates at >={msa_threshold*100:.0f}% threshold). "
            f"Try lowering --msa-threshold (e.g., 0.5-0.6) and/or --min-coverage (e.g., 0.5) "
            f"for highly variable domains, or the alignment may be too fragmented for this approach."
        )
    
    # Calculate statistics for warnings
    coverage_pct = (len(common_columns) / n_columns) * 100
    candidate_pct = (len(candidate_columns) / n_columns) * 100 if len(candidate_columns) > 0 else 0
    gap_percentage = 100 - coverage_pct
    
    print(f"\n   STEP 3: MSA Column Filtering (Correspondence Identification)")
    print(f"   {'='*70}")
    print(f"   Total MSA columns: {n_columns}")
    print(f"   Sequences matched: {len(matched_sequences)}/{len(structures)}")
    print(f"\n   3a. Candidate columns (>={msa_threshold*100:.0f}% threshold):")
    print(f"       - Found: {len(candidate_columns)} columns ({candidate_pct:.1f}% of MSA)")
    print(f"       - Dropped: {n_columns - len(candidate_columns)} columns ({100-candidate_pct:.1f}% of MSA)")
    print(f"       - Requirement: >= {min_required_candidate} sequences with residues (out of {n_matched} matched)")
    
    # Store original common_columns count before any optimization/filtering
    n_final_positions_before_complete = len(common_columns)
    
    # Optimize coverage threshold if imputation is enabled and optimization is requested
    if allow_imputation and optimize_coverage:
        optimal_threshold, optimal_positions, predicted_imputation = optimize_msa_coverage_threshold(
            alignment_dict, structure_names, structures, msa_threshold,
            candidate_columns, msa_to_structure_maps_prelim, matched_sequences
        )
        if optimal_positions is not None and len(optimal_positions) > 0:
            common_columns = optimal_positions
            min_coverage = optimal_threshold
            min_required_final = int(np.ceil(min_coverage * n_matched))
            # Update the count after optimization
            n_final_positions_before_complete = len(common_columns)
    
    coverage_pct = (len(common_columns) / n_columns) * 100
    print(f"\n   3b. Final positions (>={min_coverage*100:.0f}% coverage):")
    print(f"       - Found: {len(common_columns)} positions ({coverage_pct:.1f}% of MSA)")
    dropped_from_candidates = len(candidate_columns) - len(common_columns)
    if dropped_from_candidates > 0:
        print(f"       - Dropped from candidates: {dropped_from_candidates} positions")
        print(f"       - Reason: Structures missing coordinates at these positions")
    else:
        print(f"       - All candidate positions retained (no structures missing coordinates)")
    print(f"       - Requirement: >= {min_required_final} structures with actual coordinates (out of {n_matched} matched)")
    print(f"       - Note: Only positions where structures have actual coordinates are used")
    print(f"         (missing coordinates are excluded, not interpolated)")
    
    # Warnings for poor alignment quality
    if len(common_columns) < 10:
        print(f"\n   [WARNING] Fewer than 10 positions aligned ({len(common_columns)} positions).")
        print(f"   For highly variable domains (LRR, IDRs), use:")
        print(f"   --msa-threshold 0.5 --min-coverage 0.5")
    
    if gap_percentage > 40:
        print(f"\n   [WARNING] High gap percentage ({gap_percentage:.1f}%) detected.")
        print(f"   Consider lowering --msa-threshold and/or --min-coverage for fragmented alignments.")
    
    # Extract residues corresponding to common columns
    # CRITICAL: Only use positions where structures ACTUALLY have coordinates.
    # Do NOT interpolate missing coordinates.
    # common_columns already filtered to positions where enough structures have coordinates
    msa_to_structure_maps = {}
    homologous_structures = [None] * len(structures)  # Preserve order
    
    # Handle complete vs incomplete coverage based on allow_imputation flag
    if allow_imputation:
        print(f"\n   3c. Imputation mode (--allow-imputation):")
        if optimize_coverage:
            print(f"       - Using optimized coverage threshold: {min_coverage*100:.0f}%")
        else:
            print(f"       - Using specified coverage threshold: {min_coverage*100:.0f}%")
        print(f"       - Using {len(common_columns)} positions")
        print(f"       - Missing positions will be imputed using mean distance method")
        print(f"       - Actual imputation rate will be calculated and reported after distance feature computation")
        # Keep all common_columns (no filtering to complete coverage)
        positions_with_all_coords = common_columns
        excluded_for_incomplete = 0
    else:
        # CRITICAL: Find positions where ALL structures have coordinates
        # This ensures all structures end up with the same number of residues
        # (required for PCA which needs consistent feature dimensions)
        positions_with_all_coords = []
        for msa_col in common_columns:
            n_with_coords = 0
            for struct_idx in matched_sequences.keys():
                msa_to_struct = msa_to_structure_maps_prelim[struct_idx]
                if msa_col in msa_to_struct:
                    struct_res_idx = msa_to_struct[msa_col]
                    if struct_res_idx < structures[struct_idx].shape[0]:
                        n_with_coords += 1
            
            # Only include if ALL matched structures have coordinates at this position
            # This ensures consistent structure sizes for PCA
            if n_with_coords == len(matched_sequences):
                positions_with_all_coords.append(msa_col)
        
        # If not all positions have complete coverage, use only complete positions
        # Note: n_final_positions_before_complete was already set earlier (before optimization or here if no optimization)
        if 'n_final_positions_before_complete' not in locals():
            n_final_positions_before_complete = len(common_columns)
        excluded_for_incomplete = 0
        if len(positions_with_all_coords) < len(common_columns):
            excluded_for_incomplete = len(common_columns) - len(positions_with_all_coords)
            print(f"\n   3c. Complete coverage requirement (ALL structures must have coordinates):")
            print(f"       - Positions with complete coverage: {len(positions_with_all_coords)}")
            print(f"       - Excluded (some structures missing): {excluded_for_incomplete} positions")
            print(f"       - Reason: Some structures don't have coordinates at these positions")
            print(f"         (Likely due to structures missing different residues)")
            print(f"       - Using {len(positions_with_all_coords)} positions with complete coverage")
            print(f"         (All structures have coordinates at these positions)")
            print(f"       - Note: Use --allow-imputation to include positions with incomplete coverage")
            common_columns = positions_with_all_coords
        else:
            print(f"\n   3c. Complete coverage requirement (ALL structures must have coordinates):")
            print(f"       - All {len(common_columns)} positions have complete coverage")
            print(f"       - No positions excluded")
    
    # Summary of filtering cascade
    final_coverage_pct = (len(common_columns) / n_columns) * 100
    print(f"\n   FILTERING SUMMARY:")
    print(f"   {'='*70}")
    print(f"   Total MSA columns: {n_columns}")
    print(f"   -> Candidate columns (>={msa_threshold*100:.0f}%): {len(candidate_columns)} ({candidate_pct:.1f}%)")
    print(f"   -> Final positions (>={min_coverage*100:.0f}% coverage): {n_final_positions_before_complete} ({100*n_final_positions_before_complete/n_columns:.1f}%)")
    print(f"   -> Complete coverage (ALL structures): {len(common_columns)} ({final_coverage_pct:.1f}%)")
    print(f"\n   Total dropped:")
    print(f"   - From MSA to candidates: {n_columns - len(candidate_columns)} columns")
    print(f"   - From candidates to final: {dropped_from_candidates} positions")
    print(f"   - From final to complete: {excluded_for_incomplete} positions")
    print(f"   FINAL RESULT: {len(common_columns)} positions used for alignment")
    print(f"   {'='*70}")
    
    if len(common_columns) == 0:
        raise ValueError(
            f"No positions found where all structures have coordinates. "
            f"This may be due to structures missing different residues. "
            f"Try using --msa-threshold 0.5 --min-coverage 0.5 for highly variable domains."
        )
    
    # Extract residues - all structures will have same number of positions
    # If allow_imputation, mark missing positions as gaps (0,0,0)
    missing_masks = [] if allow_imputation else None
    
    for struct_idx in sorted(matched_sequences.keys()):
        msa_to_struct = msa_to_structure_maps_prelim[struct_idx]
        homologous_residues = []
        is_missing = [] if allow_imputation else None
        
        for msa_col in common_columns:
            # Check if structure has coordinate at this position
            has_coord = False
            if msa_col in msa_to_struct:
                struct_res_idx = msa_to_struct[msa_col]
                if struct_res_idx < structures[struct_idx].shape[0]:
                    homologous_residues.append(structures[struct_idx][struct_res_idx])
                    has_coord = True
                    if allow_imputation:
                        is_missing.append(False)
            
            if not has_coord:
                if allow_imputation:
                    # Mark as gap (will be imputed later)
                    homologous_residues.append(np.array([0.0, 0.0, 0.0]))
                    is_missing.append(True)
                else:
                    # This shouldn't happen if positions_with_all_coords is correct
                    raise ValueError(
                        f"Structure {struct_idx} missing coordinate at position {msa_col}. "
                        f"This indicates a bug in position filtering."
                    )
        
        if len(homologous_residues) == len(common_columns):
            homologous_structures[struct_idx] = np.array(homologous_residues)
            if allow_imputation:
                missing_masks.append(np.array(is_missing))
        else:
            raise ValueError(
                f"Structure {struct_idx} has {len(homologous_residues)} residues, "
                f"expected {len(common_columns)}. This indicates a bug."
            )
        
        msa_to_structure_maps[struct_idx] = msa_to_struct
    
    # Verify all structures have same number of positions
    sizes = [s.shape[0] for s in homologous_structures if s is not None]
    if len(set(sizes)) > 1:
        raise ValueError(
            f"Structures have different sizes after extraction: {set(sizes)}. "
            f"This should not happen - all structures should have {len(common_columns)} positions."
        )
    
    # For structures not in alignment, use minimum-length truncation
    # First, find the minimum size from matched structures
    matched_sizes = [s.shape[0] for s in homologous_structures if s is not None]
    if len(matched_sizes) > 0:
        min_size = min(matched_sizes)
    else:
        min_size = min(s.shape[0] for s in structures)
    
    for struct_idx in range(len(structures)):
        if homologous_structures[struct_idx] is None:
            homologous_structures[struct_idx] = structures[struct_idx][:min_size]
            print(f"   Warning: Structure {struct_idx} ({structure_names[struct_idx]}) not in alignment, using truncation to {min_size} residues")
    
    # Return missing masks if imputation is enabled
    if allow_imputation:
        return homologous_structures, msa_to_structure_maps, np.array(common_columns), missing_masks
    else:
        return homologous_structures, msa_to_structure_maps, np.array(common_columns), None

def save_coords_to_pdb(coords, filename, name="PROT"):
    """
    Save coordinate array as minimal PDB file (C-alpha atoms only).
    
    Parameters:
    -----------
    coords : array, shape (n_residues, 3)
        C-alpha coordinates
    filename : str or Path
        Output PDB file path
    name : str
        Chain identifier (default: "PROT")
    """
    with open(filename, 'w') as f:
        for i, (x, y, z) in enumerate(coords):
            # PDB format: ATOM record
            # Format: ATOM  serial  name  altLoc  resName  chainID  resSeq  iCode  x  y  z  occupancy  tempFactor
            f.write(f"ATOM  {i+1:5d}  CA  ALA {name} {i+1:4d}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           C\n")
        f.write("END\n")

def load_ca_coords_from_pdb(pdb_file):
    """
    Load C-alpha coordinates from PDB file.
    
    Parameters:
    -----------
    pdb_file : str or Path
        Path to PDB file
    
    Returns:
    --------
    coords : array, shape (n_residues, 3)
        C-alpha coordinates, or None if file cannot be read
    """
    try:
        coords = []
        with open(pdb_file, 'r') as f:
            for line in f:
                if line.startswith('ATOM') and len(line) > 15:
                    if line[13:15] == "CA":  # C-alpha atom
                        try:
                            x = float(line[30:38].strip())
                            y = float(line[38:46].strip())
                            z = float(line[46:54].strip())
                            coords.append([x, y, z])
                        except (ValueError, IndexError):
                            continue
        
        if len(coords) > 0:
            return np.array(coords)
        else:
            return None
    except Exception as e:
        print(f"Warning: Could not load coordinates from {pdb_file}: {e}")
        return None

def robust_structural_alignment_foldmason(structures, names, foldmason_executable="foldmason", threads=16):
    """
    Align ESMFold structures using FoldMason (high-speed progressive MSTA).
    
    FoldMason uses Foldseek's 3Di alphabet for fast structural alignment.
    It is highly scalable and optimized for high-speed alignment.
    """
    import tempfile
    import subprocess
    import os
    import shutil
    import threading
    import time
    import sys
    import re
    
    n_structures = len(structures)
    if n_structures < 2:
        raise ValueError("Need at least 2 structures for alignment")
    
    # Create temporary directory
    temp_dir = tempfile.mkdtemp(prefix="foldmason_")
    input_dir = os.path.join(temp_dir, "input")
    os.makedirs(input_dir, exist_ok=True)
    
    try:
        # 1. Save ESMFold coordinates as temporary PDB files
        pdb_files = []
        print(f"   Saving {n_structures} structures as temporary PDB files for FoldMason...")
        
        # Use a mapping to track which filename belongs to which structure index
        filename_to_idx = {}
        
        for i, (coords, name) in enumerate(zip(structures, names)):
            # Sanitize name but allow it to be much longer (up to 200 chars)
            # Add index i to ensure uniqueness even if names are identical
            safe_name = "".join(c if c.isalnum() or c in "._-" else "_" for c in str(name))[:200]
            filename = f"struct_{i}_{safe_name}.pdb"
            pdb_path = os.path.join(input_dir, filename)
            save_coords_to_pdb(coords, pdb_path, name="A")
            pdb_files.append(pdb_path)
            # FoldMason headers usually match the filename stem
            filename_to_idx[f"struct_{i}_{safe_name}"] = i
        
        # 2. Run FoldMason
        # Command: foldmason easy-msa <input_dir> <output_dir> <tmp_dir>
        print(f"   Running FoldMason progressive multiple structure alignment...")
        output_dir = os.path.join(temp_dir, "output")
        fm_tmp_dir = os.path.join(temp_dir, "fm_tmp")
        os.makedirs(output_dir, exist_ok=True)
        os.makedirs(fm_tmp_dir, exist_ok=True)
        
        output_prefix = os.path.join(output_dir, "aligned")
        
        cmd = [
            foldmason_executable,
            "easy-msa",
            input_dir,
            output_prefix,
            fm_tmp_dir,
            "--threads", str(threads)
        ]
        
        print(f"   FoldMason configuration:")
        print(f"     - Structures: {n_structures}")
        print(f"     - Method: FoldMason easy-msa (Foldseek-based)")
        print(f"   Estimated time: Very fast (seconds to minutes for hundreds of structures)")
        
        start_time = time.time()
        
        # Run FoldMason
        try:
            # Set a generous timeout
            timeout_seconds = 7200  # 2 hours
            
            process = subprocess.Popen(
                cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                universal_newlines=True
            )
            
            stdout_lines = []
            stderr_lines = []
            
            def read_stdout():
                for line in process.stdout:
                    stdout_lines.append(line)
                    if any(kw in line.lower() for kw in ['error', 'warning', 'alignment', 'progress']):
                        print(f"   [FoldMason] {line.strip()}")
                        sys.stdout.flush()
            
            def read_stderr():
                for line in process.stderr:
                    stderr_lines.append(line)
                    if line.strip():
                        print(f"   [FoldMason stderr] {line.strip()}")
                        sys.stdout.flush()
            
            stdout_thread = threading.Thread(target=read_stdout, daemon=True)
            stderr_thread = threading.Thread(target=read_stderr, daemon=True)
            stdout_thread.start()
            stderr_thread.start()
            
            try:
                return_code = process.wait(timeout=timeout_seconds)
            except subprocess.TimeoutExpired:
                process.kill()
                process.wait()
                raise subprocess.TimeoutExpired(cmd, timeout_seconds)
            
            stdout_thread.join(timeout=5)
            stderr_thread.join(timeout=5)
            
            if return_code != 0:
                stderr_text = ''.join(stderr_lines)
                raise subprocess.CalledProcessError(return_code, cmd, ''.join(stdout_lines), stderr_text)
            
            print(f"   FoldMason completed in {time.time() - start_time:.1f} seconds")
            
        except subprocess.CalledProcessError as e:
            print(f"   Error running FoldMason: {e}")
            raise ValueError(f"FoldMason alignment failed: {e.stderr}")
        except FileNotFoundError:
            raise ValueError(
                f"FoldMason executable not found: {foldmason_executable}\n"
                f"Please install FoldMason: conda install -c bioconda foldmason"
            )
            
        # 3. Parse alignment
        # FoldMason easy-msa produces aligned_aa.fa or aligned.fasta
        aligned_afasta = f"{output_prefix}_aa.fa"
        if not os.path.exists(aligned_afasta):
            aligned_afasta = f"{output_prefix}.fasta"
            
        if not os.path.exists(aligned_afasta):
            # Check for other possible names
            created_files = os.listdir(output_dir)
            print(f"   Output files: {created_files}")
            for f in created_files:
                if f.endswith("_aa.fa") or f.endswith(".fa") or f.endswith(".fasta"):
                    aligned_afasta = os.path.join(output_dir, f)
                    break
        
        if not os.path.exists(aligned_afasta):
            raise ValueError(f"FoldMason output FASTA file not found.")
            
        # Parse alignment from FASTA file
        print(f"   Parsing FoldMason alignment...")
        alignment_dict_afasta = {}
        current_seq_id = None
        current_seq = []
        
        with open(aligned_afasta, 'r') as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                if line.startswith('>'):
                    if current_seq_id is not None:
                        alignment_dict_afasta[current_seq_id] = ''.join(current_seq)
                    # FoldMason might include paths or .pdb in the header, clean it
                    header = line[1:].strip()
                    # Extract just the filename stem if it's a path
                    current_seq_id = os.path.basename(header).replace('.pdb', '')
                    current_seq = []
                else:
                    current_seq.append(line)
            if current_seq_id is not None:
                alignment_dict_afasta[current_seq_id] = ''.join(current_seq)
        
        # Match alignment sequences to input structures
        afasta_to_struct_idx = {}
        for afasta_id in alignment_dict_afasta.keys():
            # The afasta_id from FoldMason is usually the filename stem
            # We clean it just like we did when reading the file
            clean_id = os.path.basename(afasta_id).replace('.pdb', '')
            
            if clean_id in filename_to_idx:
                afasta_to_struct_idx[afasta_id] = filename_to_idx[clean_id]
            else:
                # Fallback fuzzy matching if headers aren't exact
                for stem, idx in filename_to_idx.items():
                    if stem in clean_id or clean_id in stem:
                        afasta_to_struct_idx[afasta_id] = idx
                        break
        
        # Check for missing structures
        matched_indices = set(afasta_to_struct_idx.values())
        missing_count = n_structures - len(matched_indices)
        if missing_count > 0:
            print(f"   Warning: {missing_count}/{n_structures} structures could not be matched to alignment.")
            for i, name in enumerate(names):
                if i not in matched_indices:
                    print(f"     - Unmatched: {name[:50]}...")
        
        if not afasta_to_struct_idx:
            raise ValueError("Could not match any FoldMason sequences to structures.")

        # Strategy: Keep ALL structures using Mean Imputation for gaps.
        # For V-domains, we want to maximize framework positions while keeping imputation low (< 5%).
        current_ids = list(afasta_to_struct_idx.keys())
        seq_length = len(next(iter(alignment_dict_afasta.values())))
        
        print(f"   Optimizing alignment for high-coverage V-domain framework...")
        
        coverage = np.zeros(seq_length)
        for sid in current_ids:
            seq = alignment_dict_afasta[sid]
            for i, char in enumerate(seq):
                if char != '-':
                    coverage[i] += 1
        
        coverage_pct = coverage / len(current_ids)
        
        # New Strategy: Prioritize high coverage (Framework) over total residue count.
        # We search for a core that keeps imputation very low.
        thresholds = [1.0, 0.99, 0.98, 0.95, 0.9, 0.8]
        best_threshold = 1.0
        n_common = 0
        
        for t in thresholds:
            mask = coverage_pct >= t
            count = np.sum(mask)
            
            # Predict imputation rate for this threshold
            # Imputed pairs roughly = total_pairs - (count_residues_at_threshold)^2
            # More accurately: track residues missing in this mask
            missing_residues_in_mask = 0
            for sid in current_ids:
                seq = alignment_dict_afasta[sid]
                for i in range(seq_length):
                    if mask[i] and seq[i] == '-':
                        missing_residues_in_mask += 1
            
            total_residue_slots = len(current_ids) * count
            predicted_rate = (missing_residues_in_mask / total_residue_slots) * 100 if total_residue_slots > 0 else 100
            
            # For V-domains, we want at least 60 residues (typical framework size)
            # and an imputation rate under 5%.
            if count >= 60 and predicted_rate <= 5.0:
                best_threshold = t
                n_common = count
                is_common = mask
                print(f"   - Found suitable core at {t*100:.0f}% coverage ({count} residues, ~{predicted_rate:.1f}% imputation)")
                break
        else:
            # Fallback to whatever gives us 60 residues if 5% is impossible
            for t in thresholds:
                mask = coverage_pct >= t
                count = np.sum(mask)
                if count >= 60:
                    best_threshold = t
                    n_common = count
                    is_common = mask
                    break
            else:
                best_threshold = thresholds[-1]
                is_common = coverage_pct >= best_threshold
                n_common = np.sum(is_common)

        print(f"   Selected {n_common} consensus positions (Framework Core, >={best_threshold*100:.0f}% coverage)")
        
        # 2. Extract coordinates and flag missing positions for imputation
        extracted_structures = []
        missing_masks = []
        final_names = []
        
        for afasta_id in current_ids:
            struct_idx = afasta_to_struct_idx[afasta_id]
            found_seq = alignment_dict_afasta[afasta_id]
            orig_coords = structures[struct_idx]
            
            struct_coords = []
            is_missing = []
            res_idx = 0
            for i in range(seq_length):
                if is_common[i]:
                    if found_seq[i] != '-':
                        struct_coords.append(orig_coords[res_idx])
                        is_missing.append(False)
                    else:
                        # Gap at a consensus position: Use a placeholder (0,0,0) for now
                        struct_coords.append(np.array([0.0, 0.0, 0.0]))
                        is_missing.append(True)
                
                if found_seq[i] != '-':
                    res_idx += 1
            
            extracted_structures.append(np.array(struct_coords))
            missing_masks.append(np.array(is_missing))
            final_names.append(names[struct_idx])
            
        # 3. Handle missing data via Distance-based Imputation
        # We fill gaps in the distance matrix rather than the coordinate matrix.
        # This is the standard "Mean Substitution" method from Polly et al. (2013).
        print(f"   Generating distance features and handling gaps...")
        
        # Calculate distance features for all pairs
        n_pairs = n_common * (n_common - 1) // 2
        distance_features_final = np.zeros((len(extracted_structures), n_pairs))
        
        # Track imputation stats
        total_elements = len(extracted_structures) * n_pairs
        imputed_count = 0
        
        for s_idx in range(len(extracted_structures)):
            coords = extracted_structures[s_idx]
            mask = missing_masks[s_idx]
            
            pair_idx = 0
            for i in range(n_common):
                for j in range(i+1, n_common):
                    if not mask[i] and not mask[j]:
                        # Both residues exist, calculate real distance
                        distance_features_final[s_idx, pair_idx] = np.linalg.norm(coords[i] - coords[j])
                    else:
                        # One or both are gaps, mark for imputation
                        distance_features_final[s_idx, pair_idx] = -1.0
                        imputed_count += 1
                    pair_idx += 1
        
        # Perform Mean Imputation on the distances
        for p_idx in range(n_pairs):
            col = distance_features_final[:, p_idx]
            valid_vals = col[col >= 0]
            if len(valid_vals) > 0:
                mean_dist = np.mean(valid_vals)
                # Replace placeholders
                distance_features_final[col < 0, p_idx] = mean_dist
            else:
                distance_features_final[:, p_idx] = 0.0 # Should not happen with high coverage

        imputation_rate = (imputed_count / total_elements) * 100
        print(f"   Distance Imputation complete:")
        print(f"   - Imputation Rate: {imputation_rate:.2f}% (Percentage of values substituted by mean)")
        print(f"   - Recommendation: < 5% is excellent, 5-10% is acceptable.")
        
        # Final dataset for return
        # For GPA, we still need coordinate arrays. We'll use the coordinate-imputed ones
        # just for the initial alignment, but the PCA will use the cleaned distance matrix.
        final_structures = [s.copy() for s in extracted_structures]
        for p in range(n_common):
            valid_coords = [final_structures[s_idx][p] for s_idx in range(len(final_structures)) if not missing_masks[s_idx][p]]
            if valid_coords:
                mean_pos = np.mean(valid_coords, axis=0)
                for s_idx in range(len(final_structures)):
                    if missing_masks[s_idx][p]:
                        final_structures[s_idx][p] = mean_pos

        print(f"   Final dataset: {len(final_structures)} structures with {n_common} common positions")
        
        # Statistics for the log
        alignment_stats = {
            'n_positions': n_common,
            'avg_pairwise_rmsd': 0.0,
            'n_structures_aligned': len(final_structures)
        }
        
        return final_structures, is_common, alignment_stats, final_names, distance_features_final
        
    finally:
        try:
            shutil.rmtree(temp_dir)
        except:
            pass

def find_optimal_correspondence(reference, target, max_distance=5.0, sequence_weight=0.1):
    """
    Find optimal one-to-one correspondence between reference and target structures.
    
    Uses a combination of:
    1. Structural distance (spatial proximity)
    2. Sequence order preference (maintains sequential ordering when possible)
    
    Parameters:
    -----------
    reference : array, shape (n_ref, 3)
        Reference structure coordinates
    target : array, shape (n_target, 3)
        Target structure coordinates
    max_distance : float
        Maximum distance (Å) for considering residues as corresponding
    sequence_weight : float
        Weight for sequence order preference (0.0 = distance only, higher = more sequence order preference)
    
    Returns:
    --------
    correspondence : list
        List of target indices corresponding to each reference residue, or None if no match
    """
    n_ref = reference.shape[0]
    n_target = target.shape[0]
    
    # Compute distance matrix
    # distances[i, j] = distance between reference[i] and target[j]
    distances = np.linalg.norm(reference[:, np.newaxis, :] - target[np.newaxis, :, :], axis=2)
    
    # Add sequence order penalty: prefer matches that maintain sequence order
    # Penalty increases with deviation from diagonal
    sequence_penalty = np.zeros_like(distances)
    for i in range(n_ref):
        for j in range(n_target):
            # Expected position if sequences were perfectly aligned
            expected_j = int(j * n_ref / n_target) if n_target > 0 else j
            # Penalty increases with distance from expected position
            sequence_penalty[i, j] = abs(i - expected_j) * sequence_weight
    
    # Combined cost: distance + sequence penalty
    cost_matrix = distances + sequence_penalty
    
    # Set distances > max_distance to a very high cost (effectively infinite)
    cost_matrix[distances > max_distance] = 1e10
    
    # Use Hungarian algorithm (optimal assignment) for one-to-one matching
    # This finds the globally optimal assignment minimizing total cost
    # while ensuring each target residue matches at most one reference residue
    
    # Hungarian algorithm works on square matrices, so we need to pad if necessary
    max_size = max(n_ref, n_target)
    
    if n_ref != n_target:
        # Pad cost matrix to square size with high cost for unmatched positions
        padded_cost = np.full((max_size, max_size), 1e10)
        padded_cost[:n_ref, :n_target] = cost_matrix
    else:
        padded_cost = cost_matrix
    
    # Find optimal assignment using Hungarian algorithm
    row_indices, col_indices = linear_sum_assignment(padded_cost)
    
    # Build correspondence map (only for valid reference residues)
    correspondence = [None] * n_ref
    for ref_idx, target_idx in zip(row_indices, col_indices):
        # Only process if within valid ranges and distance threshold
        if ref_idx < n_ref and target_idx < n_target:
            if distances[ref_idx, target_idx] <= max_distance:
                correspondence[ref_idx] = target_idx
    
    return correspondence

def structural_alignment_find_corresponding_residues(structures, reference_idx=0, max_distance=5.0, use_optimal_matching=True):
    """
    Find corresponding residues across structures using structural alignment.
    
    This function aligns structures to a reference and finds corresponding residues
    based on spatial proximity. Uses improved matching algorithm that ensures
    one-to-one correspondence and maintains sequence order when possible.
    
    Parameters:
    -----------
    structures : list of arrays, each shape (n_residues, 3)
        List of protein structures
    reference_idx : int
        Index of reference structure to align all others to
    max_distance : float
        Maximum distance (Angstroms) for considering residues as corresponding
    use_optimal_matching : bool, default=True
        If True, uses improved matching with one-to-one correspondence constraint.
        If False, uses simple nearest-neighbor matching (legacy).
    
    Returns:
    --------
    aligned_structures : list of arrays, each shape (n_common_residues, 3)
        Structures with corresponding residues extracted
    correspondence_map : list of arrays
        Maps original residue indices to common alignment positions
    """
    n_structures = len(structures)
    reference = structures[reference_idx].copy()
    
    # Align all structures to reference using Procrustes
    aligned_to_ref = []
    for i, struct in enumerate(structures):
        if i == reference_idx:
            aligned_to_ref.append(reference)
        else:
            aligned, _, _, _ = procrustes_superimpose(reference, struct, scale=False)
            aligned_to_ref.append(aligned)
    
    # Find corresponding residues using improved matching
    n_ref = reference.shape[0]
    correspondence = []  # correspondence[i][j] = index in structure i corresponding to ref residue j
    
    for struct_idx, aligned_struct in enumerate(aligned_to_ref):
        struct_correspondence = []
        
        if struct_idx == reference_idx:
            # Reference structure: identity mapping
            struct_correspondence = list(range(n_ref))
        else:
            if use_optimal_matching:
                # Use improved matching with one-to-one correspondence constraint
                struct_correspondence = find_optimal_correspondence(
                    reference, aligned_struct, 
                    max_distance=max_distance,
                    sequence_weight=0.1  # Small weight to prefer sequence order when distances are similar
                )
            else:
                # Legacy: simple nearest-neighbor matching (may have many-to-one issues)
                for ref_res_idx in range(n_ref):
                    ref_pos = reference[ref_res_idx]
                    distances = np.linalg.norm(aligned_struct - ref_pos, axis=1)
                    closest_idx = np.argmin(distances)
                    
                    # Only include if within max_distance
                    if distances[closest_idx] <= max_distance:
                        struct_correspondence.append(closest_idx)
                    else:
                        struct_correspondence.append(None)  # No corresponding residue
        
        correspondence.append(struct_correspondence)
    
    # Find residues that have correspondences in ALL structures
    n_common = n_ref
    common_positions = []
    
    for ref_pos in range(n_ref):
        # Check if this position has correspondence in all structures
        has_all = all(correspondence[i][ref_pos] is not None for i in range(n_structures))
        if has_all:
            common_positions.append(ref_pos)
    
    if len(common_positions) == 0:
        # Fallback: use minimum size if no common positions found
        print(f"   Warning: No common positions found with max_distance={max_distance}Å")
        print(f"   Falling back to minimum size approach")
        min_size = min(s.shape[0] for s in structures)
        return [s[:min_size, :] for s in structures], None
    
    # Extract corresponding residues
    aligned_structures = []
    for struct_idx in range(n_structures):
        extracted = []
        for ref_pos in common_positions:
            orig_idx = correspondence[struct_idx][ref_pos]
            extracted.append(aligned_to_ref[struct_idx][orig_idx])
        aligned_structures.append(np.array(extracted))
    
    return aligned_structures, correspondence

def structural_alignment_improved(structures, reference_idx=0, distance_threshold=3.0, preserve_positions=False):
    """
    Improved structural alignment using iterative refinement.
    
    This approach:
    1. Aligns all structures to a reference
    2. Finds corresponding residues based on spatial proximity
    3. Iteratively refines alignment using common residues
    
    Parameters:
    -----------
    structures : list of arrays
        Protein structures
    reference_idx : int
        Index of reference structure
    distance_threshold : float
        Distance threshold for residue correspondence (Angstroms)
    preserve_positions : bool, default=False
        If True, preserve all positions (for MSA-aligned structures).
        If False, may truncate to common positions.
    
    Returns:
    --------
    aligned_structures : list of arrays with same shape
        Aligned structures with corresponding residues
    """
    n_structures = len(structures)
    
    # Check if structures already have same size (likely from MSA alignment)
    structure_sizes = [s.shape[0] for s in structures]
    all_same_size = len(set(structure_sizes)) == 1
    
    # If structures have same size and preserve_positions is True, just refine alignment
    # without truncating (preserves MSA-selected positions)
    if preserve_positions and all_same_size:
        n_residues = structure_sizes[0]
        reference = structures[reference_idx].copy()
        
        # Align all structures to reference
        aligned = []
        for i, struct in enumerate(structures):
            if i == reference_idx:
                aligned.append(reference)
            else:
                aligned_struct, _, _, _ = procrustes_superimpose(reference, struct, scale=False)
                aligned.append(aligned_struct)
        
        # Iterative refinement: realign to mean, preserving all positions
        for iteration in range(5):
            # Calculate mean structure
            mean_structure = np.mean(aligned, axis=0)
            
            # Realign each structure to mean
            for struct_idx in range(n_structures):
                aligned[struct_idx], _, _, _ = procrustes_superimpose(mean_structure, aligned[struct_idx], scale=False)
            
            # Check convergence
            if iteration > 0:
                if np.linalg.norm(mean_structure - old_mean) < 0.01:
                    break
            old_mean = mean_structure.copy()
        
        return aligned
    
    # Find common positions
    reference = structures[reference_idx].copy()
    
    # Initial alignment to reference
    # First, ensure all structures have at least the reference size for alignment
    # We'll align using the minimum size, then find correspondences
    min_size = min(s.shape[0] for s in structures)
    reference_aligned = reference[:min_size, :] if reference.shape[0] > min_size else reference
    
    aligned = []
    for i, struct in enumerate(structures):
        if i == reference_idx:
            aligned.append(reference_aligned)
        else:
            # Use minimum size for alignment
            struct_aligned = struct[:min_size, :] if struct.shape[0] > min_size else struct
            aligned_struct, _, _, _ = procrustes_superimpose(reference_aligned, struct_aligned, scale=False)
            aligned.append(aligned_struct)
    
    reference = reference_aligned
    
    # Iterative refinement: find common residues, realign, repeat
    for iteration in range(5):
        # Find common residues based on distance
        n_ref = reference.shape[0]
        common_mask = np.ones(n_ref, dtype=bool)
        
        for struct_idx, aligned_struct in enumerate(aligned):
            if struct_idx == reference_idx:
                continue
            
            # For each reference residue, check if there's a close match
            struct_mask = np.zeros(n_ref, dtype=bool)
            for ref_idx in range(n_ref):
                ref_pos = reference[ref_idx]
                distances = np.linalg.norm(aligned_struct - ref_pos, axis=1)
                min_dist = np.min(distances)
                struct_mask[ref_idx] = min_dist <= distance_threshold
            
            common_mask = common_mask & struct_mask
        
        n_common = np.sum(common_mask)
        
        if n_common < 10:  # Need at least 10 common residues
            print(f"   Warning: Only {n_common} common residues found, using minimum size")
            min_size = min(s.shape[0] for s in structures)
            return [s[:min_size, :] for s in structures]
        
        # Extract common residues
        common_ref = reference[common_mask]
        common_structures = []
        for struct_idx, aligned_struct in enumerate(aligned):
            if struct_idx == reference_idx:
                common_structures.append(common_ref)
            else:
                # Find closest residues in this structure
                common_struct = []
                for ref_pos in common_ref:
                    distances = np.linalg.norm(aligned_struct - ref_pos, axis=1)
                    closest_idx = np.argmin(distances)
                    common_struct.append(aligned_struct[closest_idx])
                common_structures.append(np.array(common_struct))
        
        # Realign using common residues
        new_reference = np.mean(common_structures, axis=0)
        for struct_idx in range(n_structures):
            if struct_idx == reference_idx:
                aligned[struct_idx] = new_reference
            else:
                aligned[struct_idx], _, _, _ = procrustes_superimpose(new_reference, common_structures[struct_idx], scale=False)
        
        # Check convergence
        if iteration > 0:
            if np.linalg.norm(new_reference - old_reference) < 0.01:
                break
        old_reference = new_reference.copy()
        reference = new_reference
    
    return aligned

def procrustes_anova(structures_before, structures_after, mean_structure, n_permutations=1000):
    """
    Procrustes ANOVA: Test if GPA alignment explains significant shape variation.
    
    This is a sanity check to validate that GPA superposition significantly
    reduces shape variation compared to unaligned structures.
    
    Parameters:
    -----------
    structures_before : list of arrays
        Structures before GPA alignment
    structures_after : list of arrays
        Structures after GPA alignment
    mean_structure : array
        Mean (consensus) structure from GPA
    n_permutations : int
        Number of permutations for significance test
    
    Returns:
    --------
    results : dict
        Dictionary with ANOVA statistics
    """
    n_structures = len(structures_after)
    n_residues = structures_after[0].shape[0]
    
    # Calculate Procrustes distances (sum of squared distances to mean)
    # Before alignment: distances from original structures to their individual means
    ss_before = 0.0
    for struct in structures_before:
        struct_centered = struct - struct.mean(axis=0)
        ss_before += np.sum(struct_centered ** 2)
    
    # After alignment: distances from aligned structures to consensus mean
    ss_after = 0.0
    for struct in structures_after:
        ss_after += np.sum((struct - mean_structure) ** 2)
    
    # Total sum of squares (variation in aligned structures)
    ss_total = ss_after
    
    # Variation explained by alignment
    ss_explained = ss_before - ss_after
    
    # Proportion of variation explained
    if ss_before > 0:
        prop_explained = ss_explained / ss_before
    else:
        prop_explained = 0.0
    
    # Simple F-statistic (ratio of explained to residual variation)
    # Degrees of freedom: 
    # - Explained: 6 (translation: 3, rotation: 3) per structure
    # - Residual: n_residues * 3 - 6 per structure
    df_explained = n_structures * 6  # 6 parameters per structure (3 translation + 3 rotation)
    df_residual = n_structures * (n_residues * 3 - 6)
    
    if df_residual > 0 and ss_after > 0:
        ms_explained = ss_explained / df_explained if df_explained > 0 else 0
        ms_residual = ss_after / df_residual
        f_statistic = ms_explained / ms_residual if ms_residual > 0 else np.inf
    else:
        f_statistic = np.inf
    
    # Permutation test: randomly permute structures and calculate F-statistic
    # This tests if the observed F is significantly larger than expected by chance
    f_permuted = []
    np.random.seed(42)  # For reproducibility
    for perm in range(n_permutations):
        # Randomly shuffle structures (maintains within-structure variation)
        permuted_after = structures_after.copy()
        np.random.shuffle(permuted_after)
        
        # Calculate new mean
        permuted_mean = np.mean(permuted_after, axis=0)
        
        # Calculate SS for permuted data
        ss_permuted = 0.0
        for struct in permuted_after:
            ss_permuted += np.sum((struct - permuted_mean) ** 2)
        
        # Calculate F for permuted data
        ss_explained_perm = ss_before - ss_permuted
        ms_explained_perm = ss_explained_perm / df_explained if df_explained > 0 else 0
        ms_residual_perm = ss_permuted / df_residual if df_residual > 0 else 0
        f_perm = ms_explained_perm / ms_residual_perm if ms_residual_perm > 0 else 0
        f_permuted.append(f_perm)
    
    # P-value: proportion of permuted F-statistics >= observed F
    p_value = np.mean(np.array(f_permuted) >= f_statistic)
    
    results = {
        'ss_before': ss_before,
        'ss_after': ss_after,
        'ss_explained': ss_explained,
        'prop_explained': prop_explained,
        'f_statistic': f_statistic,
        'p_value': p_value,
        'df_explained': df_explained,
        'df_residual': df_residual
    }
    
    return results

def compute_distance_features(structures, max_distance=None, sequence_separation=None):
    """
    Compute Cα-Cα distance features instead of raw coordinates.
    
    This is more robust to global alignment issues and captures
    local geometry better. Appropriate for Procrustes-aligned structures.
    
    Parameters:
    -----------
    structures : list of arrays, each shape (n_residues, 3)
        List of protein structures (should be Procrustes-aligned)
    max_distance : float, optional
        Maximum distance (Angstroms) to include. If None, includes all pairs.
        Useful for reducing feature count in large proteins.
    sequence_separation : int, optional
        Maximum sequence separation (residue index difference) to include.
        If None, includes all pairs. Useful for focusing on local structure.
        Note: If both max_distance and sequence_separation are set, both conditions must be met.
    
    Returns:
    --------
    distance_features : array, shape (n_structures, n_pairs)
        where n_pairs = number of residue pairs meeting criteria
    pair_indices : list of tuples, optional
        List of (i, j) residue index pairs included in features
        (only returned if max_distance or sequence_separation is specified)
    """
    n_structures = len(structures)
    n_residues = structures[0].shape[0]
    
    # Determine which pairs to include
    if max_distance is None and sequence_separation is None:
        # Include all pairs (original behavior)
        n_pairs = n_residues * (n_residues - 1) // 2
        pair_indices = None
        use_filtering = False
    else:
        # Pre-compute which pairs to include
        pair_indices = []
        for i in range(n_residues):
            for j in range(i+1, n_residues):
                include = True
                
                # Check sequence separation
                if sequence_separation is not None:
                    if abs(j - i) > sequence_separation:
                        include = False
                
                # Check distance (using first structure as reference)
                if include and max_distance is not None:
                    dist = np.linalg.norm(structures[0][i] - structures[0][j])
                    if dist > max_distance:
                        include = False
                
                if include:
                    pair_indices.append((i, j))
        
        n_pairs = len(pair_indices)
        use_filtering = True
        print(f"   Filtering: {n_pairs} pairs (from {n_residues * (n_residues - 1) // 2} total)")
        if max_distance is not None:
            print(f"   Max distance threshold: {max_distance} Å")
        if sequence_separation is not None:
            print(f"   Max sequence separation: {sequence_separation} residues")
    
    distance_features = np.zeros((n_structures, n_pairs))
    
    for struct_idx, structure in enumerate(structures):
        if use_filtering:
            for pair_idx, (i, j) in enumerate(pair_indices):
                dist = np.linalg.norm(structure[i] - structure[j])
                distance_features[struct_idx, pair_idx] = dist
        else:
            # Original behavior: all pairs
            pair_idx = 0
            for i in range(n_residues):
                for j in range(i+1, n_residues):
                    dist = np.linalg.norm(structure[i] - structure[j])
                    distance_features[struct_idx, pair_idx] = dist
                    pair_idx += 1
    
    if use_filtering:
        return distance_features, pair_indices
    else:
        return distance_features

def build_phylogenetic_covariance_matrix(tree, tip_names):
    """
    Build phylogenetic covariance matrix C from tree.
    
    Under Brownian motion, C[i,j] = shared branch length between tips i and j.
    C[i,i] = total branch length from tip i to root.
    
    Parameters:
    -----------
    tree : Bio.Phylo tree object or str
        Phylogenetic tree
    tip_names : list
        List of tip names in order
    
    Returns:
    --------
    C : array, shape (n_tips, n_tips)
        Phylogenetic covariance matrix
    """
    if isinstance(tree, str):
        if Phylo is None:
            raise ImportError("BioPython required for tree parsing")
        tree = Phylo.read(tree, 'newick')
    
    n_tips = len(tip_names)
    C = np.zeros((n_tips, n_tips))
    
    # Get all terminal nodes
    terminals = {term.name: term for term in tree.get_terminals()}
    
    # Normalize names: remove dots, convert pipes to match format
    def normalize_name(name):
        """Normalize protein name for matching."""
        # Remove version numbers (e.g., .1, .2)
        # Note: re is imported at module level
        name = re.sub(r'\.(\d+)', r'\1', name)  # XP_042112342.1 -> XP_0421123421
        # Replace pipes with nothing (tree has |, NPZ doesn't)
        name = name.replace('|', '')
        # Remove colons (alignment files may have colons)
        name = name.replace(':', '')
        # Remove spaces and normalize underscores
        name = name.replace(' ', '_')
        name = re.sub(r'_{2,}', '_', name)  # Multiple underscores to single
        return name
    
    # Create normalized terminal lookup
    normalized_terminals = {}
    for term in tree.get_terminals():
        norm_name = normalize_name(term.name)
        normalized_terminals[norm_name] = term
    
    # Also create a more flexible lookup using key parts of names
    # Extract accession numbers and species names for better matching
    def extract_key_parts(name):
        """Extract key identifying parts of a name for flexible matching."""
        parts = []
        # Extract XP_/NP_/XP_ accession numbers
        acc_match = re.search(r'(XP_|NP_|XP_)(\d+)', name)
        if acc_match:
            parts.append(acc_match.group(0).replace('.', ''))
        # Extract species name (first few words before underscore)
        species_parts = name.split('_')[:3]  # First 3 parts usually species
        if len(species_parts) >= 2:
            parts.append('_'.join(species_parts[:2]))
        return parts
    
    # Create flexible lookup
    flexible_lookup = {}
    for term in tree.get_terminals():
        norm_name = normalize_name(term.name)
        key_parts = extract_key_parts(norm_name)
        for part in key_parts:
            if part not in flexible_lookup:
                flexible_lookup[part] = []
            flexible_lookup[part].append(term)
    
    # For each pair of tips, find their most recent common ancestor (MRCA)
    # and calculate shared branch length
    for i, tip1_name in enumerate(tip_names):
        norm_tip1_name = normalize_name(tip1_name)
        
        tip1 = None
        # Try exact normalized match first
        if norm_tip1_name in normalized_terminals:
            tip1 = normalized_terminals[norm_tip1_name]
        # Try exact match with original name
        elif tip1_name in terminals:
            tip1 = terminals[tip1_name]
        # Try flexible matching using key parts
        else:
            key_parts = extract_key_parts(norm_tip1_name)
            for part in key_parts:
                if part in flexible_lookup:
                    # If multiple matches, try to find best one
                    candidates = flexible_lookup[part]
                    if len(candidates) == 1:
                        tip1 = candidates[0]
                        break
                    else:
                        # Try to find best match by checking if norm_tip1_name contains candidate name or vice versa
                        for candidate in candidates:
                            cand_norm = normalize_name(candidate.name)
                            if norm_tip1_name in cand_norm or cand_norm in norm_tip1_name:
                                tip1 = candidate
                                break
                        if tip1:
                            break
                if tip1:
                    break
        
        # Final fallback: try substring matching
        if tip1 is None:
            for norm_name, term in normalized_terminals.items():
                # Check if significant portion of name matches
                if len(norm_tip1_name) > 20 and len(norm_name) > 20:
                    # Check if last 30 characters match (accession + domain info)
                    if norm_tip1_name[-30:] in norm_name or norm_name[-30:] in norm_tip1_name:
                        tip1 = term
                        break
                # Or check if first part (species) matches
                tip1_parts = norm_tip1_name.split('_')[:3]
                name_parts = norm_name.split('_')[:3]
                if len(tip1_parts) >= 2 and len(name_parts) >= 2:
                    if tip1_parts[0] == name_parts[0] and tip1_parts[1] == name_parts[1]:
                        # Species matches, check if accession matches
                        if any(part in norm_name for part in tip1_parts[2:] if len(part) > 5):
                            tip1 = term
                            break
        
        if tip1 is None:
            # Only print first few warnings to avoid spam
            if i < 5:
                print(f"Warning: {tip1_name} not found in tree (tried multiple matching strategies)")
            continue
        
        # Distance from tip1 to root
        root_to_tip1 = sum([node.branch_length or 0.0 
                           for node in tree.get_path(tip1)])
        C[i, i] = root_to_tip1
        
        for j, tip2_name in enumerate(tip_names[i+1:], start=i+1):
            norm_tip2_name = normalize_name(tip2_name)
            
            if norm_tip2_name in normalized_terminals:
                tip2 = normalized_terminals[norm_tip2_name]
            elif tip2_name in terminals:
                tip2 = terminals[tip2_name]
            else:
                # Try to find partial match
                tip2 = None
                for norm_name, term in normalized_terminals.items():
                    if tip2_name in norm_name or norm_name in tip2_name:
                        tip2 = term
                        break
                if tip2 is None:
                    continue
            
            # Find MRCA
            path1 = tree.get_path(tip1)
            path2 = tree.get_path(tip2)
            
            # Find common ancestor
            common_ancestors = set(path1) & set(path2)
            if common_ancestors:
                mrca = max(common_ancestors, key=lambda n: tree.get_path(n).__len__())
                # Distance from MRCA to root
                root_to_mrca = sum([node.branch_length or 0.0 
                                   for node in tree.get_path(mrca)])
                shared_length = root_to_mrca
            else:
                shared_length = 0.0
            
            C[i, j] = shared_length
            C[j, i] = shared_length
    
    return C

def standard_pca(X):
    """
    Standard Principal Components Analysis.
    
    Uses SVD-based approach when n_features > n_samples for numerical stability.
    
    Parameters:
    -----------
    X : array, shape (n_samples, n_features)
        Data matrix (mean-centered)
    
    Returns:
    --------
    scores : array, shape (n_samples, n_components)
        PCA scores
    eigenvectors : array, shape (n_features, n_components)
        Principal component axes
    eigenvalues : array, shape (n_components,)
        Eigenvalues
    """
    n_samples, n_features = X.shape
    n = n_samples
    
    # Use SVD-based approach when n_features >> n_samples (more numerically stable)
    # This avoids computing the large covariance matrix X.T @ X
    if n_features > n_samples:
        # SVD on X directly: X = U @ S @ V.T
        # This is more efficient and numerically stable when n_features >> n_samples
        U, s, Vt = np.linalg.svd(X, full_matrices=False)
        
        # Eigenvalues are s^2 / (n-1)
        eigenvalues = (s ** 2) / (n - 1)
        
        # Eigenvectors are V (right singular vectors)
        eigenvectors = Vt.T
        
        # Scores are U @ S (already computed)
        scores = U @ np.diag(s)
        
        # Sort by eigenvalue (descending) - already sorted by SVD
        # eigenvalues are already in descending order
        # eigenvectors are already in corresponding order
        
    else:
        # Standard approach: compute covariance matrix and eigendecompose
        P = (1.0 / (n - 1)) * (X.T @ X)
        
        # Eigen decomposition
        eigenvalues, eigenvectors = np.linalg.eigh(P)
        
        # Sort by eigenvalue (descending)
        idx = np.argsort(eigenvalues)[::-1]
        eigenvalues = eigenvalues[idx]
        eigenvectors = eigenvectors[:, idx]
        
        # Project data
        scores = X @ eigenvectors
    
    return scores, eigenvectors, eigenvalues

def phylogenetic_pca(X, C):
    """
    Phylogenetic Principal Components Analysis (Revell 2009).
    
    Parameters:
    -----------
    X : array, shape (n_samples, n_features)
        Data matrix (Procrustes coordinates)
    C : array, shape (n_samples, n_samples)
        Phylogenetic covariance matrix
    
    Returns:
    --------
    scores : array, shape (n_samples, n_components)
        pPCA scores
    eigenvectors : array, shape (n_features, n_components)
        pPCA axes
    eigenvalues : array, shape (n_components,)
        pPCA eigenvalues
    ancestral_root : array, shape (n_features,)
        Estimated ancestral root values
    """
    n = X.shape[0]
    m = X.shape[1]
    
    # Step 1: Estimate ancestral root node values (Equation 3)
    ones = np.ones(n)
    C_inv = np.linalg.pinv(C)  # Use pseudo-inverse in case C is singular
    
    # Ancestral root estimate
    denominator = ones.T @ C_inv @ ones
    if denominator == 0:
        # Fallback to mean if C is problematic
        ancestral_root = X.mean(axis=0)
    else:
        ancestral_root = ((1.0 / denominator) * ones.T @ C_inv @ X).reshape(m)
    
    # Step 2: Calculate evolutionary covariance matrix (Equation 4)
    X_centered = X - ancestral_root
    
    # Use SVD-based approach when n_features >> n_samples (more numerically stable)
    # This avoids computing the large covariance matrix X.T @ C_inv @ X
    n_samples, n_features = X_centered.shape
    
    if n_features > n_samples:
        # SVD-based approach: Transform to standard PCA form
        # P_P = (1/(n-1)) * X.T @ C_inv @ X
        # We can use SVD on: Y = sqrt(C_inv) @ X
        # But computing sqrt(C_inv) is expensive, so we use a different approach:
        # Compute Y = C_inv^(1/2) @ X using Cholesky decomposition or eigendecomposition
        
        # Method: Use eigendecomposition of C_inv to get sqrt
        # C_inv = Q @ Lambda @ Q.T, so sqrt(C_inv) = Q @ sqrt(Lambda) @ Q.T
        eigenvals_C, eigenvecs_C = np.linalg.eigh(C_inv)
        # Ensure positive eigenvalues (C_inv should be positive definite)
        eigenvals_C = np.maximum(eigenvals_C, 1e-10)
        sqrt_C_inv = eigenvecs_C @ np.diag(np.sqrt(eigenvals_C)) @ eigenvecs_C.T
        
        # Transform data: Y = sqrt(C_inv) @ X_centered
        Y = sqrt_C_inv @ X_centered
        
        # Now use SVD on Y (same as standard PCA)
        U, s, Vt = np.linalg.svd(Y, full_matrices=False)
        
        # Eigenvalues are s^2 / (n-1)
        eigenvalues = (s ** 2) / (n - 1)
        
        # Eigenvectors are V (right singular vectors)
        eigenvectors = Vt.T
        
        # Scores: Project original X_centered onto eigenvectors
        scores = X_centered @ eigenvectors
    else:
        # Standard approach: compute covariance matrix and eigendecompose
        P_P = (1.0 / (n - 1)) * (X_centered.T @ C_inv @ X_centered)
        
        # Eigen decomposition of P_P
        eigenvalues, eigenvectors = np.linalg.eigh(P_P)
        
        # Sort by eigenvalue (descending)
        idx = np.argsort(eigenvalues)[::-1]
        eigenvalues = eigenvalues[idx]
        eigenvectors = eigenvectors[:, idx]
        
        # Step 4: Project tip taxa into pPCA space (Equation 5)
        scores = X_centered @ eigenvectors
    
    return scores, eigenvectors, eigenvalues, ancestral_root

def main():
    """Main analysis pipeline."""
    import argparse
    
    parser = argparse.ArgumentParser(
        description='Phylogenetic PCA for ESMFold Protein Structures'
    )
    parser.add_argument(
        '--npz',
        type=str,
        required=True,
        help='Path to ESMFold coordinates NPZ file'
    )
    parser.add_argument(
        '--tree',
        type=str,
        required=True,
        help='Path to phylogenetic tree Newick file'
    )
    parser.add_argument(
        '--output-dir',
        type=str,
        default=None,
        help='Output directory name (default: inferred from input file)'
    )
    parser.add_argument(
        '--use-coordinate-pca',
        action='store_true',
        help='Use coordinate-based PCA instead of distance-based (NOT RECOMMENDED for Procrustes-aligned structures)'
    )
    parser.add_argument(
        '--alignment',
        type=str,
        required=False,
        default=None,
        help='Path to sequence alignment file (.aln format). If provided, uses MSA-based alignment with automatic '
             'coverage optimization and imputation (DEFAULT METHOD). This maximizes positions while minimizing '
             'imputation rate. If not provided, uses FoldMason structural alignment.'
    )
    parser.add_argument(
        '--use-structural-alignment-only',
        action='store_true',
        default=False,
        help='Use FoldMason structural alignment only (skip MSA). Finds correspondences based on 3D structure similarity. '
             'Useful for highly variable domains where sequence alignment fails. '
             'Automatically enabled if --alignment is not provided.'
    )
    parser.add_argument(
        '--no-imputation',
        action='store_true',
        default=False,
        help='Explicitly disable imputation for MSA-based alignment (this is the DEFAULT behavior). '
             'Requires 100%% coverage (all structures must have coordinates at each position). '
             'This flag is redundant since complete coverage is the default, but can be used for clarity.'
    )
    parser.add_argument(
        '--structural-refinement',
        action='store_true',
        default=True,
        help='Apply structural alignment refinement after sequence-based correspondence (default: True). '
             'When MSA alignment is good, you can skip refinement with --no-structural-refinement.'
    )
    parser.add_argument(
        '--no-structural-refinement',
        action='store_false',
        dest='structural_refinement',
        help='Skip structural alignment refinement after MSA extraction. '
             'Use this when MSA alignment is good and you want to preserve all MSA-selected positions.'
    )
    parser.add_argument(
        '--auto-skip-refinement',
        action='store_true',
        default=False,
        help='Automatically skip structural refinement if MSA alignment quality is good. '
             'Good quality = high coverage (≥80%%), all structures aligned, ≥50 positions. '
             'If quality is not good, refinement will still be applied.'
    )
    parser.add_argument(
        '--use-foldmason-refinement',
        action='store_true',
        default=False,
        help='Use FoldMason for structural refinement after MSA (instead of Procrustes). '
             'This allows imputation of missing positions, potentially capturing more positions. '
             'MSA identifies candidate positions (evolutionary guide), FoldMason verifies structural correspondence and handles gaps.'
    )
    parser.add_argument(
        '--max-distance',
        type=float,
        default=None,
        help='Maximum distance (Å) for distance features. Reduces feature count for large proteins. (default: include all pairs)'
    )
    parser.add_argument(
        '--max-sequence-separation',
        type=int,
        default=None,
        help='Maximum sequence separation (residue index difference) for distance features. Focuses on local structure. (default: include all pairs)'
    )
    parser.add_argument(
        '--correspondence-distance-threshold',
        type=float,
        default=5.0,
        help='Distance threshold (Å) for structural correspondence during refinement (default: 5.0).'
    )
    parser.add_argument(
        '--foldmason-executable',
        type=str,
        default='foldmason',
        help='Path to FoldMason executable (default: "foldmason" - assumes in PATH).'
    )
    parser.add_argument(
        '--threads',
        type=int,
        default=16,
        help='Number of threads for FoldMason alignment (default: 16).'
    )
    parser.add_argument(
        '--external-alignment',
        type=str,
        default=None,
        help='Path to directory containing pre-aligned PDB files from external tool. '
             'If provided, uses these aligned structures instead of internal structural alignment.'
    )
    parser.add_argument(
        '--msa-threshold',
        type=float,
        default=0.8,
        help='Minimum proportion of sequences that must have residues at an MSA column to be considered as candidate (default: 0.8 = 80%%). Good for most proteins. For highly variable domains (e.g., LRR), try 0.5-0.6.'
    )
    parser.add_argument(
        '--min-coverage',
        type=float,
        default=0.7,
        help='Minimum proportion of structures that must have residues at a column to be included in final alignment (default: 0.7 = 70%%). Allows some missing structures. For highly variable domains, try 0.5.'
    )
    parser.add_argument(
        '--allow-imputation',
        action='store_true',
        default=None,  # None means auto-detect: False for MSA (default), False for structural-only
        help='Allow imputation for positions with incomplete coverage (instead of requiring 100%% coverage). '
             'Uses mean distance imputation (Polly et al. 2013). When enabled, automatically optimizes coverage '
             'threshold to find best balance between number of positions and imputation rate. '
             'DEFAULT: Disabled (requires 100%% coverage). Use this flag to enable imputation. '
             'Imputation rate will be reported in output.'
    )
    parser.add_argument(
        '--no-optimize-coverage',
        action='store_true',
        default=False,
        help='Disable automatic coverage threshold optimization when --allow-imputation is used. '
             'Use the --min-coverage value directly instead of optimizing.'
    )
    parser.add_argument(
        '--pdb-dir',
        type=str,
        default=None,
        help='Optional directory containing PDB files (for reference only, not used for filtering).'
    )
    args = parser.parse_args()
    
    # File paths
    script_dir = Path(__file__).parent  # ppca/scripts/
    ppca_dir = script_dir.parent  # ppca/
    
    npz_file = Path(args.npz)
    tree_file = Path(args.tree)
    
    # Determine output directory name
    if args.output_dir:
        output_dir_path = Path(args.output_dir)
        # If absolute path provided, use it directly
        if output_dir_path.is_absolute():
            csv_dir = output_dir_path
            # For figures, create a sibling directory
            fig_dir = output_dir_path.parent / f"{output_dir_path.name}_ppca_results_figures"
            # For NPZ file, use parent directory with a simpler name
            # Extract base name without the _ppca_results_csv suffix if present
            base_name = output_dir_path.name
            if base_name.endswith('_ppca_results_csv'):
                base_name = base_name.replace('_ppca_results_csv', '')
            output_file = output_dir_path.parent / f"{base_name}_ppca_results.npz"
        else:
            # Relative path: check if it's already a full directory path or just a suffix
            output_dir_path = Path(args.output_dir)
            # If it contains path separators or already ends with _ppca_results_csv, treat as full path
            if '/' in str(args.output_dir) or '\\' in str(args.output_dir) or output_dir_path.name.endswith('_ppca_results_csv'):
                # It's a full relative path, use it directly
                csv_dir = output_dir_path
                fig_dir = output_dir_path.parent / f"{output_dir_path.name}_ppca_results_figures"
                base_name = output_dir_path.name
                if base_name.endswith('_ppca_results_csv'):
                    base_name = base_name.replace('_ppca_results_csv', '')
                output_file = output_dir_path.parent / f"{base_name}_ppca_results.npz"
            else:
                # It's just a suffix, append to ppca_dir
                output_suffix = args.output_dir
                csv_dir = ppca_dir / f"{output_suffix}_ppca_results_csv"
                fig_dir = ppca_dir / f"{output_suffix}_ppca_results_figures"
                output_file = ppca_dir / f"{output_suffix}_ppca_results.npz"
    else:
        # Infer from npz file name or path
        npz_name = npz_file.stem
        if 'tlr22' in npz_name.lower() or 'tlr22' in str(npz_file).lower():
            output_suffix = 'tlr22'
        elif 'all_tlr' in npz_name.lower() or 'all_tlr' in str(npz_file).lower():
            output_suffix = 'all_tlr'
        else:
            output_suffix = npz_name
        
        # Output directories (dataset-specific)
        csv_dir = ppca_dir / f"{output_suffix}_ppca_results_csv"
        fig_dir = ppca_dir / f"{output_suffix}_ppca_results_figures"
        output_file = ppca_dir / f"{output_suffix}_ppca_results.npz"
    
    print("="*80)
    print("Phylogenetic PCA for ESMFold Protein Structures")
    if args.output_dir and Path(args.output_dir).is_absolute():
        print(f"Output directory: {csv_dir}")
    else:
        print(f"Dataset: {output_suffix if 'output_suffix' in locals() else 'default'}")
    print("="*80)
    
    if args.use_coordinate_pca:
        print("\nWARNING: Using coordinate-based PCA (NOT RECOMMENDED)")
        print("After Procrustes alignment, coordinates are arbitrary.")
        print("Distance-based PCA is more appropriate for aligned structures.")
    else:
        print("\nUsing distance-based PCA (default, recommended for Procrustes-aligned structures)")
    
    # Load coordinates
    print("\n1. Loading ESMFold coordinates...")
    data = np.load(npz_file, allow_pickle=True)
    protein_names = list(data.keys())
    print(f"   Found {len(protein_names)} protein structures")
    
    # Load structures
    structures = []
    valid_names = []
    
    for name in protein_names:
        coords = data[name]
        if coords.shape[1] != 3:  # Skip if not 3D coordinates
            continue
        
        if coords.shape[0] > 0:  # Only add if structure has residues
            structures.append(coords)
            valid_names.append(name)
    
    print(f"   Valid structures: {len(structures)}")
    print(f"   Structure sizes: {[s.shape[0] for s in structures[:5]]}...")
    
    # Calculate average Cα-Cα distance between adjacent residues
    # This helps inform appropriate correspondence distance threshold
    adjacent_distances = []
    for struct in structures:
        if struct.shape[0] > 1:
            # Calculate distances between residue i and i+1
            for i in range(struct.shape[0] - 1):
                dist = np.linalg.norm(struct[i+1] - struct[i])
                adjacent_distances.append(dist)
    
    if len(adjacent_distances) > 0:
        avg_adjacent_dist = np.mean(adjacent_distances)
        std_adjacent_dist = np.std(adjacent_distances)
        min_adjacent_dist = np.min(adjacent_distances)
        max_adjacent_dist = np.max(adjacent_distances)
        print(f"\n   Ca-Ca distance statistics (adjacent residues):")
        print(f"   - Average: {avg_adjacent_dist:.2f} Angstrom")
        print(f"   - Std dev: {std_adjacent_dist:.2f} Angstrom")
        print(f"   - Range: {min_adjacent_dist:.2f} - {max_adjacent_dist:.2f} Angstrom")
        print(f"   - Note: Typical Ca-Ca distance for adjacent residues is ~3.8 Angstrom")
        print(f"   - Correspondence threshold ({args.correspondence_distance_threshold} Å) is {args.correspondence_distance_threshold/avg_adjacent_dist:.1f}x the average adjacent distance")
    
    # Handle variable structure sizes using hybrid sequence-structure approach
    # NOTE: In the Polly et al. (2013) paper, all specimens had the same number of 
    # landmarks (21 landmarks + 285 semi-landmarks), so no alignment was needed.
    # For protein structures with variable residue counts, we use a hybrid approach:
    # 1. Sequence alignment identifies evolutionarily homologous positions
    # 2. Structural alignment refines spatial correspondence
    # This preserves evolutionary signal while handling structural uncertainty.
    
    structure_sizes = [s.shape[0] for s in structures]
    min_size = min(structure_sizes)
    max_size = max(structure_sizes)
    print(f"\n   Structure size range: {min_size} - {max_size} residues")
    
    # Initialize distance_features_imputed variable (used if imputation is enabled)
    distance_features_imputed = None
    
    # Determine alignment mode: MSA-based with imputation (DEFAULT if alignment provided) or Structural alignment (optional)
    # NEW DEFAULT: MSA with imputation when alignment is provided
    # Structural alignment (FoldMason) is optional - use --use-structural-alignment-only
    use_structural_only = False  # DEFAULT: Use MSA if alignment provided
    if args.alignment is None or (args.alignment is not None and not Path(args.alignment).exists()):
        # No alignment provided or file not found - use structural alignment
        use_structural_only = True
        if args.alignment is not None and not Path(args.alignment).exists():
            print(f"\n   Warning: Alignment file not found: {args.alignment}")
            print(f"   Falling back to FoldMason structural alignment (default when no alignment provided).")
    elif args.use_structural_alignment_only:
        # User explicitly requested structural alignment despite providing MSA
        use_structural_only = True
        print(f"\n   Note: MSA file provided but --use-structural-alignment-only flag set.")
        print(f"   Using FoldMason structural alignment (ignoring MSA).")
    
    # Set default imputation behavior: disabled by default (complete coverage required)
    if args.allow_imputation is None:
        # Default: no imputation (requires 100% coverage) for both MSA and structural-only
        args.allow_imputation = False
    
    if use_structural_only:
        # Structural alignment: FoldMason (optional method, used when no MSA provided)
        print(f"\n2. Structural Alignment using FoldMason (OPTIONAL METHOD)")
        print(f"   Using FoldMason for progressive multiple structure alignment")
        print(f"   Method: Foldseek-based structural alphabet (3Di) -> progressive merge")
        print(f"   Benefits: Highly scalable, optimized for high-speed alignment")
        print(f"   Note: MSA-based alignment with imputation is the DEFAULT when --alignment is provided")
        
        print(f"\n   STEP 2: FoldMason Structural Alignment")
        print(f"   {'='*70}")
        
        try:
            result = robust_structural_alignment_foldmason(
                structures,
                valid_names,
                foldmason_executable=args.foldmason_executable,
                threads=args.threads
            )
            # The function now returns 5 values (distance_features added)
            structures_aligned, common_positions_mask, alignment_stats, valid_names, distance_features_imputed = result
            
            n_common = structures_aligned[0].shape[0]
            avg_structure_size = np.mean([s.shape[0] for s in structures])
            coverage_ratio = n_common / avg_structure_size if avg_structure_size > 0 else 0
            
            print(f"\n   FoldMason alignment results:")
            print(f"   - Structures aligned: {len(structures_aligned)}/{len(structures)}")
            print(f"   - Common positions: {n_common} ({coverage_ratio*100:.1f}% of average structure size)")
            print(f"   - Note: Structures are aligned by FoldMason (progressive alignment using Foldseek alphabet)")
            print(f"   - Final optimal alignment will be done by GPA (next step).")
            
            # Update structures and names
            structures = structures_aligned
            
            # Set flags
            correspondence_map = None
            reference_structure_name = None
            reference_structure_idx = None
            common_positions = None
            
            if n_common < 10:
                raise ValueError(
                    f"FoldMason alignment found too few positions ({n_common}). "
                    f"This may indicate structures are too divergent or have different domains."
                )
            
        except Exception as e:
            import traceback
            traceback.print_exc()
            print(f"   Error in FoldMason alignment: {e}")
            raise ValueError(
                "FoldMason alignment failed. "
                "Please ensure FoldMason is installed (conda install -c bioconda foldmason), "
                "or provide an MSA with --alignment for MSA-based alignment."
            )
    
    else:
        # MSA-based approach: sequence alignment (DEFAULT METHOD)
        # Determine imputation settings
        # Note: --no-imputation overrides --allow-imputation
        use_imputation = args.allow_imputation and not args.no_imputation
        optimize = use_imputation and not args.no_optimize_coverage
        
        print(f"\n2. MSA-based Alignment (DEFAULT METHOD)")
        print(f"   Using multiple sequence alignment to identify homologous positions")
        print(f"   This preserves evolutionary signal and maximizes position coverage")
        if use_imputation:
            print(f"   Method: MSA + automatic coverage optimization + imputation (OPTIONAL METHOD)")
            print(f"   This automatically finds optimal balance between positions and imputation rate")
            print(f"   Note: Use --allow-imputation to enable this method")
        else:
            print(f"   Method: MSA with strict 100%% coverage requirement (no imputation) - DEFAULT")
            print(f"   All structures must have coordinates at each position (zero imputation bias)")
        
        # Use original alignment (preserves full evolutionary signal)
        alignment_dict, _ = parse_sequence_alignment(args.alignment)
        print(f"   Loaded alignment with {len(alignment_dict)} sequences (full evolutionary MSA)")
        
        print(f"\n   STEP 2: Sequence Alignment Matching (Correspondence Identification)")
        print(f"   {'='*70}")
        print(f"   Using sequence alignment to identify homologous positions...")
        
        try:
            
            # Extract homologous positions from structures based on MSA
            # If imputation is enabled, automatically optimize coverage threshold (unless disabled)
            result = extract_homologous_positions_from_msa(
                alignment_dict, valid_names, structures, 
                msa_threshold=args.msa_threshold,
                min_coverage=args.min_coverage,
                allow_imputation=use_imputation,
                optimize_coverage=optimize
            )
            # Unpack result (always returns 4 values, but missing_masks may be None)
            structures_aligned, msa_maps, common_columns, missing_masks = result
            
            n_common = structures_aligned[0].shape[0]
            print(f"\n   Extracted {n_common} homologous positions from MSA")
            
            # If imputation is enabled, calculate distance features with imputation
            if use_imputation and missing_masks is not None:
                print(f"\n   STEP 3d: Distance-based Imputation (for incomplete coverage)")
                print(f"   {'='*70}")
                print(f"   Generating distance features and handling gaps via imputation...")
                
                # Calculate distance features for all pairs
                n_pairs = n_common * (n_common - 1) // 2
                distance_features_imputed = np.zeros((len(structures_aligned), n_pairs))
                
                # Track imputation stats
                total_elements = len(structures_aligned) * n_pairs
                imputed_count = 0
                
                for s_idx in range(len(structures_aligned)):
                    coords = structures_aligned[s_idx]
                    mask = missing_masks[s_idx]
                    
                    pair_idx = 0
                    for i in range(n_common):
                        for j in range(i+1, n_common):
                            if not mask[i] and not mask[j]:
                                # Both residues exist, calculate real distance
                                distance_features_imputed[s_idx, pair_idx] = np.linalg.norm(coords[i] - coords[j])
                            else:
                                # One or both are gaps, mark for imputation
                                distance_features_imputed[s_idx, pair_idx] = -1.0
                                imputed_count += 1
                            pair_idx += 1
                
                # Perform Mean Imputation on the distances (Polly et al. 2013 method)
                for p_idx in range(n_pairs):
                    col = distance_features_imputed[:, p_idx]
                    valid_vals = col[col >= 0]
                    if len(valid_vals) > 0:
                        mean_dist = np.mean(valid_vals)
                        # Replace placeholders
                        distance_features_imputed[col < 0, p_idx] = mean_dist
                    else:
                        distance_features_imputed[:, p_idx] = 0.0  # Should not happen with reasonable coverage
                
                imputation_rate = (imputed_count / total_elements) * 100
                print(f"   Distance Imputation complete:")
                print(f"   - Imputation Rate: {imputation_rate:.2f}% (Percentage of values substituted by mean)")
                print(f"   - Scientific guidelines:")
                print(f"     * < 5%: Excellent, results are highly reliable with negligible bias")
                print(f"     * 5-7%: Good, acceptable for most analyses")
                print(f"     * 7-10%: Moderate, consider using higher --min-coverage")
                print(f"     * > 10%: High bias risk, strongly recommend using higher --min-coverage (e.g., 0.85-0.90)")
                
                # Impute coordinates for GPA alignment (mean position of valid structures)
                for p in range(n_common):
                    valid_coords = [structures_aligned[s_idx][p] for s_idx in range(len(structures_aligned)) 
                                   if not missing_masks[s_idx][p]]
                    if valid_coords:
                        mean_pos = np.mean(valid_coords, axis=0)
                        for s_idx in range(len(structures_aligned)):
                            if missing_masks[s_idx][p]:
                                structures_aligned[s_idx][p] = mean_pos
            
            # Assess MSA alignment quality to determine if structural refinement is needed
            # Good MSA alignment indicators:
            # 1. High number of positions extracted (relative to structure sizes)
            # 2. All structures have same size (complete coverage)
            # 3. High coverage percentage
            structure_sizes_before = [s.shape[0] for s in structures]
            avg_structure_size = np.mean(structure_sizes_before)
            coverage_ratio = n_common / avg_structure_size if avg_structure_size > 0 else 0
            all_same_size = len(set([s.shape[0] for s in structures_aligned])) == 1
            
            # Auto-detect if MSA alignment is good enough to skip refinement
            msa_alignment_quality_good = (
                coverage_ratio >= 0.8 and  # At least 80% of average structure size
                all_same_size and  # All structures have same number of positions
                n_common >= 50  # At least 50 positions (reasonable minimum)
            )
            
            # Report MSA alignment quality
            print(f"\n   MSA alignment quality assessment:")
            print(f"   - Positions extracted: {n_common} ({coverage_ratio*100:.1f}% of average structure size)")
            print(f"   - All structures have same size: {all_same_size}")
            print(f"   - MSA alignment quality: {'GOOD' if msa_alignment_quality_good else 'MODERATE'}")
            
            # Step 2: Structural refinement (optional but recommended)
            # Auto-skip if MSA alignment is very good and --auto-skip-refinement is set
            if args.auto_skip_refinement and msa_alignment_quality_good:
                print(f"   - Auto-skipping structural refinement (MSA alignment quality is good)")
                should_refine = False
            else:
                should_refine = args.structural_refinement
                if msa_alignment_quality_good and args.structural_refinement:
                    print(f"   - Recommendation: Structural refinement may be skipped (MSA positions are well-aligned)")
                    print(f"   - Continuing with refinement anyway (use --no-structural-refinement or --auto-skip-refinement to skip)")
            
            if should_refine:
                print(f"\n   STEP 4: Structural Refinement (Correspondence Identification)")
                print(f"   {'='*70}")
                
                if args.use_foldmason_refinement:
                    print(f"   Using FoldMason for structural refinement (MSA + FoldMason hybrid approach)...")
                    print(f"   This allows imputation of missing positions, potentially capturing more positions.")
                    print(f"   MSA identified {n_common} candidate positions, FoldMason will verify structural correspondence.")
                    
                    try:
                        # Use FoldMason to refine MSA-selected positions
                        # This allows imputation for positions where some structures are missing
                        result = robust_structural_alignment_foldmason(
                            structures_aligned,
                            valid_names,
                            foldmason_executable=args.foldmason_executable,
                            threads=args.threads
                        )
                        structures_aligned, common_positions_mask, alignment_stats, valid_names, distance_features_imputed = result
                        n_common = structures_aligned[0].shape[0]
                        print(f"   After FoldMason refinement: {n_common} positions retained")
                        print(f"   Note: FoldMason may have found additional structurally conserved positions")
                        print(f"   or verified/imputed MSA positions with missing data.")
                    except Exception as e:
                        print(f"   Warning: FoldMason refinement failed: {e}")
                        print(f"   Falling back to Procrustes refinement...")
                        # Fall back to Procrustes
                        structures_aligned = structural_alignment_improved(
                            structures_aligned,
                            reference_idx=0,
                            distance_threshold=args.correspondence_distance_threshold,
                            preserve_positions=True
                        )
                        n_common = structures_aligned[0].shape[0]
                        print(f"   After Procrustes refinement: {n_common} positions retained")
                else:
                    print(f"   Applying Procrustes structural alignment refinement...")
                    print(f"   (This verifies spatial correspondence of homologous positions)")
                    print(f"   Note: All MSA-selected positions will be preserved (no truncation)")
                    try:
                        # Refine alignment using structural information
                        # CRITICAL: preserve_positions=True ensures we don't discard MSA-selected positions
                        structures_aligned = structural_alignment_improved(
                            structures_aligned,
                            reference_idx=0,
                            distance_threshold=args.correspondence_distance_threshold,  # Use same threshold as correspondence finding
                            preserve_positions=True  # Preserve all MSA-selected positions
                        )
                        n_common = structures_aligned[0].shape[0]
                        print(f"   After structural refinement: {n_common} positions retained (all MSA positions preserved)")
                    except Exception as e:
                        print(f"   Warning: Structural refinement failed: {e}")
                        print(f"   Continuing with sequence-based positions only...")
            else:
                if args.auto_skip_refinement and msa_alignment_quality_good:
                    print(f"\n   STEP 4: Structural Refinement (Correspondence Identification)")
                    print(f"   {'='*70}")
                    print(f"   Auto-skipping structural refinement (--auto-skip-refinement)")
                    print(f"   MSA alignment quality is good (coverage: {coverage_ratio*100:.1f}%, all structures aligned)")
                    print(f"   Using MSA-selected positions directly ({n_common} positions)")
                elif msa_alignment_quality_good:
                    print(f"\n   STEP 4: Structural Refinement (Correspondence Identification)")
                    print(f"   {'='*70}")
                    print(f"   Skipping structural refinement (--no-structural-refinement)")
                    print(f"   MSA alignment quality is good (coverage: {coverage_ratio*100:.1f}%, all structures aligned)")
                    print(f"   Using MSA-selected positions directly ({n_common} positions)")
                else:
                    print(f"\n   STEP 4: Structural Refinement (Correspondence Identification)")
                    print(f"   {'='*70}")
                    print(f"   Skipping structural refinement (--no-structural-refinement)")
                    print(f"   Note: MSA alignment quality is moderate (coverage: {coverage_ratio*100:.1f}%)")
                    print(f"   Using MSA-selected positions directly ({n_common} positions)")
            
        except Exception as e:
            print(f"   Error in hybrid approach: {e}")
            raise ValueError(
                "Hybrid sequence-structure alignment failed. "
                "Try using --use-structural-alignment-only to skip MSA, or check alignment file."
            )
    
    # Perform Generalized Procrustes Analysis
    # GPA optimally aligns structures to the mean structure
    # This is done AFTER correspondence identification (MSA or structural) to ensure:
    # 1. All structures have consistent features (same number of positions)
    # 2. Optimal alignment to consensus mean structure
    # 3. Full Procrustes (removes size differences via scaling)
    
    print("\n3. Performing Generalized Procrustes Analysis (GPA)...")
    print("   Aligning structures with consistent features to optimal mean structure...")
    aligned_structures, mean_structure = generalized_procrustes_analysis(structures_aligned, scale=True)
    print(f"   Aligned {len(aligned_structures)} structures")
    print(f"   Mean structure shape: {mean_structure.shape}")
    
    # Procrustes ANOVA validation (sanity check)
    print("\n3a. Procrustes ANOVA validation (sanity check)...")
    try:
        anova_results = procrustes_anova(structures_aligned, aligned_structures, mean_structure, n_permutations=1000)
        print(f"   Sum of squares before alignment: {anova_results['ss_before']:.4f}")
        print(f"   Sum of squares after alignment: {anova_results['ss_after']:.4f}")
        print(f"   Variation explained by GPA: {anova_results['prop_explained']*100:.2f}%")
        print(f"   F-statistic: {anova_results['f_statistic']:.4f}")
        print(f"   P-value (permutation test, n={1000}): {anova_results['p_value']:.4f}")
        
        if anova_results['p_value'] < 0.05:
            print(f"   [OK] GPA alignment is statistically significant (p < 0.05)")
        else:
            print(f"   [WARNING] GPA alignment may not be significant (p >= 0.05)")
            print(f"      This could indicate alignment issues or very similar structures")
        
        if anova_results['prop_explained'] > 0.5:
            print(f"   [OK] GPA explains >50% of shape variation (good alignment)")
        elif anova_results['prop_explained'] > 0.2:
            print(f"   [NOTE] GPA explains 20-50% of shape variation (moderate alignment)")
        else:
            print(f"   [WARNING] GPA explains <20% of shape variation (poor alignment)")
            print(f"      Note: For very similar structures (e.g., same domain), low % is expected")
    except Exception as e:
        print(f"   Warning: Procrustes ANOVA failed: {e}")
        print(f"   Continuing with analysis...")
    
    n_samples = len(aligned_structures)
    n_common_residues = aligned_structures[0].shape[0]
    
    # Choose feature type: distances (default) or coordinates (legacy)
    if args.use_coordinate_pca:
        # Coordinate-based PCA (legacy approach)
        print("\n4. Computing coordinate-based features...")
        n_features = n_common_residues * 3  # x, y, z for each residue
        X = np.array([s.flatten() for s in aligned_structures])
        X_mean = X.mean(axis=0)
        X_centered = X - X_mean
        feature_type = "coordinates"
        print(f"   Data matrix shape: {X_centered.shape}")
        print(f"   (n_samples={n_samples}, n_features={n_features})")
    else:
        # Distance-based PCA
        print("\n4. Using pre-calculated distance features (with imputation)...")
        
        # If we used FoldMason or MSA with imputation, we already have the imputed distance matrix
        if distance_features_imputed is not None:
            distance_features = distance_features_imputed
            pair_indices = None # We'll handle labels later
        else:
            # Fallback for MSA-based approach
            if args.max_distance is not None or args.max_sequence_separation is not None:
                result = compute_distance_features(
                    aligned_structures,
                    max_distance=args.max_distance,
                    sequence_separation=args.max_sequence_separation
                )
                distance_features, pair_indices = result
            else:
                distance_features = compute_distance_features(aligned_structures)
                pair_indices = None
        
        n_features = distance_features.shape[1]
        print(f"   Distance features shape: {distance_features.shape}")
        print(f"   (n_samples={n_samples}, n_pairs={n_features})")
        
        # Mean-centering distance features
        X_mean = distance_features.mean(axis=0)
        X_centered = distance_features - X_mean
        feature_type = "distances"
        print(f"   Distance features mean-centered")
    
    print(f"\n   Data matrix shape: {X_centered.shape}")
    print(f"   (n_samples={n_samples}, n_features={n_features})")
    
    # Standard PCA
    print(f"\n5. Performing standard PCA on {feature_type}...")
    pca_scores, pca_eigenvectors, pca_eigenvalues = standard_pca(X_centered)
    print(f"   PCA eigenvalues (first 10): {pca_eigenvalues[:10]}")
    print(f"   Variance explained by PC1: {pca_eigenvalues[0] / pca_eigenvalues.sum() * 100:.2f}%")
    if not args.use_coordinate_pca:
        print(f"   NOTE: PCs are based on distance features, not coordinates")
    
    # Load phylogenetic tree
    print("\n6. Loading phylogenetic tree...")
    if Phylo is None:
        print("   [ERROR] BioPython is required for phylogenetic PCA!")
        print("   Please install BioPython: pip install biopython")
        print("   or: conda install -c conda-forge biopython")
        print("\n   Without BioPython, phylogenetic covariance matrix cannot be computed.")
        print("   The analysis will fail. Please install BioPython and rerun.")
        raise ImportError(
            "BioPython is required for phylogenetic PCA. "
            "Install with: pip install biopython"
        )
    else:
        try:
            tree = Phylo.read(str(tree_file), 'newick')
            print(f"   Tree loaded with {len(list(tree.get_terminals()))} tips")
            
            # Build phylogenetic covariance matrix
            print("\n7. Building phylogenetic covariance matrix...")
            C = build_phylogenetic_covariance_matrix(tree, valid_names)
            print(f"   C matrix shape: {C.shape}")
            print(f"   C matrix diagonal range: [{C.diagonal().min():.4f}, {C.diagonal().max():.4f}]")
        except Exception as e:
            print(f"   Error loading tree: {e}")
            print("   Using placeholder covariance matrix")
            C = np.eye(n_samples) * 0.1
    
    # Phylogenetic PCA
    print(f"\n8. Performing phylogenetic PCA on {feature_type}...")
    try:
        ppca_scores, ppca_eigenvectors, ppca_eigenvalues, ancestral_root = phylogenetic_pca(X_centered, C)
        print(f"   pPCA eigenvalues (first 10): {ppca_eigenvalues[:10]}")
        print(f"   Variance explained by pPC1: {ppca_eigenvalues[0] / ppca_eigenvalues.sum() * 100:.2f}%")
        print(f"   Ancestral root shape: {ancestral_root.shape}")
        if not args.use_coordinate_pca:
            print(f"   NOTE: pPCs are based on distance features, not coordinates")
    except Exception as e:
        print(f"   Error in pPCA: {e}")
        import traceback
        traceback.print_exc()
        return
    
    # Compare PCA and pPCA
    print("\n9. Comparing PCA and pPCA results...")
    
    # Check score correlations
    if pca_scores.shape[1] > 1 and ppca_scores.shape[1] > 1:
        pca_corr = np.corrcoef(pca_scores[:, 0], pca_scores[:, 1])[0, 1]
        ppca_corr = np.corrcoef(ppca_scores[:, 0], ppca_scores[:, 1])[0, 1]
        print(f"   Correlation between PC1 and PC2: {pca_corr:.4f}")
        print(f"   Correlation between pPC1 and pPC2: {ppca_corr:.4f}")
    
    # Check variance conservation
    pca_var = pca_scores.var(axis=0).sum()
    ppca_var = ppca_scores.var(axis=0).sum()
    orig_var = X_centered.var(axis=0).sum()
    print(f"\n   Variance in original data: {orig_var:.6f}")
    print(f"   Variance in PCA scores: {pca_var:.6f}")
    print(f"   Variance in pPCA scores: {ppca_var:.6f}")
    print(f"   (Should be approximately equal)")
    
    # Save results
    print("\n10. Saving results...")
    save_dict = {
        'pca_scores': pca_scores,
        'pca_eigenvectors': pca_eigenvectors,
        'pca_eigenvalues': pca_eigenvalues,
        'ppca_scores': ppca_scores,
        'ppca_eigenvectors': ppca_eigenvectors,
        'ppca_eigenvalues': ppca_eigenvalues,
        'ancestral_root': ancestral_root,
        'mean_structure': mean_structure,
        'aligned_structures': aligned_structures,
        'protein_names': np.array(valid_names),
        'phylogenetic_covariance': C
    }
    
    # Add pair_indices if filtering was used (for distance-based features)
    if not args.use_coordinate_pca and pair_indices is not None:
        save_dict['distance_pair_indices'] = np.array(pair_indices)
    
    np.savez_compressed(output_file, **save_dict)
    print(f"   NPZ file saved to: {output_file}")
    
    # Also save as CSV files
    print("\n11. Exporting results to CSV files...")
    if not HAS_PANDAS:
        print("   Skipping CSV export (pandas not available)")
    else:
        csv_dir.mkdir(exist_ok=True)
    
        # 1. PCA scores (with protein names as index)
        pca_scores_df = pd.DataFrame(pca_scores, index=valid_names)
        pca_scores_df.columns = [f'PC{i+1}' for i in range(pca_scores.shape[1])]
        pca_scores_df.to_csv(csv_dir / "pca_scores.csv")
        print(f"   Saved: pca_scores.csv ({pca_scores.shape[0]} samples × {pca_scores.shape[1]} components)")
        
        # 2. pPCA scores (with protein names as index)
        ppca_scores_df = pd.DataFrame(ppca_scores, index=valid_names)
        ppca_scores_df.columns = [f'pPC{i+1}' for i in range(ppca_scores.shape[1])]
        ppca_scores_df.to_csv(csv_dir / "ppca_scores.csv")
        print(f"   Saved: ppca_scores.csv ({ppca_scores.shape[0]} samples × {ppca_scores.shape[1]} components)")
        
        # 3. PCA eigenvalues
        pca_eigenvalues_df = pd.DataFrame({
            'Component': [f'PC{i+1}' for i in range(len(pca_eigenvalues))],
            'Eigenvalue': pca_eigenvalues,
            'Variance_Explained': pca_eigenvalues / pca_eigenvalues.sum() * 100,
            'Cumulative_Variance': np.cumsum(pca_eigenvalues) / pca_eigenvalues.sum() * 100
        })
        pca_eigenvalues_df.to_csv(csv_dir / "pca_eigenvalues.csv", index=False)
        print(f"   Saved: pca_eigenvalues.csv")
        
        # 4. pPCA eigenvalues
        ppca_eigenvalues_df = pd.DataFrame({
            'Component': [f'pPC{i+1}' for i in range(len(ppca_eigenvalues))],
            'Eigenvalue': ppca_eigenvalues,
            'Variance_Explained': ppca_eigenvalues / ppca_eigenvalues.sum() * 100,
            'Cumulative_Variance': np.cumsum(ppca_eigenvalues) / ppca_eigenvalues.sum() * 100
        })
        ppca_eigenvalues_df.to_csv(csv_dir / "ppca_eigenvalues.csv", index=False)
        print(f"   Saved: ppca_eigenvalues.csv")
        
        # 5. PCA eigenvectors (loadings)
        if args.use_coordinate_pca:
            pca_eigenvectors_df = pd.DataFrame(pca_eigenvectors)
            pca_eigenvectors_df.columns = [f'PC{i+1}' for i in range(pca_eigenvectors.shape[1])]
            pca_eigenvectors_df.index = [f'Feature_{i+1}' for i in range(pca_eigenvectors.shape[0])]
        else:
            # Distance pairs as index - handle filtered vs unfiltered cases
            if pair_indices is not None:
                # Filtered case: use actual pair indices
                pair_labels = [f'Residue_{i+1}_Residue_{j+1}' for (i, j) in pair_indices]
            else:
                # Unfiltered case: compute all pairs
                n_residues = aligned_structures[0].shape[0]
                pair_labels = []
                for i in range(n_residues):
                    for j in range(i+1, n_residues):
                        pair_labels.append(f'Residue_{i+1}_Residue_{j+1}')
            
            pca_eigenvectors_df = pd.DataFrame(pca_eigenvectors)
            pca_eigenvectors_df.columns = [f'PC{i+1}' for i in range(pca_eigenvectors.shape[1])]
            pca_eigenvectors_df.index = pair_labels
        pca_eigenvectors_df.to_csv(csv_dir / "pca_eigenvectors.csv")
        print(f"   Saved: pca_eigenvectors.csv ({pca_eigenvectors.shape[0]} features × {pca_eigenvectors.shape[1]} components)")
        
        # 6. pPCA eigenvectors (loadings)
        if args.use_coordinate_pca:
            ppca_eigenvectors_df = pd.DataFrame(ppca_eigenvectors)
            ppca_eigenvectors_df.columns = [f'pPC{i+1}' for i in range(ppca_eigenvectors.shape[1])]
            ppca_eigenvectors_df.index = [f'Feature_{i+1}' for i in range(ppca_eigenvectors.shape[0])]
        else:
            # Distance pairs as index - reuse pair_labels from above (or recompute if needed)
            if pair_indices is not None:
                # Filtered case: use actual pair indices
                pair_labels = [f'Residue_{i+1}_Residue_{j+1}' for (i, j) in pair_indices]
            else:
                # Unfiltered case: compute all pairs
                n_residues = aligned_structures[0].shape[0]
                pair_labels = []
                for i in range(n_residues):
                    for j in range(i+1, n_residues):
                        pair_labels.append(f'Residue_{i+1}_Residue_{j+1}')
            
            ppca_eigenvectors_df = pd.DataFrame(ppca_eigenvectors)
            ppca_eigenvectors_df.columns = [f'pPC{i+1}' for i in range(ppca_eigenvectors.shape[1])]
            ppca_eigenvectors_df.index = pair_labels
        ppca_eigenvectors_df.to_csv(csv_dir / "ppca_eigenvectors.csv")
        print(f"   Saved: ppca_eigenvectors.csv ({ppca_eigenvectors.shape[0]} features × {ppca_eigenvectors.shape[1]} components)")
        
        # 7. Ancestral root structure
        if args.use_coordinate_pca:
            # Flattened coordinates
            ancestral_root_df = pd.DataFrame({
                'Coordinate': [f'Residue_{i//3+1}_{["x","y","z"][i%3]}' for i in range(len(ancestral_root))],
                'Value': ancestral_root
            })
        else:
            # Distance pairs - handle filtered vs unfiltered cases
            if pair_indices is not None:
                # Filtered case: use actual pair indices
                pair_labels = [f'Residue_{i+1}_Residue_{j+1}' for (i, j) in pair_indices]
            else:
                # Unfiltered case: compute all pairs
                n_residues = aligned_structures[0].shape[0]
                pair_labels = []
                for i in range(n_residues):
                    for j in range(i+1, n_residues):
                        pair_labels.append(f'Residue_{i+1}_Residue_{j+1}')
            
            ancestral_root_df = pd.DataFrame({
                'Distance_Pair': pair_labels,
                'Value': ancestral_root
            })
        ancestral_root_df.to_csv(csv_dir / "ancestral_root.csv", index=False)
        print(f"   Saved: ancestral_root.csv")
        
        # 8. Mean structure (consensus)
        mean_structure_df = pd.DataFrame(mean_structure, columns=['X', 'Y', 'Z'])
        mean_structure_df.index = [f'Residue_{i+1}' for i in range(mean_structure.shape[0])]
        mean_structure_df.to_csv(csv_dir / "mean_structure.csv")
        print(f"   Saved: mean_structure.csv ({mean_structure.shape[0]} residues)")
        
        # 9. Aligned structures (flattened for each structure)
        aligned_structures_array = np.array(aligned_structures)
        aligned_flat = aligned_structures_array.reshape(aligned_structures_array.shape[0], -1)
        aligned_df = pd.DataFrame(aligned_flat, index=valid_names)
        aligned_df.columns = [f'Residue_{i//3+1}_{["x","y","z"][i%3]}' 
                             for i in range(aligned_flat.shape[1])]
        aligned_df.to_csv(csv_dir / "aligned_structures.csv")
        print(f"   Saved: aligned_structures.csv ({aligned_structures_array.shape[0]} structures)")
        
        # 10. Phylogenetic covariance matrix
        C_df = pd.DataFrame(C, index=valid_names, columns=valid_names)
        C_df.to_csv(csv_dir / "phylogenetic_covariance_matrix.csv")
        print(f"   Saved: phylogenetic_covariance_matrix.csv ({C.shape[0]} × {C.shape[1]})")
    
        # 11. Protein names
        protein_names_df = pd.DataFrame({'Protein_Name': valid_names})
        protein_names_df.to_csv(csv_dir / "protein_names.csv", index=False)
        print(f"   Saved: protein_names.csv ({len(valid_names)} proteins)")
        
        # 12. Correspondence map (not applicable for structural alignment - no reference-based mapping)
        # This allows mapping PC loadings back to specific amino acids in original sequences
        # Note: Not applicable for structural alignment (no reference-based mapping)
        if 'correspondence_map' in locals() and correspondence_map is not None:
            # Create mapping table: aligned_position -> original_position for each structure
            # correspondence_map[i][ref_pos] = original residue index in structure i corresponding to reference position ref_pos
            
            # Find reference positions that were kept (have correspondences in all structures)
            n_ref = len(correspondence_map[0])
            ref_positions_kept = []
            for ref_pos in range(n_ref):
                if all(correspondence_map[i][ref_pos] is not None for i in range(len(correspondence_map))):
                    ref_positions_kept.append(ref_pos)
            
            # Build mapping table
            mapping_rows = []
            for aligned_pos, ref_pos in enumerate(ref_positions_kept):
                row = {
                    'Aligned_Position': aligned_pos + 1,  # 1-indexed
                    'Reference_Position': ref_pos + 1,    # 1-indexed (position in reference structure)
                    'Reference_Structure': reference_structure_name if 'reference_structure_name' in locals() else valid_names[0]
                }
                
                # Add original position for each structure
                for struct_idx, struct_name in enumerate(valid_names):
                    orig_idx = correspondence_map[struct_idx][ref_pos]
                    if orig_idx is not None:
                        row[f'{struct_name}_Original_Position'] = orig_idx + 1  # 1-indexed
                    else:
                        row[f'{struct_name}_Original_Position'] = None
                
                mapping_rows.append(row)
            
            if mapping_rows:
                correspondence_df = pd.DataFrame(mapping_rows)
                correspondence_df.to_csv(csv_dir / "correspondence_map.csv", index=False)
                print(f"   Saved: correspondence_map.csv ({len(mapping_rows)} aligned positions)")
                print(f"   This file maps aligned positions back to original sequence positions")
                print(f"   Use this to identify which amino acids correspond to PC loadings")
                print(f"   Example: Aligned_Position 28 -> check Original_Position columns for each structure")
        
        print(f"\n   All CSV files saved to: {csv_dir}/")
    
    # Generate figures
    print("\n12. Generating figures...")
    if not HAS_MATPLOTLIB:
        print("   Skipping figure generation (matplotlib not available)")
    else:
        fig_dir.mkdir(exist_ok=True)
        
        # 1. PCA ordination plot (PC1 vs PC2) — plot every point, verify count
        n_pca = pca_scores.shape[0]
        pca_valid = np.isfinite(pca_scores[:, 0]) & np.isfinite(pca_scores[:, 1])
        n_pca_plot = np.sum(pca_valid)
        fig, ax = plt.subplots(figsize=(10, 8))
        ax.scatter(pca_scores[pca_valid, 0], pca_scores[pca_valid, 1], alpha=0.6, s=50)
        ax.set_xlabel(f'PC1 ({pca_eigenvalues[0]/pca_eigenvalues.sum()*100:.1f}% variance)', fontsize=12)
        ax.set_ylabel(f'PC2 ({pca_eigenvalues[1]/pca_eigenvalues.sum()*100:.1f}% variance)', fontsize=12)
        ax.set_title('PCA Ordination: PC1 vs PC2', fontsize=14, fontweight='bold')
        ax.grid(True, alpha=0.3)
        ax.text(0.02, 0.98, f'N = {n_pca_plot}', transform=ax.transAxes, fontsize=11,
                verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
        plt.tight_layout()
        plt.savefig(fig_dir / "pca_ordination_pc1_pc2.png", dpi=300, bbox_inches='tight')
        plt.close()
        print(f"   Saved: pca_ordination_pc1_pc2.png ({n_pca_plot} points plotted)")
        
        # 2. pPCA ordination plot (pPC1 vs pPC2) — plot every point, verify count
        n_ppca = ppca_scores.shape[0]
        valid_mask = np.isfinite(ppca_scores[:, 0]) & np.isfinite(ppca_scores[:, 1])
        n_valid = np.sum(valid_mask)
        if n_valid < n_ppca:
            print(f"   WARNING: {n_ppca - n_valid} pPCA score(s) have NaN/Inf; plotting {n_valid} points.")
        x_plot = ppca_scores[valid_mask, 0]
        y_plot = ppca_scores[valid_mask, 1]
        fig, ax = plt.subplots(figsize=(10, 8))
        scatter = ax.scatter(x_plot, y_plot, alpha=0.6, s=50, color='orange')
        ax.set_xlabel(f'pPC1 ({ppca_eigenvalues[0]/ppca_eigenvalues.sum()*100:.1f}% variance)', fontsize=12)
        ax.set_ylabel(f'pPC2 ({ppca_eigenvalues[1]/ppca_eigenvalues.sum()*100:.1f}% variance)', fontsize=12)
        ax.set_title('Phylogenetic PCA Ordination: pPC1 vs pPC2', fontsize=14, fontweight='bold')
        ax.grid(True, alpha=0.3)
        ax.text(0.02, 0.98, f'N = {n_valid}', transform=ax.transAxes, fontsize=11,
                verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
        plt.tight_layout()
        plt.savefig(fig_dir / "ppca_ordination_ppc1_ppc2.png", dpi=300, bbox_inches='tight')
        plt.close()
        print(f"   Saved: ppca_ordination_ppc1_ppc2.png ({n_valid} points plotted)")
        
        # 3. Comparison plot: PCA vs pPCA (PC1 vs pPC1)
        fig, ax = plt.subplots(figsize=(10, 8))
        ax.scatter(pca_scores[:, 0], ppca_scores[:, 0], alpha=0.6, s=50)
        ax.set_xlabel('PC1 (Standard PCA)', fontsize=12)
        ax.set_ylabel('pPC1 (Phylogenetic PCA)', fontsize=12)
        ax.set_title('Comparison: PC1 vs pPC1', fontsize=14, fontweight='bold')
        ax.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.savefig(fig_dir / "comparison_pc1_vs_ppc1.png", dpi=300, bbox_inches='tight')
        plt.close()
        print(f"   Saved: comparison_pc1_vs_ppc1.png")
        
        # 4. Scree plot for PCA eigenvalues
        n_components_plot = min(20, len(pca_eigenvalues))
        fig, ax = plt.subplots(figsize=(10, 6))
        ax.plot(range(1, n_components_plot+1), pca_eigenvalues[:n_components_plot], 
                'o-', linewidth=2, markersize=8, label='PCA')
        ax.set_xlabel('Principal Component', fontsize=12)
        ax.set_ylabel('Eigenvalue', fontsize=12)
        ax.set_title('PCA Scree Plot', fontsize=14, fontweight='bold')
        ax.legend()
        ax.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.savefig(fig_dir / "pca_scree_plot.png", dpi=300, bbox_inches='tight')
        plt.close()
        print(f"   Saved: pca_scree_plot.png")
        
        # 5. Scree plot for pPCA eigenvalues
        n_components_plot = min(20, len(ppca_eigenvalues))
        fig, ax = plt.subplots(figsize=(10, 6))
        ax.plot(range(1, n_components_plot+1), ppca_eigenvalues[:n_components_plot], 
                'o-', linewidth=2, markersize=8, color='orange', label='pPCA')
        ax.set_xlabel('Phylogenetic Principal Component', fontsize=12)
        ax.set_ylabel('Eigenvalue', fontsize=12)
        ax.set_title('pPCA Scree Plot', fontsize=14, fontweight='bold')
        ax.legend()
        ax.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.savefig(fig_dir / "ppca_scree_plot.png", dpi=300, bbox_inches='tight')
        plt.close()
        print(f"   Saved: ppca_scree_plot.png")
        
        # 6. Variance explained comparison
        n_components_plot = min(10, len(pca_eigenvalues))
        pca_var_explained = (pca_eigenvalues / pca_eigenvalues.sum() * 100)[:n_components_plot]
        ppca_var_explained = (ppca_eigenvalues / ppca_eigenvalues.sum() * 100)[:n_components_plot]
        
        fig, ax = plt.subplots(figsize=(10, 6))
        x = np.arange(1, n_components_plot+1)
        width = 0.35
        ax.bar(x - width/2, pca_var_explained, width, label='PCA', alpha=0.8)
        ax.bar(x + width/2, ppca_var_explained, width, label='pPCA', alpha=0.8, color='orange')
        ax.set_xlabel('Component Number', fontsize=12)
        ax.set_ylabel('Variance Explained (%)', fontsize=12)
        ax.set_title('Variance Explained: PCA vs pPCA', fontsize=14, fontweight='bold')
        ax.set_xticks(x)
        ax.legend()
        ax.grid(True, alpha=0.3, axis='y')
        plt.tight_layout()
        plt.savefig(fig_dir / "variance_explained_comparison.png", dpi=300, bbox_inches='tight')
        plt.close()
        print(f"   Saved: variance_explained_comparison.png")
        
        # 7. Score correlation heatmap for pPCA (first 10 components)
        n_comp_heatmap = min(10, ppca_scores.shape[1])
        ppca_corr_matrix = np.corrcoef(ppca_scores[:, :n_comp_heatmap].T)
        
        fig, ax = plt.subplots(figsize=(10, 8))
        im = ax.imshow(ppca_corr_matrix, cmap='coolwarm', vmin=-1, vmax=1, aspect='auto')
        ax.set_xticks(range(n_comp_heatmap))
        ax.set_yticks(range(n_comp_heatmap))
        ax.set_xticklabels([f'pPC{i+1}' for i in range(n_comp_heatmap)])
        ax.set_yticklabels([f'pPC{i+1}' for i in range(n_comp_heatmap)])
        ax.set_title('pPCA Score Correlations (First 10 Components)', fontsize=14, fontweight='bold')
        cbar = plt.colorbar(im, ax=ax)
        cbar.set_label('Correlation', fontsize=11)
        
        # Add correlation values as text
        for i in range(n_comp_heatmap):
            for j in range(n_comp_heatmap):
                text = ax.text(j, i, f'{ppca_corr_matrix[i, j]:.2f}',
                             ha="center", va="center", color="black", fontsize=8)
        
        plt.tight_layout()
        plt.savefig(fig_dir / "ppca_score_correlations.png", dpi=300, bbox_inches='tight')
        plt.close()
        print(f"   Saved: ppca_score_correlations.png")
        
        # 8. PCA score correlation heatmap (should be ~0)
        n_comp_heatmap = min(10, pca_scores.shape[1])
        pca_corr_matrix = np.corrcoef(pca_scores[:, :n_comp_heatmap].T)
        
        fig, ax = plt.subplots(figsize=(10, 8))
        im = ax.imshow(pca_corr_matrix, cmap='coolwarm', vmin=-0.1, vmax=0.1, aspect='auto')
        ax.set_xticks(range(n_comp_heatmap))
        ax.set_yticks(range(n_comp_heatmap))
        ax.set_xticklabels([f'PC{i+1}' for i in range(n_comp_heatmap)])
        ax.set_yticklabels([f'PC{i+1}' for i in range(n_comp_heatmap)])
        ax.set_title('PCA Score Correlations (First 10 Components)\n(Should be ~0)', 
                    fontsize=14, fontweight='bold')
        cbar = plt.colorbar(im, ax=ax)
        cbar.set_label('Correlation', fontsize=11)
        
        # Add correlation values as text
        for i in range(n_comp_heatmap):
            for j in range(n_comp_heatmap):
                text = ax.text(j, i, f'{pca_corr_matrix[i, j]:.3f}',
                             ha="center", va="center", color="black", fontsize=8)
        
        plt.tight_layout()
        plt.savefig(fig_dir / "pca_score_correlations.png", dpi=300, bbox_inches='tight')
        plt.close()
        print(f"   Saved: pca_score_correlations.png")
        
        # 9. 3D plot: PC1 vs PC2 vs PC3
        from mpl_toolkits.mplot3d import Axes3D
        fig = plt.figure(figsize=(10, 8))
        ax = fig.add_subplot(111, projection='3d')
        ax.scatter(pca_scores[:, 0], pca_scores[:, 1], pca_scores[:, 2], 
                  alpha=0.6, s=50)
        ax.set_xlabel(f'PC1 ({pca_eigenvalues[0]/pca_eigenvalues.sum()*100:.1f}%)', fontsize=10)
        ax.set_ylabel(f'PC2 ({pca_eigenvalues[1]/pca_eigenvalues.sum()*100:.1f}%)', fontsize=10)
        ax.set_zlabel(f'PC3 ({pca_eigenvalues[2]/pca_eigenvalues.sum()*100:.1f}%)', fontsize=10)
        ax.set_title('PCA 3D Ordination: PC1 vs PC2 vs PC3', fontsize=12, fontweight='bold')
        plt.tight_layout()
        plt.savefig(fig_dir / "pca_3d_ordination.png", dpi=300, bbox_inches='tight')
        plt.close()
        print(f"   Saved: pca_3d_ordination.png")
        
        # 10. 3D plot: pPC1 vs pPC2 vs pPC3
        fig = plt.figure(figsize=(10, 8))
        ax = fig.add_subplot(111, projection='3d')
        ax.scatter(ppca_scores[:, 0], ppca_scores[:, 1], ppca_scores[:, 2], 
                  alpha=0.6, s=50, color='orange')
        ax.set_xlabel(f'pPC1 ({ppca_eigenvalues[0]/ppca_eigenvalues.sum()*100:.1f}%)', fontsize=10)
        ax.set_ylabel(f'pPC2 ({ppca_eigenvalues[1]/ppca_eigenvalues.sum()*100:.1f}%)', fontsize=10)
        ax.set_zlabel(f'pPC3 ({ppca_eigenvalues[2]/ppca_eigenvalues.sum()*100:.1f}%)', fontsize=10)
        ax.set_title('pPCA 3D Ordination: pPC1 vs pPC2 vs pPC3', fontsize=12, fontweight='bold')
        plt.tight_layout()
        plt.savefig(fig_dir / "ppca_3d_ordination.png", dpi=300, bbox_inches='tight')
        plt.close()
        print(f"   Saved: ppca_3d_ordination.png")
        
        print(f"\n   All figures saved to: {fig_dir}/")
    
    print("\n" + "="*80)
    print("Analysis complete!")
    print("="*80)
    print("\nKey findings:")
    print("- pPCA adjusts axis orientation for phylogenetic relationships")
    print("- pPCA scores are correlated (unlike PCA scores)")
    print("- Variance is conserved in both PCA and pPCA")
    print("- Use pPCA when you want to study non-phylogenetic shape variation")
    if not args.use_coordinate_pca:
        print(f"\nNOTE: This analysis used distance-based PCA (recommended for Procrustes-aligned structures)")
        print(f"      To use coordinate-based PCA, run with --use-coordinate-pca flag")

if __name__ == "__main__":
    main()

