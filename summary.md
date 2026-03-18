# Structural PCA and Evolutionary Disparity Analysis

This repository contains a complete workflow for analyzing protein structure evolution using ESMFold-predicted structures, Principal Component Analysis (PCA), phylomorphospace projection, disparity analysis, and evolutionary model fitting. This pipeline was designed to study immunoglobulin (Ig) domain structural variation of CD300 immune family receptors across cetacean and artiodactyl species.

## Overview

This analysis pipeline transforms protein sequences into structural coordinates, identifies major modes of structural variation through PCA, projects evolutionary relationships onto morphological space (phylomorphospace), tests whether structural evolution follows Brownian motion, and compares evolutionary rates between groups. The workflow consists of the following steps:

1. **Structure Prediction (ESMFold):** Predict 3D structures from protein sequences
2. **Structural pPCA:** Align structures, standardize coordinates, and perform phylogenetic PCA
3. **Ultrametric Tree Conversion:** Convert phylogenetic tree to ultrametric (time-calibrated) for time-based analyses
4. **Phylomorphospace:** Project phylogenetic tree onto PCA space

## Usage Instructions

### Prerequisites

Python packages:

```
pip install numpy pandas scipy scikit-learn matplotlib seaborn biopython
```

R packages:

```r
install.packages(c("ape", "phytools", "dplyr"))
```

---

### Step 1: Predict Structures with ESMFold

```bash
python3 ESMFold_analysis/Scripts/esmfold/esmfold_predict.py \
    <input_sequences.fasta> \
    --output-dir ESMFold_analysis/ESMFold_output/structures \
    --batch-size 1
```

Options:

- `--force`: Re-run predictions even if PDB files exist
- `--batch-size N`: Process N sequences per batch (default: 1)
- `--output-dir PATH`: Directory for output PDB files

---

### Step 2: Perform Structural pPCA

This step performs **phylogenetic PCA (pPCA)** on the extracted structural coordinates. Unlike standard PCA, pPCA accounts for the non-independence of species observations due to shared evolutionary history by computing principal components from a phylogenetically corrected variance-covariance matrix (Revell 2009). This ensures that the resulting PC axes reflect axes of independent evolutionary variation rather than axes inflated by phylogenetic signal, making downstream model fitting and disparity comparisons statistically valid.

```bash
python3 ESMFold_analysis/Scripts/ppca/ppca_protein_structures.py \
    ESMFold_analysis/ESMFold_output/structures/esmfold_coordinates_no_PIGR_single_domain.npz \
    ESMFold_analysis/ESMFold_output/structures/coordinate_metadata_no_PIGR_single_domain.csv \
    --output-prefix ESMFold_analysis/ESMFold_output/pca_results_no_PIGR_single_domain
```

Options:

- `--n-components N`: Number of principal components to compute
- `--min-length N`: Filter structures shorter than N residues
- `--max-length N`: Pad/truncate all structures to N residues
- `--output-dir PATH`: Directory for output files

---

### Step 3: Build Tree and Convert to Ultrametric

**Important:** For time-based evolutionary analyses (phylomorphospace, evolutionary models, phenograms, etc.), the tree must be ultrametric (all tips equidistant from root). IQ-TREE produces phylograms (branch lengths = substitutions), which need to be converted to chronograms (branch lengths = time).

**Alignment (FAMSA):**

```bash
sbatch Tree/Scripts/famsa.sh
```

**Tree inference (IQ-TREE):**

```bash
sbatch Tree/Scripts/iqtree.sh
```

**Ultrametric conversion — using SLURM (recommended):**

```bash
sbatch Tree/Scripts/make_ultrametric_tree.sh
```

**Ultrametric conversion — manual execution:**

```bash
Rscript Tree/Scripts/make_ultrametric_tree.R \
    Tree/Tree_for_Analysis/Arteodactyl_ig_hits_no_PIGR_single_domain_iqtree.treefile \
    Tree/Tree_for_Analysis/Arteodactyl_ig_hits_no_PIGR_single_domain_iqtree_ultrametric.nwk
```

**What it does:**

Uses a hierarchical approach to convert the tree to ultrametric:
1. `chronos()` with relaxed clock (most accurate, maximum likelihood dating)
2. `chronos()` with strict clock (if relaxed fails)
3. `force.ultrametric(method="nnls")` (preserves relative branch lengths)
4. `force.ultrametric(method="extend")` (last resort)

Scales tree to desired root height (default: 1.0) and produces a chronogram where branch lengths represent time.

**Why ultrametric trees are needed:**
- Phylograms (from IQ-TREE): Branch lengths = substitutions per site
- Chronograms (ultrametric): Branch lengths = time (all tips equidistant from root)
- Time-based models (OU, EB) and phenograms require chronograms

**Output:** `Tree/Tree_for_Analysis/Arteodactyl_ig_hits_no_PIGR_single_domain_iqtree_ultrametric.nwk`

---

### Step 4: Create Phylomorphospace Plot

**Colored by species group (2D, recommended):**

```bash
Rscript ESMFold_analysis/Scripts/ppca/phylomorphospace_colored.R \
    ESMFold_analysis/ESMFold_output/pca_results_no_PIGR_single_domain_pca_scores.csv \
    Tree/Tree_for_Analysis/Arteodactyl_ig_hits_no_PIGR_single_domain_iqtree_ultrametric.nwk \
    ESMFold_analysis/ESMFold_output/phylomorphospace_colored_no_PIGR_single_domain_ultrametric \
    PC1 PC2
```

**3D phylomorphospace:**

```bash
Rscript ESMFold_analysis/Scripts/ppca/create_3d_phylomorphospace_universal.R \
    ESMFold_analysis/ESMFold_output/pca_results_no_PIGR_single_domain_pca_scores.csv \
    Tree/Tree_for_Analysis/Arteodactyl_ig_hits_no_PIGR_single_domain_iqtree_ultrametric.nwk \
    ESMFold_analysis/ESMFold_output/phylomorphospace_3d_no_PIGR_single_domain \
    PC1 PC2 PC3
```

**Arguments:**

- `pca_scores_csv`: Path to PCA scores CSV output from Step 2
- `tree_file`: Path to ultrametric Newick tree file
- `output_prefix`: Prefix for output files
- `PC_x`, `PC_y` (optional): PCs for axes (default: PC1, PC2)

---

## Citation and References

If you use this pipeline, please cite:

- **ESMFold:** Lin et al. (2023) "Evolutionary-scale prediction of atomic-level protein structure with a language model"
- **IQ-TREE:** Minh et al. (2020) "IQ-TREE 2: New models and efficient methods for phylogenetic inference"
- **phytools:** Revell (2012) "phytools: an R package for phylogenetic comparative biology (and other things)"
- **Phylogenetic PCA:** Revell (2009) "Size-correction and principal components for interspecific comparative studies." *Evolution*, 63, 3258–3268
- **FAMSA:** Deorowicz et al. (2016) "FAMSA: Fast and accurate multiple sequence alignment of huge protein families"
