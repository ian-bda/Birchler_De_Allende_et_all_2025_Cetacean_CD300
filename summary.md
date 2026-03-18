# Structural PCA and Evolutionary Disparity Analysis

This repository contains a complete workflow for analyzing protein structure evolution using ESMFold-predicted structures, Principal Component Analysis (PCA), phylomorphospace projection, disparity analysis, and evolutionary model fitting. This pipeline was designed to study immunoglobulin (Ig) domain structural variation of CD300 immune family receptors across cetacean and artiodactyl species.

---
## Overview

This analysis pipeline transforms protein sequences into structural coordinates, identifies major modes of structural variation through PCA, projects evolutionary relationships onto morphological space (phylomorphospace), tests whether structural evolution follows Brownian motion, and compares evolutionary rates between groups. The workflow consists of the following steps:

1. **Structure Prediction** (ESMFold): Predict 3D structures from protein sequences
2. **Structural PCA**: Align structures, standardize coordinates, and perform dimensionality reduction
3. **Ultrametric Tree Conversion**: Convert phylogenetic tree to ultrametric (time-calibrated) for time-based analyses
4. **Phylomorphospace**: Project phylogenetic tree onto PCA space

---
## Usage Instructions

### Prerequisites

**Python packages:**
```bash
pip install numpy pandas scipy scikit-learn matplotlib seaborn biopython
```

**R packages:**
```R
install.packages(c("ape", "phytools", "dplyr"))
# Optional: geiger package for additional evolutionary models
# install.packages("geiger")
```

### Step 1: Predict Structures with ESMFold

**Using SLURM (recommended for large datasets):**
```bash
cd ESMFold_analysis/Scripts/esmfold
sbatch run_esmfold.sh
```

**Manual execution:**
```bash
python3 Scripts/esmfold/esmfold_predict.py \
    Arteodactyl_ig_hits.fasta \
    --output-dir ESMFold_output/structures \
    --batch-size 1
```

**Options:**
- `--force`: Re-run predictions even if PDB files exist
- `--batch-size N`: Process N sequences per batch (default: 1)
- `--output-dir PATH`: Directory for output files

### Step 2: Perform Structural PCA

```bash
python3 Scripts/pca/structural_pca.py \
    ESMFold_output/structures/esmfold_coordinates_no_PIGR_single_domain.npz \
    ESMFold_output/structures/coordinate_metadata_no_PIGR_single_domain.csv \
    --output-prefix ESMFold_output/pca_results_no_PIGR_single_domain \
    --n-components 10
```

**Options:**
- `--n-components N`: Number of principal components to compute (default: 10)
- `--min-length N`: Filter structures shorter than N residues
- `--max-length N`: Pad/truncate all structures to N residues
- `--output-dir PATH`: Directory for output files

### Step 3: Make Tree and Convert Tree to Ultrametric (Time-Calibrated)

**Important**: For time-based evolutionary analyses (phylomorphospace, evolutionary models, phenograms, etc.), the tree must be ultrametric (all tips equidistant from root). IQ-TREE produces phylograms (branch lengths = substitutions), which need to be converted to chronograms (branch lengths = time).

**Famsa (for alignment):**
```bash
sbatch Tree/Scripts/famsa.sh
```

**IQtree:**
```bash
sbatch Tree/Scripts/iqtree.sh
```

**Make Ultrametric:**
**Using SLURM (recommended):**
```bash
sbatch Tree/Scripts/make_ultrametric_tree.sh
```

**Manual execution:**
```bash
Rscript Tree/scripts/make_ultrametric_tree.R \
    Tree/Tree_for_Analysis/Arteodactyl_ig_hits_no_PIGR_single_domain_iqtree.treefile \
    Tree/Tree_for_Analysis/Arteodactyl_ig_hits_no_PIGR_single_domain_iqtree_ultrametric.nwk
```

**What it does:**
- Uses a hierarchical approach to convert tree to ultrametric:
  1. **`chronos()` with relaxed clock** (most accurate, maximum likelihood dating)
  2. **`chronos()` with strict clock** (if relaxed fails)
  3. **`force.ultrametric(method="nnls")`** (preserves relative branch lengths)
  4. **`force.ultrametric(method="extend")`** (last resort)
- Scales tree to desired root height (default: 1.0)
- Produces a chronogram where branch lengths represent time

**Why ultrametric trees are needed:**
- **Phylograms** (from IQ-TREE): Branch lengths = substitutions per site (evolutionary distance)
- **Chronograms** (ultrametric): Branch lengths = time (all tips equidistant from root)
- Time-based models (OU, EB) require chronograms to model rate variation through time
- Phenograms and time-based visualizations require chronograms
- Early Burst model explicitly models rate decay over time, which requires time-calibrated branches

**Output:**
- `{tree_name}_ultrametric.nwk`: Ultrametric tree in Newick format

### Step 4: Create Phylomorphospace Plot

**Colored by species group (recommended, using ultrametric tree):**
```bash
Rscript Scripts/pca/phylomorphospace_colored.R \
    ESMFold_output/pca_results_no_PIGR_single_domain_new_tree_pca_scores.csv \
    Tree/Tree_for_Analysis/Arteodactyl_ig_hits_no_PIGR_single_domain_iqtree_ultrametric.nwk \
    ESMFold_output/phylomorphospace_colored_no_PIGR_single_domain_ultrametric \
    PC1 PC2
```

**Basic phylomorphospace (no colors):**
```bash
Rscript Scripts/pca/phylomorphospace_from_pca.R \
    ESMFold_output/pca_results_no_PIGR_single_domain_pca_scores.csv \
    Tree/Arteodactyl_ig_hits_no_PIGR_single_domain_iqtree.treefile \
    ESMFold_output/phylomorphospace_no_PIGR_single_domain \
    PC1 PC2
```

**Arguments:**
- `pca_scores_csv`: Path to PCA scores CSV
- `tree_file`: Path to Newick tree file (must match sequences in PCA data)
- `output_prefix`: Prefix for output files
- `PC_x` (optional): PC for X-axis (default: PC1)
- `PC_y` (optional): PC for Y-axis (default: PC2)


## Citation and References

If you use this pipeline, please cite:
- **ESMFold**: Lin et al. (2023) "Evolutionary-scale prediction of atomic-level protein structure with a language model"
- **IQ-TREE**: Minh et al. (2020) "IQ-TREE 2: New models and efficient methods for phylogenetic inference"
- **phytools**: Revell (2012) "phytools: an R package for phylogenetic comparative biology"

---
