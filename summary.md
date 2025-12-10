# Structural PCA and Evolutionary Disparity Analysis

This repository contains a complete workflow for analyzing protein structure evolution using ESMFold-predicted structures, Principal Component Analysis (PCA), phylomorphospace projection, disparity analysis, and evolutionary model fitting. This pipeline was designed to study immunoglobulin (Ig) domain structural variation of CD300 immune family receptors across cetacean and artiodactyl species.

---
## Overview

This analysis pipeline transforms protein sequences into structural coordinates, identifies major modes of structural variation through PCA, projects evolutionary relationships onto morphological space (phylomorphospace), tests whether structural evolution follows Brownian motion, and compares evolutionary rates between groups. The workflow consists of the following steps:

1. **Structure Prediction** (ESMFold): Predict 3D structures from protein sequences
2. **Structural PCA**: Align structures, standardize coordinates, and perform dimensionality reduction
3. **Ultrametric Tree Conversion**: Convert phylogenetic tree to ultrametric (time-calibrated) for time-based analyses
4. **Phylomorphospace**: Project phylogenetic tree onto PCA space
5. **Disparity Analysis**: Test the relationship between morphological and phylogenetic distances
6. **Evolutionary Models**: Fit evolutionary models (BM, OU, EB) and test for rate differences between groups
7. **Convex Hull Analysis**: Calculate and compare morphological space occupied by cetaceans vs non-cetaceans
8. **EB Simulation Analysis**: Test if observed morphological patterns deviate from Early Burst model predictions
9. **Phenograms**: Visualize trait evolution over time along the phylogeny
10. **BAMM Correlation Analysis**: Test correlation between diversification rates and shape evolution rates

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

**Note**: For this analysis, use the filtered dataset excluding PIGR sequences and multi-domain proteins:

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

### Step 3: Convert Tree to Ultrametric (Time-Calibrated)

**Important**: For time-based evolutionary analyses (phylomorphospace, evolutionary models, phenograms, etc.), the tree must be ultrametric (all tips equidistant from root). IQ-TREE produces phylograms (branch lengths = substitutions), which need to be converted to chronograms (branch lengths = time).

**Using SLURM (recommended):**
```bash
sbatch Tree/Scripts/make_ultrametric_tree.sh
```

**Manual execution:**
```bash
Rscript BiFrost/scripts/make_ultrametric_tree.R \
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

**Note**: The colored version automatically excludes humans and colors points by family (Cetacean, Terrestrial Artiodactyl) with different shapes for each family.

### Step 5: Perform Disparity Analysis

**Using ultrametric tree (recommended):**
```bash
sbatch ESMFold_analysis/Scripts/pca/run_disparity_analysis.sh
```

**Manual execution:**
```bash
python3 Scripts/pca/disparity_analysis.py \
    ESMFold_output/pca_results_no_PIGR_single_domain_new_tree_pca_scores.csv \
    Tree/Tree_for_Analysis/Arteodactyl_ig_hits_no_PIGR_single_domain_iqtree_ultrametric.nwk \
    --output-prefix ESMFold_output/disparity_no_PIGR_single_domain_ultrametric \
    --exclude-humans
```

**Options:**
- `--output-prefix PREFIX`: Prefix for output files (default: "disparity")
- `--pc-components N M ...`: Which PC components to use for morphological distance (default: all)
- `--exclude-humans`: Exclude human (Homo) sequences from analysis (recommended for consistency)

**Output:**
- `{prefix}_pairwise_distances.csv`: All pairwise distances
- `{prefix}_disparity_summary.csv`: Summary statistics

### Step 6: Fit Evolutionary Models and Test Rate Differences

**Manual execution:**
```bash
Rscript Scripts/pca/evolutionary_models.R \
    ESMFold_output/pca_results_no_PIGR_single_domain_new_tree_pca_scores.csv \
    Tree/Tree_for_Analysis/Arteodactyl_ig_hits_no_PIGR_single_domain_iqtree_ultrametric.nwk \
    ESMFold_output/evolutionary_models_no_PIGR_single_domain_ultrametric \
    PC1,PC2
```

**Arguments:**
- `pca_scores_csv`: Path to PCA scores CSV
- `tree_file`: Path to Newick tree file (must match sequences in PCA data)
- `output_prefix`: Prefix for output files
- `PC_components`: Comma-separated list of PCs to analyze (default: PC1,PC2)

**What it does:**
- Fits evolutionary models (BM, OU, EB) and compares their fit using AIC
- Calculates Blomberg's K and Pagel's λ to measure phylogenetic signal
- Compares evolutionary rates between cetaceans and terrestrial artiodactyls
- Tests if different groups evolve at different rates
- Automatically excludes humans from analysis

**Output:**
- `{prefix}_no_humans_model_comparison.csv`: Model AIC comparisons
- `{prefix}_no_humans_rate_comparison.csv`: Rate differences between groups (cetaceans vs terrestrial)
- `{prefix}_no_humans_phylogenetic_signal.csv`: Blomberg's K and Pagel's λ statistics

### Step 7: Calculate Convex Hull Analysis

**Question**: How much morphological space do cetaceans vs non-cetaceans occupy, and do they overlap?

**Using SLURM:**
```bash
sbatch ESMFold_analysis/Scripts/pca/run_convex_hull_analysis.sh
```

**Manual execution:**
```bash
Rscript Scripts/pca/convex_hull_analysis.R \
    ESMFold_output/pca_results_no_PIGR_single_domain_new_tree_pca_scores.csv \
    Tree/Tree_for_Analysis/Arteodactyl_ig_hits_no_PIGR_single_domain_iqtree_ultrametric.nwk \
    ESMFold_output/convex_hull_no_PIGR_single_domain \
    PC1 PC2
```

**What it does:**
- Calculates convex hull (smallest convex polygon containing all points) for cetaceans and non-cetaceans in PCA space
- Calculates hull area for each group
- Calculates overlap area between groups
- Tests if groups occupy significantly different morphological spaces
- Automatically excludes humans

**Output:**
- `{prefix}_convex_hull_results.csv`: Hull areas, overlap, and statistics
- `{prefix}_convex_hull_plot.pdf`: Visualization of convex hulls

**Interpretation:**
- **Large hull area**: Group occupies more morphological space (higher diversity)
- **Small hull area**: Group occupies less morphological space (lower diversity)
- **High overlap**: Groups share similar morphologies
- **Low overlap**: Groups occupy distinct morphological regions

### Step 8: EB Simulation Analysis

**Question**: Does the observed morphological pattern deviate from Early Burst (EB) model predictions?

**Using SLURM:**
```bash
sbatch ESMFold_analysis/Scripts/pca/run_eb_simulation_analysis.sh
```

**Manual execution:**
```bash
Rscript Scripts/pca/eb_simulation_analysis.R \
    ESMFold_output/pca_results_no_PIGR_single_domain_new_tree_pca_scores.csv \
    Tree/Tree_for_Analysis/Arteodactyl_ig_hits_no_PIGR_single_domain_iqtree_ultrametric.nwk \
    ESMFold_output/eb_simulation_no_PIGR_single_domain \
    1000 \
    PC1 PC2
```

**What it does:**
1. Fits EB models to PC1 and PC2 separately
2. Simulates PC1 and PC2 under EB model 1000 times
3. For each simulation, calculates convex hull metrics (cetacean area, non-cetacean area, overlap, area ratio)
4. Compares observed metrics to simulated distribution
5. Calculates p-values to test if observed data falls outside EB model predictions

**Output:**
- `{prefix}_eb_simulation_results.csv`: Observed vs simulated metrics with p-values
- `{prefix}_eb_simulation_plots.pdf`: Distribution plots comparing observed to simulated

**Interpretation:**
- **p < 0.05**: Observed pattern significantly deviates from EB model (other processes likely important)
- **p ≥ 0.05**: Observed pattern consistent with EB model (early burst explains the data)
- **Observed < simulated**: Observed metric is smaller than predicted (e.g., less morphological diversity than expected)
- **Observed > simulated**: Observed metric is larger than predicted (e.g., more morphological diversity than expected)

  ---
## Workflow Explained

### ESMFold Structure Prediction

**Script**: `Scripts/esmfold/esmfold_predict.py`  
**Input**: FASTA file with protein sequences  
**Output**: PDB files (atomic coordinates) + NPZ file (C-alpha coordinates)

**What it does:**
- Loads protein sequences from a FASTA file
- Uses ESMFold (evolutionary scale modeling) to predict 3D protein structures
- Saves structures as PDB files (standard format for atomic coordinates)
- Extracts C-alpha (CA) atom coordinates from each structure
- Compiles all CA coordinates into a single compressed NumPy file (`esmfold_coordinates.npz`)
- Creates a metadata CSV mapping sequence IDs to coordinate arrays

**Why C-alpha atoms?**
- C-alpha atoms form the protein backbone
- They represent the overall fold structure consistently
- Using all atoms would introduce noise from side-chain variation
- C-alpha coordinates are sufficient for structural comparison

**Resume capability:**
- Skips sequences that already have PDB files (unless `--force` is used)
- Can extract coordinates from existing PDB files if NPZ is missing
- Ensures efficient processing of large datasets

### Standardization (Z-score Normalization)

**Why standardization is essential:**
Before PCA, all coordinate dimensions must have equal weight. Without standardization:
- Coordinates with larger numerical ranges (e.g., structures with larger absolute positions) would dominate PCA
- The first principal component would simply reflect the dimension with the largest variance, not the most biologically meaningful variation
- PCA results would be biased by arbitrary coordinate values

**How z-score normalization works:**
For each coordinate dimension (column) in the matrix:
1. **Calculate mean**: `μ = mean(all values in column)`
2. **Calculate standard deviation**: `σ = std(all values in column)`
3. **Transform each value**: `z = (x - μ) / σ`

**Result:**
- Each column has mean = 0
- Each column has standard deviation = 1
- All dimensions contribute equally to PCA
- Variance captured by PCA reflects true structural variation, not scale differences

**Mathematical formulation:**
```
X_scaled = (X - μ) / σ
```
Where:
- `X` is the original coordinate matrix
- `μ` is a vector of column means
- `σ` is a vector of column standard deviations
- Division is element-wise (broadcasting)

**After standardization:**
- Mean of each column ≈ 0 (centered)
- Standard deviation of each column = 1 (unit variance)
- All coordinate dimensions are on the same scale

---

### Disparity Analysis

**Script**: `Scripts/pca/disparity_analysis.py`  
**Question**: Does structural evolution follow Brownian motion?

**What is Brownian Motion?**
Brownian motion (also called a random walk) is a null model for trait evolution:
- Traits evolve randomly with no directional trend
- Expected trait distance increases proportionally with phylogenetic distance
- More distantly related species should be more different (on average)
- No selection, no constraint, just random drift

**Mathematical model:**
Under Brownian motion:
```
Morphological Distance = α × Phylogenetic Distance + ε
```
Where:
- `α` (slope) = rate of morphological divergence
- `ε` = random error (deviations from perfect linear relationship)

**What the analysis does:**
1. **Calculates morphological distances**: Euclidean distance in PCA space between all pairs of sequences
   - Uses specified PC components (default: all PCs)
   - Formula: `d = √(Σ(PC_i_diff)²)` for each PC
   
2. **Calculates phylogenetic distances**: Branch length distance between all pairs in the tree
   - Sum of branch lengths along the shortest path between two tips
   - Represents evolutionary time/divergence

3. **Matches distances**: Pairs sequences that appear in both PCA data and tree

4. **Tests Brownian motion**:
   - **Linear regression**: `morph_dist ~ slope × phylo_dist + intercept`
   - **Pearson correlation**: Tests linear relationship
   - **Spearman correlation**: Tests monotonic relationship (robust to non-linearity)
   - **R-squared**: Proportion of variance in morphological distance explained by phylogenetic distance

**Output statistics:**
- **Slope (α)**: Rate of morphological divergence per unit of phylogenetic distance
- **Intercept**: Expected morphological distance at zero phylogenetic distance (should be ~0)
- **R²**: Proportion of variance explained (0-1, where 1 = perfect linear relationship)
- **Pearson r**: Linear correlation coefficient (-1 to 1)
- **Pearson p**: Statistical significance of correlation
- **Spearman ρ**: Rank correlation coefficient (monotonic relationship)
- **Spearman p**: Statistical significance of rank correlation

**Interpretation:**
- **R² > 0.5**: Strong relationship → Brownian motion explains >50% of morphological variation
- **R² 0.3-0.5**: Moderate relationship → Brownian motion explains ~30-50% of variation
- **R² < 0.3**: Weak relationship → Brownian motion explains <30% of variation
  - Other factors (selection, constraints, drift) likely more important
  
- **p < 0.05**: Significant relationship (phylogenetic distance predicts morphological distance)
- **p ≥ 0.05**: No significant relationship (morphology not proportional to phylogeny)

**What weak R² means:**
- Selection may be acting on structure (directional or stabilizing)
- Structural constraints may limit variation
- Neutral drift alone cannot explain most morphological differences
- Different lineages may evolve at different rates

---

### Evolutionary Model Fitting and Rate Comparison

**Script**: `Scripts/pca/evolutionary_models.R`  
**Question**: Do cetaceans and terrestrial artiodactyls evolve at different rates?

**What it does:**
1. **Fits single-rate Brownian Motion (BM)**: Assumes all lineages evolve at the same rate
   - Estimates evolutionary rate (σ²) for each PC dimension
   - Calculates AIC for model comparison

2. **Fits multi-rate Brownian Motion**: Tests if different groups have different rates
   - Fits BM separately to cetaceans and terrestrial artiodactyls
   - Calculates rate ratio (cetacean rate / terrestrial rate)
   - Compares AIC to single-rate model

3. **Calculates phylogenetic signal statistics**: Measures how much trait variation is explained by phylogeny
   - Blomberg's K: Compares observed phylogenetic signal to Brownian motion expectation
   - Pagel's λ: Measures phylogenetic signal on a 0-1 scale
   - Both include significance tests (p-values)

4. **Compares trait means and variances**: Descriptive statistics for each group
   - Mean trait values (position in PCA space)
   - Variance in trait values (spread in PCA space)
   - Difference in means between groups

**Evolutionary Models:**
- **Brownian Motion (BM)**: Traits evolve as a random walk
  - Single rate: All lineages evolve at rate σ²
  - Multi-rate: Different groups have different rates (σ²_cetacean vs σ²_terrestrial)
- **Ornstein-Uhlenbeck (OU)**: Traits evolve toward an optimum (requires geiger package)
- **Early Burst (EB)**: Evolution is faster early in the tree (requires geiger package)

**Model Comparison:**
- Uses AIC (Akaike Information Criterion) to compare models
- Lower AIC = better model fit (accounts for model complexity)
- AIC weights show relative support for each model

**Rate Ratio Interpretation:**
- **Rate ratio = 1.0**: Groups evolve at the same rate
- **Rate ratio > 1.0**: Cetaceans evolve faster (e.g., 1.3 = 30% faster)
- **Rate ratio < 1.0**: Cetaceans evolve slower (e.g., 0.7 = 30% slower)

**Output statistics:**
- **σ² (sigma-squared)**: Evolutionary rate parameter
  - Higher σ² = faster evolution (more change per unit time)
  - Lower σ² = slower evolution (less change per unit time)
- **Rate ratio**: Cetacean rate / Terrestrial rate
- **AIC**: Model fit statistic (lower is better)
- **Mean difference**: Difference in trait means between groups

**Phylogenetic Signal Statistics:**

The script calculates two measures of phylogenetic signal to assess how much trait variation is explained by phylogeny:

1. **Blomberg's K**:
   - Measures the amount of phylogenetic signal relative to Brownian motion
   - **K = 1.0**: Traits evolve exactly as expected under Brownian motion
   - **K < 1.0**: Less phylogenetic signal than Brownian motion (traits are more similar among close relatives than expected, or more similar among distant relatives)
   - **K > 1.0**: More phylogenetic signal than Brownian motion (traits are more different among close relatives than expected)
   - **P-value < 0.05**: Significant phylogenetic signal (traits are not randomly distributed on the tree)
   - **P-value ≥ 0.05**: No significant phylogenetic signal (traits appear randomly distributed)

2. **Pagel's λ**:
   - Measures phylogenetic signal on a scale from 0 to 1
   - **λ = 0**: No phylogenetic signal (traits are independent of phylogeny, evolution is completely convergent/random)
   - **λ = 1**: Full phylogenetic signal (traits follow Brownian motion, closely related species are similar)
   - **0 < λ < 1**: Intermediate phylogenetic signal (some phylogenetic structure, but less than Brownian motion)
   - **P-value < 0.05**: Significant phylogenetic signal
   - **P-value ≥ 0.05**: No significant phylogenetic signal

**Interpreting Phylogenetic Signal:**
- **High signal (K ≈ 1, λ ≈ 1)**: Traits are strongly conserved phylogenetically. Closely related species are similar, suggesting traits are inherited from common ancestors.
- **Low signal (K < 1, λ ≈ 0)**: Traits show little phylogenetic structure. Species are similar/different regardless of relatedness, suggesting convergent evolution, strong selection, or rapid evolution.
- **Intermediate signal (0 < λ < 1)**: Some phylogenetic structure exists, but traits have evolved more independently than expected under pure Brownian motion.

### EB Simulation Analysis

**Script**: `Scripts/pca/eb_simulation_analysis.R`  
**Question**: Does observed morphology deviate from Early Burst model predictions?

**What is the Early Burst (EB) model?**
- Models adaptive radiation: Evolution is faster early in the tree and decelerates over time
- Rate formula: `rate(t) = σ² × exp(a × t)`
  - `a < 0`: Rate decreases over time (early burst)
  - `a = 0`: Rate constant (reduces to Brownian Motion)
  - `a > 0`: Rate increases over time (late burst, rare)
- Captures pattern of rapid diversification followed by stabilization

**What the analysis does:**
1. **Fits EB models**: Estimates EB parameters (σ², a) for PC1 and PC2 separately
2. **Simulates traits**: Simulates PC1 and PC2 under fitted EB model 1000 times
3. **Calculates metrics**: For each simulation, calculates:
   - Cetacean convex hull area
   - Non-cetacean convex hull area
   - Overlap area
   - Area ratio (cetacean / non-cetacean)
4. **Compares to observed**: Compares observed metrics to simulated distribution
5. **Calculates p-values**: Tests if observed data falls outside EB model predictions

**Statistical testing:**
- **P-value < 0.05**: Observed pattern significantly deviates from EB model
  - Other evolutionary processes (selection, constraints) likely important
  - EB model alone cannot explain the data
- **P-value ≥ 0.05**: Observed pattern consistent with EB model
  - Early burst pattern explains the data
  - No evidence for additional processes

**Interpretation of results:**
- **Observed < simulated (p < 0.05)**: Observed metric is smaller than predicted
  - Example: Cetacean hull area smaller than expected under EB
  - Suggests additional constraints or processes beyond early burst
- **Observed > simulated (p < 0.05)**: Observed metric is larger than predicted
  - Example: Cetacean hull area larger than expected under EB
  - Suggests additional diversification or relaxed constraints

## Citation and References

If you use this pipeline, please cite:
- **ESMFold**: Lin et al. (2023) "Evolutionary-scale prediction of atomic-level protein structure with a language model"
- **IQ-TREE**: Minh et al. (2020) "IQ-TREE 2: New models and efficient methods for phylogenetic inference"
- **phytools**: Revell (2012) "phytools: an R package for phylogenetic comparative biology"

---

## Troubleshooting

**ESMFold prediction fails:**
- Check GPU availability (ESMFold requires GPU for fast inference)
- Reduce batch size if memory issues occur
- Check FASTA file format (valid sequences required)

**PCA alignment fails:**
- Ensure all structures have valid coordinates
- Check for structures with zero or very few residues
- Verify NPZ file is not corrupted

**Phylomorphospace plot is blank:**
- Check that tree tip labels match PCA score sequence IDs
- Verify tree is rooted and dichotomous
- Check for NaN values in PCA scores

**Disparity analysis has low R²:**
- This is a biological result, not an error
- Low R² indicates weak Brownian motion (selection/constraints may be important)
- Check that tree and PCA data contain the same sequences

**Evolutionary models analysis fails:**
- Check that tree and PCA data contain the same sequences
- Ensure tree is rooted and has branch lengths
- **Important**: Use ultrametric tree for OU and EB models (required for time-based models)
- Verify that both cetaceans and terrestrial artiodactyls are present in the dataset
- If geiger package is not available, only Brownian Motion models will be fitted

**Ultrametric tree conversion fails:**
- Check that input tree is rooted and has branch lengths
- If `chronos()` fails, try using `force.ultrametric(method="nnls")` directly
- Zero-length branches may cause issues; script should handle these automatically
- Ensure tree file is in valid Newick format

**EB simulation analysis fails:**
- Requires `geiger` package: `install.packages("geiger")`
- Ensure ultrametric tree is used (required for EB model)
- Check that both cetaceans and non-cetaceans are present in the dataset
- Large number of simulations (1000) may take several hours; reduce if needed
