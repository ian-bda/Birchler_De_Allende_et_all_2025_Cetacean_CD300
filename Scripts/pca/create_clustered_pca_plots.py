#!/usr/bin/env python3
"""
Create PC1 vs PC2 plots for clustered PCA results from CSV files.
Colors and shapes match phylomorphospace_colored.R exactly.
"""

import argparse
import os
import sys
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

def classify_family(seq_id):
    """
    Classify family from sequence ID (same as phylomorphospace_colored.R).
    """
    seq_id_str = str(seq_id)
    
    # Extract genus (first part before underscore or pipe)
    if "_" in seq_id_str:
        genus = seq_id_str.split("_")[0]
    elif "|" in seq_id_str:
        genus = seq_id_str.split("|")[0]
    else:
        genus = seq_id_str
    
    # Family classifications based on Artiodactyla taxonomy
    # Camelidae (Tylopoda)
    if genus in ["Camelus", "Vicugna"]:
        return "Camelidae"
    
    # Suidae (pigs)
    if genus in ["Sus", "Phacochoerus"]:
        return "Suidae"
    
    # Cervidae (deer)
    if genus in ["Cervus", "Dama", "Muntiacus", "Odocoileus", "Rangifer"]:
        return "Cervidae"
    
    # Bovidae (cattle, sheep, antelopes)
    if genus in ["Bison", "Bos", "Bubalus", "Budorcas", "Capra", 
                 "Capricornis", "Oryx", "Ovibos", "Ovis", "Pantholops"]:
        return "Bovidae"
    
    # Moschidae (musk deer)
    if genus == "Moschus":
        return "Moschidae"
    
    # Hippopotamidae (hippos)
    if genus == "Hippopotamus":
        return "Hippopotamidae"
    
    # Cetaceans (whales and dolphins)
    if genus in ["Balaenoptera", "Delphinapterus", "Delphinus", "Eschrichtius", 
                 "Eubalaena", "Globicephala", "Kogia", "Lagenorhynchus", 
                 "Lipotes", "Monodon", "Neophocaena", "Orcinus", "Phocoena", 
                 "Physeter", "Pseudorca", "Tursiops", "Sousa", "Mesoplodon"]:
        return "Cetacea"
    
    # Hominidae (humans) - though these should be excluded
    if genus == "Homo":
        return "Hominidae"
    
    # Unknown/Other
    return "Unknown"

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

def create_cluster_plot(pca_scores_file, output_file, cluster_id):
    """
    Create PC1 vs PC2 scatter plot for a cluster with colors and shapes 
    matching phylomorphospace_colored.R exactly.
    """
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
        
        # Classify species groups (for colors) and families (for shapes)
        df['species_group'] = df[id_col].apply(classify_species_group)
        df['family'] = df[id_col].apply(classify_family)
        
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
        # Colors as RGB tuples (alpha handled separately in scatter)
        group_colors_map = {
            "Cetacean": (0.27450980392156865, 0.5098039215686274, 0.7058823529411765),  # steelblue
            "Terrestrial_Artiodactyl": (1.0, 0.5490196078431373, 0.0),  # darkorange
            "Other": (0.5, 0.5, 0.5)  # gray50
        }
        alpha_value = 0.7  # Match phylomorphospace_colored.R
        
        # Shape map (same as phylomorphospace_colored.R)
        # R pch values: 21=filled circle, 22=filled square, 23=filled diamond, 
        #               24=filled triangle up, 25=filled triangle down
        # Matplotlib markers: 'o'=circle, 's'=square, 'D'=diamond, 
        #                     '^'=triangle up, 'v'=triangle down
        family_shapes = {
            "Camelidae": 'o',        # Filled circle
            "Suidae": 's',          # Filled square
            "Cervidae": 'D',        # Filled diamond
            "Bovidae": '^',         # Filled triangle up
            "Moschidae": 'v',       # Filled triangle down
            "Hippopotamidae": 'o',  # Filled circle
            "Cetacea": 'o',         # Filled circle
            "Hominidae": 'o',       # Filled circle (yellow color override, but humans excluded)
            "Unknown": 'o'          # Filled circle
        }
        
        # Create plot
        fig = plt.figure(figsize=(10, 8))
        ax = plt.gca()
        
        # Get unique families in this cluster
        unique_families = df['family'].unique()
        unique_families = [f for f in unique_families if pd.notna(f) and f != "Unknown"]
        
        # Plot each family separately with appropriate colors and shapes
        for fam in unique_families:
            family_mask = df['family'] == fam
            if family_mask.sum() == 0:
                continue
            
            # Get the species group for this family (for color)
            family_group = df.loc[family_mask, 'species_group'].iloc[0]
            
            # Special case: Hominidae gets yellow color (though should be excluded)
            if fam == "Hominidae":
                fill_color = (1.0, 1.0, 0.0)  # yellow
                use_alpha = 1.0  # no transparency for yellow
            else:
                fill_color = group_colors_map.get(family_group, group_colors_map["Other"])
                use_alpha = alpha_value
            
            # Get shape for this family
            shape_marker = family_shapes.get(fam, 'o')
            
            # Get PC1 and PC2 values for this family
            family_pc1 = df.loc[family_mask, 'PC1'].values
            family_pc2 = df.loc[family_mask, 'PC2'].values
            
            # All shapes are filled (using facecolor) with black border
            ax.scatter(family_pc1, family_pc2, 
                      marker=shape_marker,
                      facecolor=fill_color,
                      edgecolor='black',  # black border
                      s=50,
                      alpha=use_alpha,
                      linewidth=1.2,
                      label=fam)
        
        # Set labels and title
        ax.set_xlabel(f'PC1 ({pc1_var:.1f}% variance)', fontsize=12)
        ax.set_ylabel(f'PC2 ({pc2_var:.1f}% variance)', fontsize=12)
        ax.set_title(f'Group {cluster_id + 1}: PC1 vs PC2 Structural Morphospace', fontsize=14)
        ax.legend(loc='best', fontsize=9, title='Family')
        ax.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close(fig)
        
        return True
    except Exception as e:
        print(f"  Error creating plot: {e}")
        import traceback
        traceback.print_exc()
        if 'fig' in locals():
            try:
                plt.close(fig)
            except:
                pass
        return False

def main():
    parser = argparse.ArgumentParser(
        description="Create PC1 vs PC2 plots for clustered PCA results with colors and shapes matching phylomorphospace_colored.R"
    )
    
    parser.add_argument(
        "clustered_pca_dir",
        help="Directory containing clustered PCA results (e.g., ESMFold_analysis/ESMFold_output/clustered_pca/)"
    )
    
    args = parser.parse_args()
    
    if not os.path.exists(args.clustered_pca_dir):
        print(f"Error: Directory not found: {args.clustered_pca_dir}")
        sys.exit(1)
    
    print("="*60)
    print("Creating Cluster PCA Plots")
    print("="*60)
    print(f"Directory: {args.clustered_pca_dir}")
    print("Colors and shapes match phylomorphospace_colored.R")
    print()
    
    # Find all cluster PCA scores files
    created_plots = 0
    for cluster_id in range(3):
        scores_file = os.path.join(
            args.clustered_pca_dir,
            f"clustered_pca_no_PIGR_single_domain_cluster_{cluster_id}_pca_scores.csv"
        )
        
        if os.path.exists(scores_file):
            output_file = scores_file.replace('_pca_scores.csv', '_PC1_PC2_scatter.pdf')
            print(f"Creating plot for Group {cluster_id + 1} (Cluster {cluster_id})...")
            if create_cluster_plot(scores_file, output_file, cluster_id):
                created_plots += 1
                print(f"  [OK] Created: {output_file}")
            else:
                print(f"  [ERROR] Failed to create plot for Group {cluster_id + 1}")
        else:
            print(f"Group {cluster_id + 1} scores file not found: {scores_file}")
    
    print()
    print("="*60)
    print(f"Created {created_plots} plots")
    print("="*60)
    print("\nPlot details:")
    print("  - Colors: Cetacean=steelblue, Terrestrial_Artiodactyl=darkorange, Other=gray50")
    print("  - Shapes: Camelidae=circle, Suidae=square, Cervidae=diamond,")
    print("            Bovidae=triangle_up, Moschidae=triangle_down, Cetacea=circle")
    print("  - Alpha: 0.7 (70% opacity) for overlapping points")

if __name__ == "__main__":
    main()

