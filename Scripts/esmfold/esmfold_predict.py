#!/usr/bin/env python3
"""
ESMFold Protein Structure Prediction Script

This script takes a FASTA file as input and generates PDB files for each sequence
using ESMFold protein structure prediction.

Usage:
    python esmfold_predict.py input.fasta [options]

Requirements:
    - torch
    - biotite
    - esm
"""

import argparse
import os
import sys
from pathlib import Path
from Bio import SeqIO
from Bio.PDB import PDBParser
import torch
import numpy as np
import pandas as pd
from transformers import EsmForProteinFolding, AutoTokenizer


def setup_esmfold(device='cpu'):
    """Initialize ESMFold model"""
    print("Loading ESMFold model...")
    model = EsmForProteinFolding.from_pretrained("facebook/esmfold_v1")
    model.eval()
    model.to(device)
    
    # Load tokenizer
    tokenizer = AutoTokenizer.from_pretrained("facebook/esmfold_v1")
    
    return model, tokenizer


def predict_structure(model, tokenizer, sequence, sequence_id, device='cpu'):
    """Predict protein structure for a single sequence"""
    print(f"Predicting structure for {sequence_id}...")
    
    # Tokenize sequence
    tokenized = tokenizer(sequence, return_tensors="pt", add_special_tokens=False)
    tokenized = {k: v.to(device) for k, v in tokenized.items()}
    
    # Run ESMFold prediction
    with torch.no_grad():
        output = model(**tokenized)
    
    # Extract PDB format from the output
    # ESMFold returns coordinates and confidence, we need to convert to PDB
    pdb_output = model.infer_pdb(sequence)
    
    # Extract coordinates and confidence scores from the model output
    # ESMFold output.positions shape: [8, 1, L, 14, 3] - 8 models, 1 sequence, L residues, 14 atoms, 3 coords
    # We'll use the first model and extract CA coordinates (atom index 1)
    coords = output.positions[0, 0, :, 1, :].cpu().numpy()  # [L, 3] - CA coordinates
    confidence = output.plddt[0, 0, :].cpu().numpy()  # [L] - per-residue confidence
    
    return {
        'pdb': pdb_output,
        'coordinates': coords,  # [L, 3] CA coordinates
        'confidence': confidence,  # [L] per-residue confidence
        'sequence_id': sequence_id,
        'sequence': sequence
    }


def get_pdb_path(sequence_id, output_dir):
    """Get the expected PDB file path for a sequence ID"""
    clean_id = "".join(c for c in sequence_id if c.isalnum() or c in ('-', '_'))
    pdb_filename = f"{clean_id}.pdb"
    return os.path.join(output_dir, pdb_filename)


def save_pdb(output_dict, sequence_id, output_dir):
    """Save PDB structure to file"""
    pdb_path = get_pdb_path(sequence_id, output_dir)
    
    # Write PDB file
    with open(pdb_path, 'w') as f:
        f.write(output_dict['pdb'])
    
    print(f"Saved structure to {pdb_path}")
    return pdb_path


def extract_coordinates_from_pdb(pdb_path):
    """Extract CA (alpha carbon) coordinates from a PDB file"""
    try:
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure('protein', pdb_path)
        
        # Extract CA coordinates from all chains
        ca_coords = []
        for model in structure:
            for chain in model:
                for residue in chain:
                    # Check if residue has CA atom
                    if 'CA' in residue:
                        ca_atom = residue['CA']
                        ca_coords.append(ca_atom.coord)
        
        if len(ca_coords) == 0:
            return None
        
        # Convert to numpy array
        coords_array = np.array(ca_coords)
        return coords_array
    except Exception as e:
        print(f"  Warning: Could not extract coordinates from PDB: {e}")
        return None


def load_existing_coordinates(output_dir):
    """Load existing coordinates from NPZ file if it exists"""
    npz_path = os.path.join(output_dir, "esmfold_coordinates.npz")
    if os.path.exists(npz_path):
        try:
            coords_dict = np.load(npz_path)
            # Convert to regular dict (np.load returns a NpzFile object)
            coords_dict = {key: coords_dict[key] for key in coords_dict.files}
            print(f"Loaded {len(coords_dict)} existing coordinate sets from {npz_path}")
            return coords_dict
        except Exception as e:
            print(f"Warning: Could not load existing coordinates: {e}")
            return {}
    return {}


def save_coordinates(all_coordinates, sequence_id_mapping, output_dir):
    """Save all coordinates in a single file for PCA analysis"""
    # Create output directory if needed
    os.makedirs(output_dir, exist_ok=True)
    
    # Save as NPZ (NumPy compressed format) - efficient for large datasets
    npz_path = os.path.join(output_dir, "esmfold_coordinates.npz")
    np.savez_compressed(npz_path, **all_coordinates)
    print(f"\nSaved all coordinates to {npz_path}")
    
    # Also save metadata as CSV with original sequence IDs
    metadata = []
    for clean_id, coords in all_coordinates.items():
        original_id = sequence_id_mapping.get(clean_id, clean_id)
        metadata.append({
            'clean_id': clean_id,  # Key used in NPZ file
            'original_sequence_id': original_id,  # Original FASTA sequence ID
            'n_residues': coords.shape[0],
            'coords_shape': f"{coords.shape[0]}x{coords.shape[1]}"
        })
    
    metadata_df = pd.DataFrame(metadata)
    metadata_path = os.path.join(output_dir, "coordinate_metadata.csv")
    metadata_df.to_csv(metadata_path, index=False)
    print(f"Saved coordinate metadata to {metadata_path}")
    
    return npz_path, metadata_path


def process_fasta(input_fasta, output_dir, max_sequences=None, device='cpu', force=False):
    """Process FASTA file and generate PDB files"""
    
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Read FASTA file
    sequences = list(SeqIO.parse(input_fasta, "fasta"))
    total_sequences = len(sequences)
    
    if max_sequences:
        sequences = sequences[:max_sequences]
        print(f"Processing first {len(sequences)} sequences (limited by --max-sequences)")
    else:
        print(f"Processing all {total_sequences} sequences")
    
    # Load existing coordinates if available
    existing_coords = load_existing_coordinates(output_dir) if not force else {}
    
    # Check if we need to load the model (only if we need to make new predictions)
    need_model = False
    for record in sequences:
        sequence_id = record.id
        pdb_path = get_pdb_path(sequence_id, output_dir)
        if force or not os.path.exists(pdb_path):
            need_model = True
            break
    
    # Load ESMFold model only if needed
    model = None
    tokenizer = None
    if need_model:
        model, tokenizer = setup_esmfold(device)
    
    # Process each sequence
    successful_predictions = 0
    failed_predictions = 0
    skipped_predictions = 0
    all_coordinates = {}  # Dictionary to store all coordinates: {clean_id: coords_array}
    sequence_id_mapping = {}  # Map clean_id to original sequence_id
    
    for i, record in enumerate(sequences, 1):
        sequence_id = record.id
        sequence = str(record.seq)
        clean_id = "".join(c for c in sequence_id if c.isalnum() or c in ('-', '_'))
        pdb_path = get_pdb_path(sequence_id, output_dir)
        
        # Check if PDB file already exists
        if not force and os.path.exists(pdb_path):
            print(f"\n[{i}/{len(sequences)}] Skipping {sequence_id} (PDB file already exists)")
            # Try to load coordinates from existing NPZ file
            if clean_id in existing_coords:
                all_coordinates[clean_id] = existing_coords[clean_id]
                sequence_id_mapping[clean_id] = sequence_id
                skipped_predictions += 1
            else:
                # Try to extract coordinates from the PDB file
                print(f"  Extracting coordinates from PDB file...")
                coords = extract_coordinates_from_pdb(pdb_path)
                if coords is not None:
                    all_coordinates[clean_id] = coords
                    sequence_id_mapping[clean_id] = sequence_id
                    print(f"  Successfully extracted {coords.shape[0]} CA coordinates from PDB")
                else:
                    print(f"  Warning: Could not extract coordinates from PDB file")
                skipped_predictions += 1
            continue
        
        print(f"\n[{i}/{len(sequences)}] Processing {sequence_id}")
        print(f"Sequence length: {len(sequence)} amino acids")
        
        try:
            # Predict structure (returns dict with pdb, coordinates, confidence)
            output = predict_structure(model, tokenizer, sequence, sequence_id, device)
            
            # Save PDB file
            save_pdb(output, sequence_id, output_dir)
            
            # Store coordinates (using clean sequence ID as key for NPZ compatibility)
            all_coordinates[clean_id] = output['coordinates']
            sequence_id_mapping[clean_id] = sequence_id  # Store mapping
            
            successful_predictions += 1
            
        except Exception as e:
            print(f"Error processing {sequence_id}: {str(e)}")
            failed_predictions += 1
            continue
    
    # Save all coordinates in a single file
    if all_coordinates:
        save_coordinates(all_coordinates, sequence_id_mapping, output_dir)
    
    # Summary
    print(f"\n{'='*50}")
    print(f"ESMFold Prediction Summary")
    print(f"{'='*50}")
    print(f"Total sequences processed: {len(sequences)}")
    print(f"Successful predictions: {successful_predictions}")
    print(f"Skipped (already exist): {skipped_predictions}")
    print(f"Failed predictions: {failed_predictions}")
    print(f"Output directory: {output_dir}")
    
    if all_coordinates:
        print(f"\nCoordinates saved to:")
        print(f"- {os.path.join(output_dir, 'esmfold_coordinates.npz')} ({len(all_coordinates)} coordinate sets)")
        print(f"- {os.path.join(output_dir, 'coordinate_metadata.csv')} (metadata)")
    
    if successful_predictions > 0:
        print(f"\nPDB files saved in: {output_dir}")
        print("You can view these structures using:")
        print("- PyMOL: pymol *.pdb")
        print("- ChimeraX: chimerax *.pdb")
        print("- VMD: vmd *.pdb")


def main():
    parser = argparse.ArgumentParser(
        description="Predict protein structures using ESMFold",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    # Process all sequences in a FASTA file
    python esmfold_predict.py input.fasta -o structures/
    
    # Process only first 10 sequences
    python esmfold_predict.py input.fasta -o structures/ --max-sequences 10
    
    # Process with custom output directory
    python esmfold_predict.py my_proteins.fasta -o my_structures/
        """
    )
    
    parser.add_argument(
        "input_fasta",
        help="Input FASTA file containing protein sequences"
    )
    
    parser.add_argument(
        "-o", "--output-dir",
        default="esmfold_structures",
        help="Output directory for PDB files (default: esmfold_structures)"
    )
    
    parser.add_argument(
        "--max-sequences",
        type=int,
        help="Maximum number of sequences to process (useful for testing)"
    )
    
    parser.add_argument(
        "--gpu",
        action="store_true",
        help="Use GPU if available (requires CUDA)"
    )
    
    parser.add_argument(
        "--force",
        action="store_true",
        help="Force re-prediction even if PDB files already exist"
    )
    
    args = parser.parse_args()
    
    # Check if input file exists
    if not os.path.exists(args.input_fasta):
        print(f"Error: Input file '{args.input_fasta}' not found")
        sys.exit(1)
    
    # Check if input is a valid FASTA file
    try:
        sequences = list(SeqIO.parse(args.input_fasta, "fasta"))
        if len(sequences) == 0:
            print(f"Error: No sequences found in '{args.input_fasta}'")
            sys.exit(1)
    except Exception as e:
        print(f"Error reading FASTA file: {str(e)}")
        sys.exit(1)
    
    # Set device
    if args.gpu and torch.cuda.is_available():
        device = torch.device("cuda")
        print(f"Using GPU: {torch.cuda.get_device_name()}")
    else:
        device = torch.device("cpu")
        print("Using CPU")
    
    # Process the FASTA file
    try:
        process_fasta(args.input_fasta, args.output_dir, args.max_sequences, device, args.force)
    except KeyboardInterrupt:
        print("\n\nPrediction interrupted by user")
        sys.exit(1)
    except Exception as e:
        print(f"\nError during prediction: {str(e)}")
        sys.exit(1)


if __name__ == "__main__":
    main()
