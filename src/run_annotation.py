#!/usr/bin/env python
# run_annotation.py

import os
import argparse
import subprocess
from pathlib import Path

# Species database mapping
SPECIES_DB = {
    "1": {"name": "Human", "orgdb": "org.Hs.eg.db", "kegg": "hsa", "version": "3.13.0"},
    "2": {"name": "Mouse", "orgdb": "org.Mm.eg.db", "kegg": "mmu", "version": "3.13.0"},
    "3": {"name": "Rat", "orgdb": "org.Rn.eg.db", "kegg": "rno", "version": "3.13.0"},
    "4": {"name": "Arabidopsis", "orgdb": "org.At.tair.db", "kegg": "ath", "version": "3.18.1"},
    "5": {"name": "Rice", "orgdb": "org.Osativa.eg.db", "kegg": "osa", "version": "3.18.0"},
    "6": {"name": "Maize", "orgdb": "org.Zmays.eg.db", "kegg": "zma", "version": "3.18.1"},
    "7": {"name": "Fruit Fly", "orgdb": "org.Dm.eg.db", "kegg": "dme", "version": "3.13.0"},
    "8": {"name": "Nematode", "orgdb": "org.Ce.eg.db", "kegg": "cel", "version": "3.13.0"},
}


def find_recirc_file(output_dir):
    """Finds the *_recirc.txt file in the given output directory"""
    output_dir = Path(output_dir)
    recirc_files = list(output_dir.glob("*_recirc.txt"))

    if not recirc_files:
        print(f"ERROR: No *_recirc.txt file found in {output_dir}")
        return None
    if len(recirc_files) > 1:
        print(f"ERROR: Multiple *_recirc.txt files found: {recirc_files}")
        return None

    return recirc_files[0]


def validate_input(input_file):
    """Validates the input file format (must contain 'geneName' column)"""
    path = Path(input_file)

    if not path.exists():
        print(f"ERROR: File not found: {input_file}")
        return False
    if path.stat().st_size == 0:
        print(f"ERROR: File is empty: {input_file}")
        return False

    # Check if 'geneName' column exists
    try:
        with open(path, 'r') as f:
            header = f.readline().strip().split('\t')
        if 'geneName' not in header:
            print("ERROR: 'geneName' column is missing in the input file.")
            print("Detected columns:", ', '.join(header))
            return False
    except Exception as e:
        print(f"ERROR: Failed to read file: {str(e)}")
        return False

    return True


def select_species():
    """Prompts the user to select a species"""
    print("\nPlease select a species for annotation:")
    for key, info in SPECIES_DB.items():
        print(f"  [{key}] {info['name']:12} | OrgDb: {info['orgdb']} (v{info['version']}) | KEGG: {info['kegg']}")

    while True:
        choice = input("\nEnter species number (1-8) or 'q' to quit: ").strip()
        if choice.lower() == 'q':
            print("Exiting annotation process.")
            exit(0)
        if choice in SPECIES_DB:
            return SPECIES_DB[choice]
        print("Invalid selection. Please enter a number between 1 and 8.")


def run_r_analysis(species_info, input_file, output_dir):
    """Runs the R script for functional enrichment analysis"""
    r_script = Path(__file__).parent / "anno/annotation.R"

    if not r_script.exists():
        print(f"ERROR: R script not found at {r_script}")
        exit(1)

    cmd = [
        "Rscript", str(r_script),
        "--input", str(input_file),
        "--orgdb", species_info["orgdb"],
        "--kegg", species_info["kegg"],
        "--output", str(output_dir)
    ]

    print("\nStarting functional enrichment analysis...")
    try:
        subprocess.run(cmd, check=True)
        print("\nAnnotation completed. Results saved in:", output_dir)
    except subprocess.CalledProcessError as e:
        print(f"ERROR: R script execution failed: {e}")
        exit(1)


def main(output_dir=None):
    """Main function to execute the annotation process"""
    if output_dir is None:
        output_dir = "detection_output"  # Default path if none provided

    output_dir = Path(output_dir)

    # 1. Locate the input file
    input_file = find_recirc_file(output_dir)
    if not input_file:
        exit(1)

    # 2. Validate input file format
    if not validate_input(input_file):
        exit(1)

    # 3. Let the user select a species
    species_info = select_species()
    print(f"\nSelected species: {species_info['name']}")

    # 4. Prepare output directory
    annotation_output = output_dir / "annotation_results"
    annotation_output.mkdir(parents=True, exist_ok=True)

    # 5. Run enrichment analysis
    run_r_analysis(species_info, input_file, annotation_output)


if __name__ == "__main__":
    import sys

    if len(sys.argv) > 1:
        main(sys.argv[1])
    else:
        main()  # Use default path if no argument provided
