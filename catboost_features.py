# -*- coding: utf-8 -*-
"""
Created on Thu Jun  5 21:46:04 2025

@author: jbren

"""
import os
import sys
import subprocess
import csv
import numpy as np
import pandas as pd
import argparse
from Bio.PDB import PDBParser
from Bio.PDB.DSSP import DSSP

# =============================================================================
# 1. CONFIGURATION - PLEASE EDIT THESE PATHS
# =============================================================================
# --- Dataset Path Configurations ---
# Define the file paths for each dataset. The script will choose which to use
# based on the command-line argument ('cancer' or 'tsu').
DATASET_PATHS = {
    'cancer': {
        'input_csv': r"C:\Users\jbren\Documents\Documents\phosphomutations\final_cancer.csv",
        'pdb_path': r'C:\Users\jbren\Documents\Documents\phosphomutations\protein\structural_information\pdb_cancer',
        'output_csv': r"C:\Users\jbren\Documents\Documents\phosphomutations\final_cancer_foldx_ss.csv"
    },
    'tsu': {
        'input_csv': r"C:\Users\jbren\Documents\Documents\phosphomutations\tsui.csv",
        'pdb_path': r'C:\Users\jbren\Documents\Documents\phosphomutations\protein\structural_information\pdb_files',
        'output_csv': r"C:\Users\jbren\Documents\Documents\phosphomutations\tsui_foldx_ss.csv"
    }
}

# --- Executable and Data Paths ---
# Path to the DSSP executable
DSSP_EXE = r"C:\Users\jbren\Documents\Documents\phosphomutations\dssp-2.2.1-win64.exe"
# Path to the FoldX executable
FOLDX_EXE = r"C:\path\to\your\foldx_20251231.exe" # <-- IMPORTANT: UPDATE THIS PATH

# --- FoldX Specific Paths ---
# These are relative to where you run the script or can be absolute paths.
REPAIRED_PDB_PATH = './foldx_repair'
MUTANT_PDB_PATH = './Mutant_Structures'
INDIVIDUAL_LIST_FILENAME = 'individual_list.txt' # Temporary file for FoldX mutations

# --- Analysis Flags ---
# Set to True to run the FoldX repair command (can be slow)
RUN_FOLDX_REPAIR = False
# Set to True to run the FoldX BuildModel command
RUN_FOLDX_BUILDMODEL = True


# =============================================================================
# 2. ATOM CLASSIFICATION DEFINITIONS (from position.py)
# =============================================================================
# Based on YRB (PMC4602141) and general chemical properties
YRB_POSITIVE_CHARGE_SOURCES = {
    ('LYS', 'NZ'), ('ARG', 'NE'), ('ARG', 'NH1'), ('ARG', 'NH2')
}
YRB_NEGATIVE_CHARGE_SOURCES = {
    ('ASP', 'OD1'), ('ASP', 'OD2'), ('GLU', 'OE1'), ('GLU', 'OE2')
}
YRB_HYDROPHOBIC_CARBONS = {
    ('ALA', 'CB'), ('ARG', 'CB'), ('ARG', 'CG'), ('ASN', 'CB'), ('ASP', 'CB'),
    ('CYS', 'CB'), ('GLN', 'CB'), ('GLN', 'CG'), ('GLU', 'CB'), ('GLU', 'CG'),
    ('HIS', 'CB'), ('HIS', 'CG'), ('HIS', 'CD2'), ('HIS', 'CE1'), ('ILE', 'CB'),
    ('ILE', 'CG1'), ('ILE', 'CG2'), ('ILE', 'CD1'), ('LEU', 'CB'), ('LEU', 'CG'),
    ('LEU', 'CD1'), ('LEU', 'CD2'), ('LYS', 'CB'), ('LYS', 'CG'), ('LYS', 'CD'),
    ('MET', 'CB'), ('MET', 'CG'), ('MET', 'CE'), ('PHE', 'CB'), ('PHE', 'CG'),
    ('PHE', 'CD1'), ('PHE', 'CD2'), ('PHE', 'CE1'), ('PHE', 'CE2'), ('PHE', 'CZ'),
    ('PRO', 'CB'), ('PRO', 'CG'), ('PRO', 'CD'), ('THR', 'CG2'), ('TRP', 'CB'),
    ('TRP', 'CG'), ('TRP', 'CD2'), ('TRP', 'CE3'), ('TRP', 'CZ2'), ('TRP', 'CZ3'),
    ('TRP', 'CH2'), ('TYR', 'CB'), ('TYR', 'CG'), ('TYR', 'CD1'), ('TYR', 'CD2'),
    ('TYR', 'CE1'), ('TYR', 'CE2'), ('TYR', 'CZ'), ('VAL', 'CB'), ('VAL', 'CG1'),
    ('VAL', 'CG2')
}
SIDECHAIN_POLAR_ATOMS_SPECIFIC = {
    ('ASN', 'OD1'), ('ASN', 'ND2'), ('CYS', 'SG'), ('GLN', 'OE1'), ('GLN', 'NE2'),
    ('HIS', 'ND1'), ('HIS', 'NE2'), ('MET', 'SD'), ('SER', 'OG'), ('THR', 'OG1'),
    ('TRP', 'NE1'), ('TYR', 'OH')
}


# =============================================================================
# 3. CLASSES AND HELPER FUNCTIONS
# =============================================================================

class ProteinSecondaryStructureAnalyzer:
    """
    Analyzes protein secondary structure, solvent accessibility, and local
    environment for a specific residue. (Derived from position.py)
    """
    def __init__(self, pdb_file, dssp_executable):
        self.pdb_file = pdb_file
        if not os.path.exists(pdb_file):
            raise FileNotFoundError(f"PDB file not found: {pdb_file}")
        
        parser = PDBParser(QUIET=True)
        structure_id = os.path.basename(pdb_file).replace('.pdb', '')
        self.structure = parser.get_structure(structure_id, pdb_file)
        self.model = self.structure[0]
        self.dssp = DSSP(self.model, pdb_file, dssp=dssp_executable)
        self.position = None
        self.alt_position = None
        self.is_beta_hairpin = False
        self.ss = None
        self.RelSASA = np.nan

    # ... (All helper methods from the original class are included below) ...
    def _get_ca_coordinates(self, residue_number):
        try:
            return self.model['A'][residue_number]['CA'].get_coord()
        except KeyError:
            return None

    def _is_current_sheet_hairpin(self, sheet_start, sheet_end):
        if sheet_end - sheet_start < 2: return False
        c_term_residues = [sheet_end - 2, sheet_end - 1, sheet_end]
        search_end = min(sheet_end + 20, sheet_end + 50)
        for search_res in range(sheet_end + 1, search_end):
            try:
                search_key = ('A', (' ', search_res, ' '))
                if self.dssp[search_key][2] in ['E', 'B']:
                    if self._check_hairpin_distance(c_term_residues, search_res):
                        return True
            except KeyError:
                continue
        return False

    def _check_hairpin_distance(self, c_term_residues, other_sheet_residue):
        other_coord = self._get_ca_coordinates(other_sheet_residue)
        if other_coord is None: return False
        for res_id in c_term_residues:
            c_term_coord = self._get_ca_coordinates(res_id)
            if c_term_coord is not None:
                if np.linalg.norm(c_term_coord - other_coord) <= 5.0:
                    return True
        return False

    def _is_n_terminal_hairpin_sheet(self, sheet_start, sheet_end):
        search_end = min(sheet_end + 50, sheet_end + 100)
        for search_res in range(sheet_end + 1, search_end):
            try:
                search_key = ('A', (' ', search_res, ' '))
                if self.dssp[search_key][2] in ['E', 'B']:
                    return True
            except KeyError:
                continue
        return False

    def analyze_residue(self, residue_number):
        residue_number = int(residue_number)
        try:
            dssp_key = ('A', (' ', residue_number, ' '))
            self.ss, self.RelSASA = self.dssp[dssp_key][2], self.dssp[dssp_key][3]
        except KeyError:
            return

        if self.ss == '-':
            self.is_beta_hairpin = False
            return

        if self.ss in ['H', 'G', 'I', 'T', 'S']:
            position, current_residue = 1, residue_number - 1
            while current_residue > 0:
                try:
                    if self.dssp[('A', (' ', current_residue, ' '))][2] == self.ss:
                        position += 1
                        current_residue -= 1
                    else: break
                except KeyError: break
            self.is_beta_hairpin = False
            self.position = position
            self.alt_position = min(self.position - 1, 4)

        elif self.ss in ['E', 'B']:
            sheet_start, sheet_end = residue_number, residue_number
            # Find N-terminal end
            res_cursor = residue_number - 1
            while res_cursor > 0:
                try:
                    if self.dssp[('A', (' ', res_cursor, ' '))][2] in ['E', 'B']:
                        sheet_start = res_cursor
                        res_cursor -= 1
                    else: break
                except KeyError: break
            # Find C-terminal end
            res_cursor = residue_number + 1
            while True:
                try:
                    if self.dssp[('A', (' ', res_cursor, ' '))][2] in ['E', 'B']:
                        sheet_end = res_cursor
                        res_cursor += 1
                    else: break
                except KeyError: break
            
            self.is_beta_hairpin = self._is_current_sheet_hairpin(sheet_start, sheet_end)
            if not self.is_beta_hairpin:
                self.position = sheet_end - residue_number + 1
            else:
                if self._is_n_terminal_hairpin_sheet(sheet_start, sheet_end):
                    self.position = sheet_end - residue_number + 1
                else:
                    self.position = residue_number - sheet_start + 1

    def count_atoms_by_property_in_radius(self, target_res_num, radius):
        counts = {'hydrophobic_carbons': 0, 'positive_charge_source_atoms': 0,
                  'negative_charge_source_atoms': 0, 'polar_atoms': 0, 'total_atoms': 0}
        try:
            target_residue = self.model['A'][(' ', int(target_res_num), ' ')]
        except KeyError:
            return counts
        
        target_coords = [atom.get_coord() for atom in target_residue]
        if not target_coords: return counts
        
        center = np.mean(target_coords, axis=0)
        
        for residue in self.model.get_residues():
            if residue.id[0] != ' ': continue # Skip HETATMs
            res_name = residue.get_resname()
            for atom in residue:
                if np.linalg.norm(atom.get_coord() - center) <= radius:
                    counts['total_atoms'] += 1
                    atom_id, atom_elem = atom.get_id()[0], atom.element.upper()
                    classified = False
                    if (res_name, atom_id) in YRB_POSITIVE_CHARGE_SOURCES:
                        counts['positive_charge_source_atoms'] += 1; classified = True
                    elif (res_name, atom_id) in YRB_NEGATIVE_CHARGE_SOURCES:
                        counts['negative_charge_source_atoms'] += 1; classified = True
                    elif atom_elem == 'C' and (res_name, atom_id) in YRB_HYDROPHOBIC_CARBONS:
                        counts['hydrophobic_carbons'] += 1; classified = True
                    
                    if not classified:
                        if atom_id in ('N', 'O') or (res_name, atom_id) in SIDECHAIN_POLAR_ATOMS_SPECIFIC:
                            counts['polar_atoms'] += 1
        return counts

def run_foldx_command(command_string):
    """A helper to run a FoldX command and handle errors."""
    print(f"Running FoldX command: {command_string}")
    try:
        result = subprocess.run(command_string, shell=False, capture_output=True, text=True, check=True)
        print("FoldX ran successfully.")
        # print(result.stdout) # Uncomment for verbose output
    except subprocess.CalledProcessError as e:
        print(f"ERROR: FoldX command failed.")
        print(f"Stderr: {e.stderr}")
        print(f"Stdout: {e.stdout}")
        raise  # Re-raise the exception to stop the script if a step fails
    except FileNotFoundError:
        print(f"ERROR: Could not find FoldX executable. Please check FOLDX_EXE path.")
        raise

def repair_pdb(pdb_id, pdb_path):
    """
    Runs the FoldX RepairPDB command on a given PDB file.
    (Derived from foldx_stability.py)
    """
    if not os.path.isdir(REPAIRED_PDB_PATH):
        print(f"Creating directory for repaired PDBs: {REPAIRED_PDB_PATH}")
        os.makedirs(REPAIRED_PDB_PATH)

    repaired_pdb_file = f"{pdb_id}_Repair.pdb"
    
    # Run RepairPDB
    cmd_repair = (
        f"{FOLDX_EXE} --command=RepairPDB --pdb-dir={pdb_path} "
        f"--pdb={pdb_id}.pdb --output-dir={REPAIRED_PDB_PATH} "
        f"--output-file={repaired_pdb_file} --water=IGNORE"
    )
    run_foldx_command(cmd_repair)
    
    # Run Optimize on the repaired structure
    cmd_optimize = (
        f"{FOLDX_EXE} --command=Optimize --pdb-dir={REPAIRED_PDB_PATH} "
        f"--pdb={repaired_pdb_file} --output-dir={REPAIRED_PDB_PATH} "
        f"--water=IGNORE"
    )
    run_foldx_command(cmd_optimize)
    print(f"Finished repairing and optimizing {pdb_id}")
    return f"Optimized_{repaired_pdb_file}"


def build_model_and_analyze(repaired_pdb_filename, mutation_str):
    """
    Runs FoldX BuildModel and parses the output to get stability changes.
    (Derived from foldx_stability.py)
    """
    if not os.path.isdir(MUTANT_PDB_PATH):
        print(f"Creating directory for mutant PDBs: {MUTANT_PDB_PATH}")
        os.makedirs(MUTANT_PDB_PATH)

    # Create the individual list file for the mutation
    # Format: OriginalResidue Chain ResidueNumber NewResidue; e.g., SA123p;
    # Assuming Chain A and phosphoserine/threonine mutation format
    original_res = mutation_str[0]
    res_num = mutation_str[1:-1]
    new_res_code = 'p' if original_res in ['S', 'T'] else 'a' # p for phospho, a for alanine
    foldx_mutation = f"{original_res}A{res_num}{new_res_code};"

    with open(INDIVIDUAL_LIST_FILENAME, "w") as f:
        f.write(foldx_mutation)

    # Run BuildModel
    output_prefix = f"{os.path.splitext(repaired_pdb_filename)[0]}_{mutation_str}"
    cmd_build = (
        f"{FOLDX_EXE} --command=BuildModel --pdb-dir={REPAIRED_PDB_PATH} "
        f"--pdb={repaired_pdb_filename} --mutant-file={INDIVIDUAL_LIST_FILENAME} "
        f"--output-dir={MUTANT_PDB_PATH} --output-file={output_prefix} "
        f"--water=IGNORE --pH=7.3 --ionStrength=0.15 --out-pdb=true"
    )
    run_foldx_command(cmd_build)

    # Parse the Dif_ file for energy terms
    dif_file_path = os.path.join(MUTANT_PDB_PATH, f"Dif_{output_prefix}.fxout")
    if not os.path.exists(dif_file_path):
        print(f"Warning: FoldX output file not found: {dif_file_path}")
        return {}

    with open(dif_file_path, 'r') as file:
        lines = file.readlines()
    
    if len(lines) < 10:
        print(f"Warning: FoldX output file is malformed: {dif_file_path}")
        return {}

    header = [h.strip() for h in lines[8].split('\t')]
    values = [v.strip() for v in lines[9].split('\t')]
    
    # Prefixing column names with 'FoldX_' to avoid conflicts
    foldx_results = {f"FoldX_{h}": v for h, v in zip(header, values)}
    return foldx_results


# =============================================================================
# 4. MAIN EXECUTION
# =============================================================================
def main():
    """
    Main function to orchestrate the analysis pipeline.
    """
    # --- Command-line Argument Parsing ---
    parser = argparse.ArgumentParser(description="Run protein analysis pipeline.")
    subparsers = parser.add_subparsers(dest='mode', required=True, help="Operating mode")

    # Create the parser for the "dataset" command
    parser_dataset = subparsers.add_parser('dataset', help='Process a full dataset from a predefined CSV file.')
    parser_dataset.add_argument('name', choices=DATASET_PATHS.keys(), 
                                help=f"The name of the dataset to process: {list(DATASET_PATHS.keys())}")

    # Create the parser for the "single" command
    parser_single = subparsers.add_parser('single', help='Process a single PDB file.')
    parser_single.add_argument('--pdb_file', required=True, help='Full path to the input PDB file.')
    parser_single.add_argument('--mutation', required=True, help="Mutation to analyze, e.g., 'S123A'.")
    parser_single.add_argument('--output_csv', required=True, help='Path for the output CSV file.')

    args = parser.parse_args()

    # --- Set up DataFrame and paths based on mode ---
    if args.mode == 'dataset':
        paths = DATASET_PATHS[args.name]
        INPUT_CSV_PATH = paths['input_csv']
        PDB_PATH = paths['pdb_path']
        OUTPUT_CSV_PATH = paths['output_csv']
        
        print(f"--- Running analysis for dataset: {args.name.upper()} ---")
        print(f"Input CSV: {INPUT_CSV_PATH}")
        
        if not os.path.exists(INPUT_CSV_PATH):
            print(f"FATAL: Input CSV not found at {INPUT_CSV_PATH}")
            sys.exit(1)
        df = pd.read_csv(INPUT_CSV_PATH)
        print(f"Loaded {len(df)} rows from {INPUT_CSV_PATH}")

    elif args.mode == 'single':
        PDB_FILE_PATH = args.pdb_file
        PDB_PATH = os.path.dirname(PDB_FILE_PATH)
        pdb_id = os.path.basename(PDB_FILE_PATH).replace('.pdb', '')
        mutation = args.mutation
        OUTPUT_CSV_PATH = args.output_csv

        print("--- Running analysis for single PDB file ---")
        print(f"PDB file: {PDB_FILE_PATH}")
        print(f"Mutation: {mutation}")

        if not os.path.exists(PDB_FILE_PATH):
            print(f"FATAL: PDB file not found at {PDB_FILE_PATH}")
            sys.exit(1)
            
        # Create a DataFrame with a single row for processing
        df = pd.DataFrame([{'pdb': pdb_id, 'mutation': mutation}])
    
    print(f"PDB Path: {PDB_PATH}")
    print(f"Output CSV: {OUTPUT_CSV_PATH}")

    # --- Prepare for Results ---
    all_results = []

    # --- Loop and Process Each Row ---
    for index, row in df.iterrows():
        print("-" * 50)
        print(f"Processing row {index+1}/{len(df)}...")

        pdb_id = row.get('pdb')
        mutation = row.get('mutation')

        if pd.isna(pdb_id) or pd.isna(mutation):
            print(f"  Skipping row {index} due to missing PDB ID or mutation.")
            all_results.append({}) # Append empty dict to maintain row alignment
            continue

        print(f"  PDB: {pdb_id}, Mutation: {mutation}")
        current_row_results = {}
        
        try:
            # --- Part 1: Position and Structure Analysis ---
            pdb_file_path = os.path.join(PDB_PATH, f"{pdb_id}.pdb")
            residue_number = int(mutation[1:-1])
            residue_type = mutation[0]

            analyzer = ProteinSecondaryStructureAnalyzer(pdb_file_path, DSSP_EXE)
            analyzer.analyze_residue(residue_number)

            sasa = np.nan
            if residue_type == 'Y': sasa = analyzer.RelSASA * 222
            elif residue_type == 'T': sasa = analyzer.RelSASA * 142
            elif residue_type == 'S': sasa = analyzer.RelSASA * 130

            atom_5A = analyzer.count_atoms_by_property_in_radius(residue_number, 5.0)
            atom_10A = analyzer.count_atoms_by_property_in_radius(residue_number, 10.0)

            current_row_results.update({
                'SS': analyzer.ss, 'RelSASA': analyzer.RelSASA, 'SASA': sasa,
                'Position': analyzer.position, 'AltPosition': analyzer.alt_position,
                'Hairpin': analyzer.is_beta_hairpin,
                'Polar_5A': atom_5A['polar_atoms'],
                'Hydrophobic_5A': atom_5A['hydrophobic_carbons'],
                'Positive_5A': atom_5A['positive_charge_source_atoms'],
                'Negative_5A': atom_5A['negative_charge_source_atoms'],
                'Total_5A': atom_5A['total_atoms'],
                'Polar_10A': atom_10A['polar_atoms'],
                'Hydrophobic_10A': atom_10A['hydrophobic_carbons'],
                'Positive_10A': atom_10A['positive_charge_source_atoms'],
                'Negative_10A': atom_10A['negative_charge_source_atoms'],
                'Total_10A': atom_10A['total_atoms'],
            })
            print(f"  Structure analysis complete. SS: {analyzer.ss}, Position: {analyzer.position}")

            # --- Part 2: FoldX Stability Analysis ---
            repaired_pdb_filename = f"Optimized_{pdb_id}_Repair.pdb"
            repaired_pdb_full_path = os.path.join(REPAIRED_PDB_PATH, repaired_pdb_filename)

            if RUN_FOLDX_REPAIR or not os.path.exists(repaired_pdb_full_path):
                print("  Running FoldX RepairPDB...")
                repaired_pdb_filename = repair_pdb(pdb_id, PDB_PATH)
            else:
                print("  Skipping FoldX repair, using existing file.")
            
            if RUN_FOLDX_BUILDMODEL:
                print("  Running FoldX BuildModel...")
                foldx_energy_terms = build_model_and_analyze(repaired_pdb_filename, mutation)
                current_row_results.update(foldx_energy_terms)
                print(f"  FoldX analysis complete. Found {len(foldx_energy_terms)} energy terms.")

        except Exception as e:
            print(f"  !! An error occurred while processing row {index}: {e}")
        
        all_results.append(current_row_results)

    # --- Combine and Save Results ---
    print("-" * 50)
    print("Combining results with original DataFrame...")
    
    # For single mode, the original df is just a placeholder. 
    # For dataset mode, we want to append to the original.
    if args.mode == 'single':
        final_df = pd.DataFrame(all_results)
    else: # dataset mode
        results_df = pd.DataFrame(all_results, index=df.index)
        final_df = pd.concat([df, results_df], axis=1)

    # Save to new CSV
    final_df.to_csv(OUTPUT_CSV_PATH, index=False)
    print(f"Processing complete. All results saved to: {OUTPUT_CSV_PATH}")


if __name__ == "__main__":
    main()
