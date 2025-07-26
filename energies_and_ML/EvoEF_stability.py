import os
import shutil
import subprocess
import csv




# def run_command(command):
#     """Utility function to run a command and handle exceptions."""
#     try:
#         subprocess.run(command, check=True, capture_output=True)
#     except subprocess.CalledProcessError as e:
#         print(f"An error occurred: {e}")

def process_pdb(pdb_path, evoef_path, output_directory,mutant,individual_file_path):
    """Process a single PDB file through the repair, mutation, and ΔΔG calculation steps."""
    repair_pdb(pdb_path, evoef_path,output_directory)
    return deltadeltaG(pdb_path, evoef_path,output_directory,mutant,individual_file_path)


def repair_pdb(pdb_path, evoef_path, output_directory):
    # Ensure the output directory exists
    os.makedirs(output_directory, exist_ok=True)
    
    # Prepare the output file name
    base_filename = os.path.basename(pdb_path).replace('.pdb', '')
    repaired_pdb = f"{output_directory}/{base_filename}_Repair.pdb"
    
    # Check if the file has already been repaired and moved
    if os.path.isfile(repaired_pdb):
        print("Already Prepared")
    else:
        # Command to repair the structure
        repair_command = [
            evoef_path, "--command=RepairStructure", f"--pdb={pdb_path}"
        ]
        print(f"Repairing {pdb_path}...")
        
        # Run the command and capture output
        process = subprocess.run(repair_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        if process.returncode == 0:
            print("Repair completed successfully.")
        else:
            print(f"Error in repairing PDB: {process.stderr}")
        move_command = ["mv",f"{base_filename}_Repair.pdb",f"{output_directory}"]
        # print(move_command)
        process = subprocess.run(move_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        
        # Optionally handle the output
        if process.returncode == 0:
            print("Move completed successfully.")
        else:
            print(f"Error in repairing PDB: {process.stderr}")

def deltadeltaG(pdb_path, evoef_path,output_directory,mutant,individual_file_path):
    base_filename = os.path.basename(pdb_path).replace('.pdb', '')
    repaired_pdb = os.path.join(f"{output_directory}", f"{base_filename}_Repair.pdb")
    mutant_list = individual_file_path # Path to the file containing mutation specifications
    new_directory = "./EvoEF-mutant-cancer0"
    WT_directory = "./EvoEF-WT-cancer0"
    os.makedirs(new_directory, exist_ok=True)
    os.makedirs(WT_directory, exist_ok=True)
    # Step 2: Build mutants
    build_mutant_command = [
        evoef_path, "--command=BuildMutant", f"--pdb={repaired_pdb}", f"--mutant_file={mutant_list}"
    ]
    print(f"Building mutants for {repaired_pdb}...")
    process = subprocess.run(build_mutant_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if process.returncode == 0:
            print("Build mutant completed successfully.")
    else:
        print(f"Error in mutating PDB: {process.stderr}")
        return None
    mutant_pdb = f"{base_filename}_Repair_Model_0001.pdb"
    new_path = os.path.join(new_directory, f"{base_filename}_{mutant}.pdb")
    shutil.move(mutant_pdb, new_path)  # Move and rename the file
    WT_pdb = f"{base_filename}_Repair_Model_0001_WT.pdb"
    new_path_wt = os.path.join(WT_directory, f"{base_filename}_{mutant}_WT.pdb")
    shutil.move(WT_pdb, new_path_wt)  # Move and rename the file

    # Step 3: Calculate ΔΔGstability
    # Calculate stability for the repaired structure
    # WT_pdb_path = os.path.join(output_directory, base_filename+"_Repair_Model_0001.pdb")
    compute_stability_ref_command = [
        evoef_path, "--command=ComputeStability", f"--pdb={repaired_pdb}"
    ]
    ref_output = subprocess.check_output(compute_stability_ref_command, text=True)
    # print(len(ref_output))
    compute_stability_mut_command = [
                evoef_path, "--command=ComputeStability", f"--pdb={new_path}"
            ]
    mut_output = subprocess.check_output(compute_stability_mut_command, text=True)
    ΔGstability_ref = extract_stability_score(ref_output)
    ΔGstability_mut = extract_stability_score(mut_output)
    ΔΔGstability = ΔGstability_mut - ΔGstability_ref
    print(f"ΔΔGstability for {base_filename}: {ΔΔGstability}")
    return ΔΔGstability


def extract_stability_score(output):
    """Extract the stability score from EvoEF output."""
    for line in output.splitlines():
        if "Total" in line:
            return float(line.split()[-1])
    return None

def main():
    evoef_path = "./EvoEF/EvoEF"
    pdb_directory = "./pdb_cancer"
    output_directory = "./EvoEF-Repair_pdb_cancer0"
    csv_file_path = 'cancer0_part_1.csv'
    os.makedirs(output_directory, exist_ok=True)
    csv_path = "./EvoEF-Repair_pdb_cancer0/results.csv"

    with open(csv_file_path, mode='r', encoding='utf-8') as file:
        csv_reader = csv.reader(file)
        # Convert the reader object to a list (if direct conversion is sufficient)
        data_list = list(csv_reader)
    
    
    # Open the CSV file once, outside the loop
    with open(csv_path, 'w', newline='') as file:
        csv_writer = csv.writer(file)
        csv_writer.writerow(["PDB Name", "Mutant", "ΔΔGstability"])  # Write header

        # Process each entry in data_list (skipping header row)
        individual_file_path = os.path.join(output_directory,'individual.txt')
        for line in data_list[1:]:
            if os.path.isfile(individual_file_path):
                os.remove(individual_file_path)
            with open(individual_file_path, "w+") as file:
                mutant = line[1] + 'A' + line[2] + 'E'
                file.write(mutant + ';\n')

            # Determine PDB path
 
            pdb_p = os.path.join(pdb_directory, line[0])

            pdb_path = pdb_p + '.pdb'

            # Write the computed data to the already open CSV file
            csv_writer.writerow([pdb_p, mutant, process_pdb(pdb_path, evoef_path, output_directory, mutant,individual_file_path)])
            
    # clean_up_command = "mv *.pdb ./EvoEF-repair_pdb"
    # run_command(clean_up_command)
    print("All PDB files have been processed.")

if __name__ == "__main__":
    main()
