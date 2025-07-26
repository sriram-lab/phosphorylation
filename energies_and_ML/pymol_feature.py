from os.path import isfile
import pymol.cmd as cmd
import csv



# Open the output CSV file where results will be written
with open('output_features_tsu(all).csv', 'w', newline='') as output_file:
    writer = csv.writer(output_file)
    # Write headers for the CSV file
    writer.writerow(['Area', 'Number of Residues', 'Number of Atoms', 'Number of COOH Atoms', 'pdb', 'Residue Number', 'ss', 'phi', 'psi', 'Length of PDB', 'Mut Res/Length','Exp','Mut','EvoEF','Foldx'])

    # Load the input CSV file
    with open('phoslist_tsu_new.csv', "r") as file:
        lines = file.readlines()

    for line in lines[1:]:
        data = line.strip().split(',')
        pdbfile = "pdb_files/" + data[1] # Assuming the second column is the PDB file name
        resinum = int(float(data[4]))  # Assuming the third column is the residue number

        if isfile(pdbfile):
            print(f"Processing {pdbfile} for residue {resinum}")
            # cmd.reinitialize()
            cmd.load(pdbfile, "prot")
            cmd.set('dot_solvent', 'on')
            
            # Calculate the total number of residues in the PDB
            total_residues = cmd.count_atoms('name CA')  # Count alpha carbons, assumes one per residue
            
            # Selecting the residue of interest
            try:
                cmd.select('sel', f'resi {resinum} and chain A')
            except:
                cmd.select('sel', f'resi {resinum}')

            # Area calculation
            area = cmd.get_area('sel')
            cmd.select('sel2', f'all within 5 of sel')
            natom = cmd.count_atoms('sel2')
            cmd.select('sel2', 'byres sel2 and not sel')
            nres = cmd.count_atoms('sel2 and name CA')
            ncooh = cmd.count_atoms('sel2 and name CA and (resname GLU or resname ASP)')

            # Secondary structure calculation
            cmd.select('helix', 'sel and ss h')
            cmd.select('sheet', 'sel and ss s')
            ss = 0
            if cmd.count_atoms('helix') > 0:
                ss = 1
            elif cmd.count_atoms('sheet') > 0:
                ss = 2

            # Compute the dihedral angles phi and psi
            phi = None
            psi = None
            chains = cmd.get_chains("all")
            try:
                if len(chains) > 1:
                    chain = 'A'  # Assuming chain A, adjust as needed
                    if cmd.count_atoms(f'resi {resinum} and chain {chain}') > 0:
                        atom1 = f"name C and resi {resinum-1} and chain {chain}"
                        atom2 = f"name N and resi {resinum} and chain {chain}"
                        atom3 = f"name CA and resi {resinum} and chain {chain}"
                        atom4 = f"name C and resi {resinum} and chain {chain}"
                        phi = cmd.get_dihedral(atom1, atom2, atom3, atom4)
                    
                    if cmd.count_atoms(f'resi {resinum+1} and chain {chain}') > 0:
                        atom1 = f"name N and resi {resinum} and chain {chain}"
                        atom2 = f"name CA and resi {resinum} and chain {chain}"
                        atom3 = f"name C and resi {resinum} and chain {chain}"
                        atom4 = f"name N and resi {resinum+1} and chain {chain}"
                        psi = cmd.get_dihedral(atom1, atom2, atom3, atom4)
                else:
                    if cmd.count_atoms(f'resi {resinum}') > 0:
                        atom1 = f"name C and resi {resinum-1}"
                        atom2 = f"name N and resi {resinum}"
                        atom3 = f"name CA and resi {resinum}"
                        atom4 = f"name C and resi {resinum}"
                        phi = cmd.get_dihedral(atom1, atom2, atom3, atom4)
                    
                    if cmd.count_atoms(f'resi {resinum+1}') > 0:
                        atom1 = f"name N and resi {resinum}"
                        atom2 = f"name CA and resi {resinum}"
                        atom3 = f"name C and resi {resinum}"
                        atom4 = f"name N and resi {resinum+1}"
                        psi = cmd.get_dihedral(atom1, atom2, atom3, atom4)
            except:
                print(f"Error calculating dihedrals for residue {resinum}")

            # Calculate the mutation residue ratio
            mut_res_length_ratio = resinum / total_residues if total_residues else 0

            # Write the computed data to the CSV file
            writer.writerow([area, nres, natom, ncooh, data[1], resinum, ss, phi, psi, total_residues, mut_res_length_ratio,data[0],data[2],data[6],data[7]])
            cmd.delete('all')
        else:
            print(f"File not found: {pdbfile}")
            writer.writerow([''] * 13)  # Write empty fields if file not found

# print("All data processed and saved to output_features_tsu.csv.")

