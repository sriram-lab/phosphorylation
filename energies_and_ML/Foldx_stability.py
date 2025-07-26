import csv
import subprocess
import os

base_dir = './cancer_Y_p'
FOLDX_loc="./foldx_20241231"
PDB_loc="./pdb_cancer" #must be given as relative path from FOLDX executable
RepairedPDBs_loc='./merged_cancer_foldx_repair' #must be given as relative path from FOLDX executable
MutantPDBs_loc=os.path.join(base_dir,'Mutant_Structures') #must be given as relative path from FOLDX executable
indivdual_list_loc=os.path.join(base_dir,'individual_list.txt')

# Specify the path to your CSV file
csv_file_path = 'circuit_final1.csv'

def repair_pdb(csv_file_path):
	os.makedirs(RepairedPDBs_loc, exist_ok=True)
	# Open the CSV file and read its contents into a list
	with open(csv_file_path, mode='r', encoding='utf-8') as file:
		csv_reader = csv.reader(file)
        # Convert the reader object to a list (if direct conversion is sufficient)
		data_list = list(csv_reader)
	pdb_list = []
	for data in data_list[1:]:
		pdb_list.append(data[1])
    # uniq_pdbs=['2MH8.pdb']
    # Print the list to verify its contents
	i = 0
	while i<len(pdb_list):
		repaired_pdb = RepairedPDBs_loc+"/"+pdb_list[i]+"_Repair.pdb"
		if os.path.isfile(repaired_pdb):
			print("Already Prepared")
		else:
		#os.mkdir(path)
			try:
				# First command: RepairPDB
				cmd_str = f"{FOLDX_loc} --command=RepairPDB --water=IGNORE --output-dir={RepairedPDBs_loc} --pdb-dir={PDB_loc} --pdb={pdb_list[i]}.pdb"
				print(cmd_str)
				result = subprocess.run(cmd_str, shell=True, capture_output=True)
				result.check_returncode()  # This will raise an exception if the command failed					
				# Second command: Optimize
				cmd_str = f"{FOLDX_loc} --command=Optimize --water=IGNORE --output-dir={RepairedPDBs_loc} --pdb-dir={RepairedPDBs_loc} --pdb={pdb_list[i]}_Repair.pdb"
				print(cmd_str)
				result = subprocess.run(cmd_str, shell=True, capture_output=True)
				result.check_returncode()  # This will raise an exception if the command failed					
				print(f"Finished repairing {pdb_list[i]}")

			except subprocess.CalledProcessError as e:
				print(f"An error occurred while processing {pdb_list[i]}: {e}")
				print(f"Command output: {e.output.decode()}")  # This will display the output of the command which failed
			except Exception as e:
				print(f"An unexpected error occurred: {e}")
		i+=1	
        

        

def buildmodel(csv_file_path):
	os.makedirs(MutantPDBs_loc, exist_ok=True)
	with open(csv_file_path, mode='r', encoding='utf-8') as file:
		csv_reader = csv.reader(file)
		data_list = list(csv_reader)
	#pdb_name = data_list[1][1]
	for data in data_list[727:]:
		if os.path.isfile(indivdual_list_loc):
			os.remove(indivdual_list_loc)
		file = open(indivdual_list_loc,"w+")
		mutations = 'YA'+data[2]+'p'
		file.write('YA'+data[2]+'p;')
		file.close()
		pdb_name = data[1]
		cmd_list = [FOLDX_loc, '--command=BuildModel', '--mutant-file=' + indivdual_list_loc, '--water=IGNORE', '--output-dir=' + MutantPDBs_loc, 
					'--pdb-dir=' + RepairedPDBs_loc, '--pdb=Optimized_' + pdb_name + '_Repair.pdb', '--output-file=' + pdb_name + '_' + mutations,
					'--pdbHydrogens=true', '--pH=7.3', '--ionStrength=0.15', '--out-pdb=true']

		print("Running command:", ' '.join(cmd_list))
		try:
			result = subprocess.run(cmd_list, capture_output=True, text=True)
			# print("STDOUT:", result.stdout)
			# print("STDERR:", result.stderr)
		except Exception as e:
			print("Failed to run command:", e)
			print("STDOUT:", result.stdout)
			print("STDERR:", result.stderr)


def dif_file(csv_file_path):
	with open(csv_file_path, mode='r', encoding='utf-8') as file:
		csv_reader = csv.reader(file)
		data_list = list(csv_reader)
	# Define the name of your CSV file
	filename = os.path.join(base_dir,'cal_fold.csv')

	# Open the file with write ('w') mode
	with open(filename, mode='w', newline='') as file:
		writer = csv.writer(file)
		writer.writerow(['pdb_name', 'mutation','energy'])
		for data in data_list[1:]:
			data2 = mutations = 'YA'+data[2]+'p'
			data1 = pdb_name = data[1]
			fxout = "Dif_"+pdb_name+"_"+mutations+"_Optimized_"+pdb_name+"_Repair.fxout"
			try:
				data3 = extract_deltaG(fxout)
			except:
				data3 = 'error'
			writer.writerow([data1, data2, data3])
	print(f'CSV file "{filename}" created successfully.')

	

def extract_deltaG(file_name):
	file = open(MutantPDBs_loc + '/'+ file_name)
	pdb = False
	total_energy = 0
	for line in file:
		if line.startswith('Pdb'):
			pdb = True
		if pdb :
			values = line.strip().split('\t')
			total_energy = values[1]
	return total_energy

repair_pdb(csv_file_path)
buildmodel(csv_file_path)
dif_file(csv_file_path)