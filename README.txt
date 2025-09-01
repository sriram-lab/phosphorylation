PhosphoDDG: A program to predict protein stability change due to phosphorylation
For the most recent version, see Github: https://github.com/sriram-lab/phosphorylation

Pre-Machine Learning
--------------------

prepare_tsu1.m
	Creates a list of PDB IDs, energies, and residues for the Tsuboyama data
	input: dG_non_redundant_natural_Fig5.csv
	output: phoslist_tsu.csv

prepare_tsu2.m
	Aligns sequences for Tsuboyama proteins and adds new residue numbers to list
	input: phoslist_tsu.csv
	output: phoslist_tsu_new.csv

prepare_pancancer1.m
	Creates a list of Uniprot IDs and residues for the cancer dataset
	input: mmc1.xlsx, mmc2.xlsx
 	output: phostable.csv

prepare_pancancer2.m
	Reads alignments for cancer proteins and adds new residue numbers to list
	input: phostable.csv, aligned sequences
	output: phoslist_with_conversion.csv

catboost_features.py
	Calculates FoldX and structural features
	Usage: python catboost_features.py tsu to calculate the features for the Tsuboyama dataset
	       python catboost_features.py cancer to calculate the features for the cancer dataset
		python your_script_name.py single --pdb_file <pdbfile> --mutation <ex. "S123A"> --output_csv <outfile> 
Machine Learning
--------------------
CatBoostRegression.py
	Reads the feature file from the external and Tsuiboyama datasets and constructs the cataboost model
	input: external and tsuiboyama datasets as csv
	ouutput: Catboost model 

Post-Machine Learning
---------------------
remove_duplicates.m
	Removes rows with the same uniprot ID and mutated residue number
	input: Dataset_6.csv (from machine learning)
	output: edited Dataset_6.csv

centralities.m
	Determines whether network locations with certain centrality measures are prone to harbor destabilizing phosphorylations
	input: Dataset_6.csv, mmc2.xlsx
 	output: box plots and p values for centralities with low and high ddG

search_humsavar.m
	Searches the humsavar database for tyrosine to aspartate mutations at the same residue positions as phosphorylations from the cancer dataset
	input: humsavar2.txt, Dataset_6.csv
	output: rows from cancer output data matching to humsavar

search_clinvar.m
	Searches the clinvar database for tyrosine to aspartate mutations at the same residue positions as phosphorylations from the cancer dataset
	input: variant_summary.txt, Dataset_6.csv
	output: rows from cancer output data matching to clinvar

find_cancer_rapid.m 
	Search for energy of mutation within all alphafold2 mutations, for the cancer ptm dataset
	input: alpha_pedictions.csv, Dataset_6.csv
	output: column of energies

psp_search.m
	Search phophosite plus for energy of mutation within alphafold2 mutations
	input: alpha_predictions.csv, posit
	output: energies

The list of all 300k tyrosine predictions from alphafold2 are available on Zenodo. The XGB_Predictions column contains predictions referenced in the publication. 

---------------------
The table column headings in Tables S5 and S6 are as follows:
uniprot: The UniProt ID
Mut_res: phosphorylated/phosphomimetic mutated residue 
ML predictions: Delta-Delta-G prediction of the Catboost full method (kcal/mol)
NCBI_Gene_ID: NCBI Gene ID
role_cancer: Tumor suppressor, oncogene, or driver annotations from CancerMine
gene_hugo_id: HGNC gene ID
phosphorylation FoldX: Delta-Delta-G of the FoldX direct phosphorylation model (kcal/mol)
Number of Residues: Number of residues with any atom within 5 Angstroms of the phosphorylated residue
Number of Atoms: Number of atoms within 5 Angstroms
Number of COOH atoms: Number of GLU or ASP side chain oxygens within 5 Angstroms
phi: residue phi angle (degrees)
psi: residue psi angle (degrees)
Length of PDB: Number of residues in protein
FoldX: FoldX phosphomimetic Delta-Delta-G (kcal/mol)
SS: Secondary structure (DSSP)
RelSASA: Relative Solvent Accessible Surface Area (DSSP)
AltPosition: residue’s location within its secondary structure element (see publication Methods)
Hairpin: presence within a beta hairpin
The rest of the contact terms quantify short-range or long-range contacts within 10 Angstroms, of various types
SASA: Absolute solvent accessible surface area
P: Number of local (residue-based) parallel relations
X: Number of local cross relations 
IP: Number of local inverse parallel relations
S: Number of local series relations

 