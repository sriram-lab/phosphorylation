# Prediction of ddG for Tyrosine phosphomimetic mutation of protein 

We used FoldX, EvoEF and structural information to predict the ddG of phosphormimetic protein. We trained GLM(generalized linear model), GBM(Gradient boosted machine) and XGboost on the tsuboyama et al. dataset with the specific data of Tyrosine mutating to Glutamine. 

## Programs

- Use FoldX to calculate ddG given a csv file of protein and mutations
- Use EvoEF to calculate ddG given a csv file of protein and mutations
- Use Pymol to extract protein residue features like solvent accessible surface are, number of residue in contact etc.
- Use Various method to predict the ddG of mutations

## Prerequisites

Before you begin, ensure you have met the following requirements:
- FoldX
- EvoEF
- Pymol


