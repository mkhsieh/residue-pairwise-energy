# Residue-Pairwise-Energy

This analysis uses atomic contact information of interest to filter residue pairs and provides additional options for backbone-backbone and sidechain-sidechain residual pairwise energies. Currently, the script requires specifying the atoms to be selected (i.e. protein and name CA). The user can define a pairwise cutoff distance and a pairwise filtering threshold (based on the percentage of pairwise contact frequencies within the cutoff distance). For processing, the user can decide the number of CPU cores to use for NAMD2 individual pairwise energy calculations (default: 1 core per task), and the user can also decide the total number of CPU cores will be used in this analysis (default: total number of CPU cores - 1 ).

# Usage: 

pairwiseEng_namd2.py [-h] [--pdb PDB] [--psf PSF] [--traj TRAJ] [--sel1 SEL1 [SEL1 ...]] [--sel2 SEL2 [SEL2 ...]] [--cutoff CUTOFF] [--DistanceMap] [--ContactMap] [--calc] [--intengbackbone] [--intengsidechain] [--exe NAMD2EXE] [--namd2NumCores NAMD2NUMCORES] [--pairfilterthreshold pairfilterthreshold] [--outfolder OUTFOLDER] [--parafolder PARAFOLDER] [--numCores NUMCORES]
