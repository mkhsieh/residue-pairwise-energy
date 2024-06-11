# Residue-Pairwise-Energy

The analysis takes user interested atom contact information to filter residue pairs and provides extra options for backbone and sidechain interresidual pairwise energies. Current, the script requires to specify the atom for selection (i.e., protein and name CA). Users could define the pairwise cutoff distance and pairfilterthreshold (the percentage based on the frequency of pairwise contact within the cutoff distance). For processing, user could decide number of CPU cores for inidividual pairwise energy calculation using NAMD2 (defualtL 1 core per task) and total number of CPU cores will be used in this analysis (defualt: total CPU cores - 1).  


Usage: 

pairwiseEng_namd2.py [-h] [--pdb PDB] [--psf PSF] [--traj TRAJ] [--sel1 SEL1 [SEL1 ...]] [--sel2 SEL2 [SEL2 ...]] [--cutoff CUTOFF] [--DistanceMap] [--ContactMap] [--calc] [--intengbackbone] [--intengsidechain] [--exe NAMD2EXE] [--namd2NumCores NAMD2NUMCORES] [--pairfilterthreshold pairfilterthreshold] [--outfolder OUTFOLDER] [--parafolder PARAFOLDER] [--numCores NUMCORES]
