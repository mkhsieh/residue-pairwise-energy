# Residue-Pairwise-Energy

The analysis takes user interested atom contact information to filter residue pairs and provides extra options for backbone and sidechain interresidual pairwise energies. Current, the script requires to specify the atom for selection (i.e., protein and name CA). Users could define the pairwise cutoff distance and pairfilterthreshold (the percentage based on the frequency of pairwise contact within the cutoff distance). For processing, user could decide number of CPU cores for inidividual pairwise energy calculation using NAMD2 (defualtL 1 core per task) and total number of CPU cores will be used in this analysis (defualt: total CPU cores - 1).  


Usage: 

pairwiseEng_namd2.py [-h] [--pdb PDB] [--psf PSF] [--traj TRAJ] [--sel1 SEL1 [SEL1 ...]] [--sel2 SEL2 [SEL2 ...]] [--cutoff CUTOFF] [--DistanceMap] [--ContactMap] [--calc] [--intengbackbone] [--intengsidechain] [--exe NAMD2EXE] [--namd2NumCores NAMD2NUMCORES] [--pairfilterthreshold pairfilterthreshold] [--outfolder OUTFOLDER] [--parafolder PARAFOLDER] [--numCores NUMCORES]


User
--namd2NumCores NAMD2NUMCORES
Number of CPU cores to be employed for interaction energy calculation by NAMD2 executable in a single subprocess. If not specified, it defaults to 1, and NUMCORES subprocesses are used.

--pairfilterpercentage PAIRFILTERPERCENTAGE
When given, residues whose centers of masses are within the PAIRFILTERCUTOFF distance from each other for at least PAIRFILTERPERCENTAGE percent of the trajectory will be taken into account in further evaluations. When not given, it defaults to 75%

--outfolder OUTFOLDER
Folder path for storing calculation results.

--parafolder PARAFOLDER
Folder path for storing calculation results.

--numCores NUMCORES Number of CPU cores to be employed. If not specified, it defaults to the number of cpu cores present in your computer.
