# residue-pairwise-energy

usage: pairwiseEng_namd2.py [-h] [--pdb PDB] [--psf PSF] [--traj TRAJ] [--sel1 SEL1 [SEL1 ...]] [--sel2 SEL2 [SEL2 ...]] [--cutoff CUTOFF] [--DistanceMap] [--ContactMap] [--calc] [--exe NAMD2EXE] [--namd2NumCores NAMD2NUMCORES] [--pairfilterpercentage PAIRFILTERPERCENTAGE] [--outfolder OUTFOLDER]
[--parafolder PARAFOLDER] [--numCores NUMCORES]

options:
-h, --help show this help message and exit
--pdb PDB the corresponding PDB file of the DCD trajectory.
--psf PSF the corresponding PSF file of the DCD trajectory.
--traj TRAJ Name of the trajectory file
--sel1 SEL1 [SEL1 ...]
A atom selection (mdtraj selection style) string which determines the first group of selected residues.
--sel2 SEL2 [SEL2 ...]
A atom selection (mdtraj selection style) string which determines the second group of selected residues.
--cutoff CUTOFF Non-bonded interaction distance cutoff (Angstroms) for pairfilter which generates frequency matrix of contacts.
--DistanceMap generation of residue pairwise (mean) distance map.
--ContactMap generation of residue contact frequency map with .
--calc performing residue pairwise energy calculation.
--exe NAMD2EXE Path to the namd2 executable. (assumes namd2 is in the executable search path.)
--namd2NumCores NAMD2NUMCORES
Number of CPU cores to be employed for interaction energy calculation by NAMD2 executable in a single subprocess. If not specified, it defaults to 1, and NUMCORES subprocesses are used.
--pairfilterpercentage PAIRFILTERPERCENTAGE
When given, residues whose centers of masses are within the PAIRFILTERCUTOFF distance from each other for at least PAIRFILTERPERCENTAGE percent of the trajectory will be taken into account in further evaluations. When not given, it defaults to 75%
--outfolder OUTFOLDER
Folder path for storing calculation results.
--parafolder PARAFOLDER
Folder path for storing calculation results.
--numCores NUMCORES Number of CPU cores to be employed. If not specified, it defaults to the number of cpu cores present in your computer.
