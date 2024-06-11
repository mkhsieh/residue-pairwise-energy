import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pyemma, itertools
import mdtraj as md
import sys, os, argparse, glob, re
import multiprocessing, subprocess
from multiprocessing import Process

def dis_matrix(args):
    # define residues pairwise maxtrix
    print(f'read pdb file {args.pdb} for topology setup')
    args.file = md.formats.PDBTrajectoryFile(args.pdb, mode='r')
    top=args.file.topology
    ## alpha carbon is utilized for protein residue
    ## N1 atom is utilized for dna base (residue)
    ## TO DO:
    ##   add option for COM of residue
    print(f'the frist selection: {" ".join(args.sel1)}')
    print(f'the second selection: {" ".join(args.sel2)}')
    indices1 = top.select(' '.join(args.sel1))
    indices2 = top.select(' '.join(args.sel2))
    n1 = len(indices1)
    n2 = len(indices2)

    pairs = list(itertools.product(indices1,indices2))
    feat = pyemma.coordinates.featurizer(top)
    feat.add_distances(pairs, periodic=True)

    # import data from traj files
    print(f'read trajectory files: {args.traj}')
    pairdist_traj = args.traj
    pairdist = pyemma.coordinates.load(pairdist_traj, feat)
    pairdist = np.array(pairdist)
    pairdist = pairdist*10
    ## for single traj file
    ## data shape [ n_frame_per_file, n_sel1 x n_sel2 ]
    #  for mutiple traj files
    ## data shape [ n_traj_file, n_frame_per_file, n_sel1 x n_sel2 ]
    #print(pairdist.shape)
    #print(pairdist.ndim)

    if len(args.traj) == 1 : pairdist = data.reshape(pairdist.shape[0],pairdist.shape[-1])
    else : pairdist = pairdist.reshape(pairdist.shape[0]*pairdist.shape[1],pairdist.shape[-1])
    #if pairdist.ndim == 3: pairdist = pairdist.reshape(pairdist.shape[0]*pairdist.shape[1],pairdist.shape[-1])
    #if pairdist.ndim == 2: pairdist = pairdist.reshape(pairdist.shape[0],pairdist.shape[-1])
    pairdist = pairdist.reshape(len(pairdist),n1,n2)

    dis = np.mean(pairdist,axis=0)
    #se = np.std(pairdist, axis=0, ddof=1) / np.sqrt(np.size(pairdist, axis=1))
    frq = np.where(data < args.cutoff, 1, 0)
    #print("frq shape:", frq.shape)
    freq = np.mean(frq,axis=0)
    with open(args.outfolder+'/contact_dis.npy','wb') as f:
        np.save(f,dis)

    with open(args.outfolder+'/contact_freq.npy','wb') as f:
        np.save(f,freq)

    if args.DistanceMap:
        generate_distanacemap(args,dis)

    if args.ContactMap:
        generate_contactmap(args,freq)

    return None

def generate_distanacemap(args,dis):
    xini=1
    yini=1
    plt.figure(figsize=(48,48))
    im = plt.imshow(dis,cmap='Blues', interpolation='nearest',origin='upper')
    ax = plt.gca();
    '''
    # Major ticks
    ax.set_xticks(np.arange(0, n2, 1))
    ax.set_yticks(np.arange(0, n1, 1))

    # Labels for major ticks
    ax.set_xticklabels(np.arange(xini, n2+xini, 1))
    ax.set_yticklabels(np.arange(yini, n1+yini, 1))

    # Minor ticks
    ax.set_xticks(np.arange(-.5, n2, 1), minor=True)
    ax.set_yticks(np.arange(-.5, n1, 1), minor=True)

    '''
    # Gridlines based on minor ticks
    ax.grid(which='minor', color='black', linestyle='-', linewidth=2)
    plt.rcParams["font.weight"]= 'bold'
    plt.rcParams["axes.labelweight"]='bold'
    plt.rcParams["figure.autolayout"]=True
    plt.grid(c='black')
    cbar =plt.colorbar(im)
    cbar.ax.tick_params(labelsize=25)
    plt.savefig(os.path.join(args.outfolder,'/DistanceMap.png'))
    plt.close()

def generate_contactmap(args,freq):
    xini=1
    yini=1
    plt.figure(figsize=(48,48))
    im = plt.imshow(freq,cmap='Greens', interpolation='nearest',origin='upper')
    ax = plt.gca();
    '''
    # Major ticks
    ax.set_xticks(np.arange(0, n2, 1))
    ax.set_yticks(np.arange(0, n1, 1))

    # Labels for major ticks
    ax.set_xticklabels(np.arange(xini, n2+xini, 1))
    ax.set_yticklabels(np.arange(yini, n1+yini, 1))

    # Minor ticks
    ax.set_xticks(np.arange(-.5, n2, 1), minor=True)
    ax.set_yticks(np.arange(-.5, n1, 1), minor=True)

    '''
    # Gridlines based on minor ticks
    ax.grid(which='minor', color='black', linestyle='-', linewidth=2)
    plt.rcParams["font.weight"]= 'bold'
    plt.rcParams["axes.labelweight"]='bold'
    plt.rcParams["figure.autolayout"]=True
    plt.grid(c='black')
    cbar =plt.colorbar(im)
    cbar.ax.tick_params(labelsize=25)
    plt.savefig(os.path.join(args.outfolder,'/ContactMap.png'))
    plt.close()

def pairwise_energy_calc(args,freq):
    top=args.file.topology
    natom=top.n_atoms
    # residue check
    #for res in top.residues:
    #    print(res,res.name,res.index)

    ### filtering residue pairs for energy calc

    indices1 = top.select(' '.join(args.sel1))
    indices2 = top.select(' '.join(args.sel2))
    pairs = np.argwhere(freq > args.pairfilterpercentage/100.0 )
    args.filtered_pairs = []

    if len(pairs) > 0:
        for a, b in pairs:
            if a < b:
                #print('pair mextri index: ',a,b)
                #print('indexes: ',indices1[a],indices2[b])
                for atom in top.atoms:
                    if atom.index == indices1[a]:
                        a_res =atom.residue
                        a_resid =atom.residue.index # residue.index starts from 1
                        #print(atom.name, atom.index, atom.residue, atom.residue.index)
                    if atom.index == indices2[b]:
                        b_res =atom.residue
                        b_resid =atom.residue.index
                        #print(atom.name, atom.index, atom.residue, atom.residue.index)
                args.filtered_pairs.append([a_resid, b_resid])
                #print('resid: ',a_resid,b_resid)
        print(args.filtered_pairs)
        if len(args.filtered_pairs) == 0:
            print(f'there is no residue pairs over the pairfilterpercentage {args.pairfilterpercentage} %%.')
        else:
            calcEnergiesNAMD(args)
            collectEnergiesNAMD(args)
    else:
        print(f'there is no residue pairs over the pairfilterpercentage {args.pairfilterpercentage} %%.')


def collectEnergiesNAMD(args):
    filePaths = glob.glob(args.outfolder+'/*.log')
    #print(filePaths)
    energiesDict = dict()
    for filePath in filePaths:
        # Get the interaction residues
        matches = re.search('(\d+)_(\d+)_energies.log',filePath)
        if not matches:
            continue

        # Get residue indices
        res1 = int(matches.groups()[0])
        res2 = int(matches.groups()[1])

        # Read in the first line (header) output file and count number of total lines.
        f = open(filePath,'r')
        lines = f.readlines()

        # Ignore lines not starting with ENERGY:
        lines = [line for line in lines if line.startswith('ENERGY:')]
        f.close()

        lines = [line.strip('\n').split() for line in lines if line.startswith('ENERGY:')]
        lines = [[float(integer) for integer in line[1:]] for line in lines]

        # Frame: 0, Elec: 5, VdW: 6, Total: 10
        headers = ['Frame','Elec','VdW','Total']
        headerColumns = [0,5,6,10] # Order in which the headers appear in NAMD2 log

        numTitles = len(headers)
        # Assign each column into appropriate key's value in energyOutput dict.
        energyOutput = dict()
        for i in range(0,numTitles):
            energyOutput[headers[i]] = [line[headerColumns[i]] for line in lines]

        # Puts this energyOutput dict into energies dict with keys as residue ids
        energiesDict[str(res1)+'-'+str(res2)] = energyOutput
        # Also store it as res2,res1 (it is the same thing after all)
        energiesDict[str(res2)+'-'+str(res1)] = energyOutput
        #print(energiesDict)

    # Prepare a pandas data table from parsed energies, write it to new files depending on type of energy
    df_total = pd.DataFrame()
    df_elec = pd.DataFrame()
    df_vdw = pd.DataFrame()
    for key,value in list(energiesDict.items()):
        # Try-check block to allow collecting even when parse of pair IE fails.
        try:
            df_total[key] = value['Total']
            df_elec[key] = value['Elec']
            df_vdw[key] = value['VdW']
        except:
            print('Failed to parse interaction data for pair '+str(key))

    df_total.to_csv(os.path.join(args.outfolder,'energies_intEnTotal.csv'))
    df_elec.to_csv(os.path.join(args.outfolder,'energies_intEnElec.csv'))
    df_vdw.to_csv(os.path.join(args.outfolder,'energies_intEnVdW.csv'))

    return None

def calcEnergiesNAMD(args):
    # use Process for multi-task distributing to CPU cores
    pairs=np.array_split(args.filtered_pairs,args.numCores)
    procs = []
    for pairs_singleCore in pairs:
        proc = Process(target=run_single_namd2,args=(args,pairs_singleCore))
        procs.append(proc)
        proc.start()

    # complete the processes
    for proc in procs:
        proc.join()

def run_single_namd2(args, filtered_pairs_singleCore):
    if not os.path.exists(args.outfolder):
        os.makedirs(args.outfolder)

    #prepare_tmp_pdb
    natom=args.file.topology.n_atoms
    args.filtered_pairs_singleCore = filtered_pairs_singleCore
    for a_resid, b_resid in args.filtered_pairs_singleCore:
        bfactor=bfactor_tmp_pdb(args,a_resid,b_resid)
        pairPDB = '%s/%s_%s-tmp.pdb' % (args.outfolder,a_resid,b_resid)
        outf = md.Trajectory(args.file.positions, args.file.topology)
        outf.save_pdb(pairPDB,bfactors=bfactor)

        write_namdconf(args,a_resid,b_resid)
        namdConf = '%s/%s_%s-tmp.namd' % (args.outfolder,a_resid,b_resid)

        # Run namd2 to compute the energies
        try:
            pid_namd2 = subprocess.Popen([args.namd2exe,'+p%i' % args.namd2NumCores,namdConf],
                stdout=open(
                    os.path.join(args.outfolder,'%i_%i_energies.log' %
                        (a_resid,b_resid)),'w'),stderr=subprocess.PIPE)
            _,error = pid_namd2.communicate()

        except KeyboardInterrupt:
            print('Keyboard interrupt detected.')
            sys.exit(0)

        if error:
            print('Error while calling NAMD executable:\n'+str(error))
            error = error.decode().split('\n')
            fatalErrorLine = None

            for i in range(0,len(error)):
                if 'FATAL ERROR:' in error[i]:
                    fatalErrorLine = error[i]
                    continue

            if fatalErrorLine:
                return fatalErrorLine

        pid_namd2.wait()

def bfactor_tmp_pdb(args,a_resid,b_resid):
    natom=args.file.topology.n_atoms
    bfactor = np.zeros(natom)
    if args.intengbackbone:
        for atom in args.file.topology:
            if atom.residue.index == a_resid:
                for atom.residue.index in args.file.topology.select("protein and backbone"): bfactor[atom.index] = 1
            if atom.residue.index == b_resid:
                for atom.residue.index in args.file.topology.select("protein and backbone"): bfactor[atom.index] = 2
                    
    elif args.intengsidechain:
        for atom in args.file.topology:
            if atom.residue.index == a_resid:
                for atom.residue.index in args.file.topology.select("protein and sidechain"): bfactor[atom.index] = 1
            if atom.residue.index == b_resid:
                for atom.residue.index in args.file.topology.select("protein and sidechain"): bfactor[atom.index] = 2
    else:
        for atom in args.file.topology.atoms:
            if atom.residue.index == a_resid:
               bfactor[atom.index] = 1
            if atom.residue.index == b_resid:
                bfactor[atom.index] = 2
    return bfactor

def prepare_tmp_pdb(args,file,residue_pairs):
    for a_resid, b_resid in residue_pairs:
        bfactor=_bfactor_tmp_pdb(file.topology,a_resid,b_resid)
        pairPDB = '%s/%s_%s-tmp.pdb' % (args.outfolder,a_resid,b_resid)
        outf = md.Trajectory(file.positions, file.topology)
        outf.save_pdb(pairPDB,bfactors=bfactor)

def write_namdconf(args,a_resid,b_resid):
    cutoff = 12.0
    switchdist = 10.0
    pairPDB = '%s/%s_%s-tmp.pdb' % (args.outfolder,a_resid,b_resid)
    namdConf = '%s/%s_%s-tmp.namd' % (args.outfolder,a_resid,b_resid)
    f = open(namdConf,'w')
    f.write('structure %s\n' % (os.path.dirname(__file__)+'/'+args.psf))
    f.write('paraTypeCharmm on\n')
    f.write('parameters %s\n' % (os.path.dirname(__file__)+'/'+args.parafolder+'/par_all36m_prot.prm'))
    f.write('parameters %s\n' % (os.path.dirname(__file__)+'/'+args.parafolder+'/par_all36_na.prm'))
    f.write('parameters %s\n' % (os.path.dirname(__file__)+'/'+args.parafolder+'/par_all36_carb.prm'))
    f.write('parameters %s\n' % (os.path.dirname(__file__)+'/'+args.parafolder+'/par_all36_lipid.prm'))
    f.write('parameters %s\n' % (os.path.dirname(__file__)+'/'+args.parafolder+'/par_all36_cgenff.prm'))
    f.write('parameters %s\n' % (os.path.dirname(__file__)+'/'+args.parafolder+'/par_interface.prm'))
    f.write('parameters %s\n' % (os.path.dirname(__file__)+'/'+args.parafolder+'/toppar_all36_moreions.str'))
    f.write('parameters %s\n' % (os.path.dirname(__file__)+'/'+args.parafolder+'/toppar_all36_nano_lig.str'))
    f.write('parameters %s\n' % (os.path.dirname(__file__)+'/'+args.parafolder+'/toppar_all36_nano_lig_patch.str'))
    f.write('parameters %s\n' % (os.path.dirname(__file__)+'/'+args.parafolder+'/toppar_all36_synthetic_polymer.str'))
    f.write('parameters %s\n' % (os.path.dirname(__file__)+'/'+args.parafolder+'/toppar_all36_synthetic_polymer_patch.str'))
    f.write('parameters %s\n' % (os.path.dirname(__file__)+'/'+args.parafolder+'/toppar_all36_polymer_solvent.str'))
    f.write('parameters %s\n' % (os.path.dirname(__file__)+'/'+args.parafolder+'/toppar_water_ions.str'))
    f.write('parameters %s\n' % (os.path.dirname(__file__)+'/'+args.parafolder+'/toppar_dum_noble_gases.str'))
    f.write('parameters %s\n' % (os.path.dirname(__file__)+'/'+args.parafolder+'/toppar_ions_won.str'))
    f.write('parameters %s\n' % (os.path.dirname(__file__)+'/'+args.parafolder+'/toppar_all36_prot_arg0.str'))
    f.write('parameters %s\n' % (os.path.dirname(__file__)+'/'+args.parafolder+'/toppar_all36_prot_c36m_d_aminoacids.str'))
    f.write('parameters %s\n' % (os.path.dirname(__file__)+'/'+args.parafolder+'/toppar_all36_prot_fluoro_alkanes.str'))
    f.write('parameters %s\n' % (os.path.dirname(__file__)+'/'+args.parafolder+'/toppar_all36_prot_heme.str'))
    f.write('parameters %s\n' % (os.path.dirname(__file__)+'/'+args.parafolder+'/toppar_all36_prot_na_combined.str'))
    f.write('parameters %s\n' % (os.path.dirname(__file__)+'/'+args.parafolder+'/toppar_all36_prot_retinol.str'))
    f.write('parameters %s\n' % (os.path.dirname(__file__)+'/'+args.parafolder+'/toppar_all36_prot_model.str'))
    f.write('parameters %s\n' % (os.path.dirname(__file__)+'/'+args.parafolder+'/toppar_all36_prot_modify_res.str'))
    f.write('parameters %s\n' % (os.path.dirname(__file__)+'/'+args.parafolder+'/toppar_all36_na_nad_ppi.str'))
    f.write('parameters %s\n' % (os.path.dirname(__file__)+'/'+args.parafolder+'/toppar_all36_na_rna_modified.str'))
    f.write('parameters %s\n' % (os.path.dirname(__file__)+'/'+args.parafolder+'/toppar_all36_lipid_sphingo.str'))
    f.write('parameters %s\n' % (os.path.dirname(__file__)+'/'+args.parafolder+'/toppar_all36_lipid_archaeal.str'))
    f.write('parameters %s\n' % (os.path.dirname(__file__)+'/'+args.parafolder+'/toppar_all36_lipid_bacterial.str'))
    f.write('parameters %s\n' % (os.path.dirname(__file__)+'/'+args.parafolder+'/toppar_all36_lipid_cardiolipin.str'))
    f.write('parameters %s\n' % (os.path.dirname(__file__)+'/'+args.parafolder+'/toppar_all36_lipid_cholesterol.str'))
    f.write('parameters %s\n' % (os.path.dirname(__file__)+'/'+args.parafolder+'/toppar_all36_lipid_dag.str'))
    f.write('parameters %s\n' % (os.path.dirname(__file__)+'/'+args.parafolder+'/toppar_all36_lipid_inositol.str'))
    f.write('parameters %s\n' % (os.path.dirname(__file__)+'/'+args.parafolder+'/toppar_all36_lipid_lnp.str'))
    f.write('parameters %s\n' % (os.path.dirname(__file__)+'/'+args.parafolder+'/toppar_all36_lipid_lps.str'))
    f.write('parameters %s\n' % (os.path.dirname(__file__)+'/'+args.parafolder+'/toppar_all36_lipid_mycobacterial.str'))
    f.write('parameters %s\n' % (os.path.dirname(__file__)+'/'+args.parafolder+'/toppar_all36_lipid_miscellaneous.str'))
    f.write('parameters %s\n' % (os.path.dirname(__file__)+'/'+args.parafolder+'/toppar_all36_lipid_model.str'))
    f.write('parameters %s\n' % (os.path.dirname(__file__)+'/'+args.parafolder+'/toppar_all36_lipid_prot.str'))
    f.write('parameters %s\n' % (os.path.dirname(__file__)+'/'+args.parafolder+'/toppar_all36_lipid_tag.str'))
    f.write('parameters %s\n' % (os.path.dirname(__file__)+'/'+args.parafolder+'/toppar_all36_lipid_yeast.str'))
    f.write('parameters %s\n' % (os.path.dirname(__file__)+'/'+args.parafolder+'/toppar_all36_lipid_hmmm.str'))
    f.write('parameters %s\n' % (os.path.dirname(__file__)+'/'+args.parafolder+'/toppar_all36_lipid_detergent.str'))
    f.write('parameters %s\n' % (os.path.dirname(__file__)+'/'+args.parafolder+'/toppar_all36_lipid_ether.str'))
    f.write('parameters %s\n' % (os.path.dirname(__file__)+'/'+args.parafolder+'/toppar_all36_carb_glycolipid.str'))
    f.write('parameters %s\n' % (os.path.dirname(__file__)+'/'+args.parafolder+'/toppar_all36_carb_glycopeptide.str'))
    f.write('parameters %s\n' % (os.path.dirname(__file__)+'/'+args.parafolder+'/toppar_all36_carb_imlab.str'))
    f.write('parameters %s\n' % (os.path.dirname(__file__)+'/'+args.parafolder+'/toppar_all36_label_spin.str'))
    f.write('parameters %s\n' % (os.path.dirname(__file__)+'/'+args.parafolder+'/toppar_all36_label_fluorophore.str'))
    f.write('numsteps 1\n')
    f.write('exclude scaled1-4\n')
    f.write('outputname %i_%i-tmp\n' % (a_resid, b_resid))
    f.write('temperature 0\n')
    f.write('COMmotion yes\n')
    f.write('cutoff %f\n' % cutoff)
    if switchdist:
        f.write('switching on\n')
        f.write('switchdist %f\n' % switchdist)
    else:
        f.write('switching off\n')

    f.write('pairInteraction on\n')
    f.write('pairInteractionFile %s\n' % (os.path.dirname(__file__)+'/'+pairPDB))
    f.write('pairInteractionGroup1 1\n')
    f.write('pairInteractionGroup2 2\n')
    f.write('coordinates %s\n' % (os.path.dirname(__file__)+'/'+pairPDB))
    f.write('set ts 0\n')
    #f.write('coorfile open dcd %i_%i-temp.dcd\n' % (pair[0],pair[1]))
    f.write('coorfile open dcd %s\n' % (os.path.dirname(__file__)+'/'+args.traj))
    f.write('while { ![coorfile read] } {\n')
    f.write('\tfirstTimeStep $ts\n')
    f.write('\trun 0\n')
    f.write('\tincr ts 1\n')
    #f.write('\tcoorfile skip\n') # Don't need it once you apply stride to tray_dry.dcd in outputfolder
    f.write('}\n')
    f.write('coorfile close')

def prepare_namdconf(args,residue_pairs):
    for a_resid, b_resid in residue_pairs:
        _write_namdconf(args,a_resid,b_resid)

def main():
    parser = argparse.ArgumentParser(
            description='calc contact map.'
    )
    parser.add_argument('--pdb', dest='pdb', type=str, default=[False],
                        help='the corresponding PDB file of the DCD trajectory. '
    )
    parser.add_argument('--psf', dest='psf', type=str, default=[False],
                        help='the corresponding PSF file of the DCD trajectory. '
    )
    #parser.add_argument('--traj', dest='traj', type=str, default=[False], nargs='+',
    parser.add_argument('--traj', dest='traj', type=str, default=[False],
                        help='Name of the trajectory file'
    )
    parser.add_argument('--sel1', dest='sel1', type=str, default=['all'], nargs='+',
                        help='A atom selection (mdtraj selection style) string which determines the first group of selected '
                             'residues. '
    )
    parser.add_argument('--sel2', dest='sel2', type=str, default=['all'], nargs='+',
                        help='A atom selection (mdtraj selection style) string which determines the second group of selected '
                             'residues.'
    )
    parser.add_argument('--cutoff', dest='cutoff', type=float, default=12, nargs=1,
                        help='Non-bonded interaction distance cutoff (Angstroms) for pairfilter which generates '
                             'frequency matrix of contacts.'
    )
    parser.add_argument(
            '--DistanceMap', dest='DistanceMap', action='store_true',
            help='generation of residue pairwise (mean) distance map. ',
            default=False,
    )
    parser.add_argument(
            '--ContactMap', dest='ContactMap', action='store_true',
            help='generation of residue contact frequency map with . ',
            default=False,
    )
    parser.add_argument(
            '--calc', dest='calc', action='store_true',
            help='performing residue pairwise energy calculation. ',
            default=False,
    )
    parser.add_argument(
            '--intengbackbone', dest='intengbackbone', action='store_true',
            help='performing residue pairwise backbone energy calculation. ',
            default=False,
    )
    parser.add_argument(
            '--intengsidechain', dest='intengsidechain', action='store_true',
            help='performing residue pairwise sidechain energy calculation. ',
            default=False,
    )
    parser.add_argument('--exe', dest='namd2exe', default=None, type=str,
                        help='Path to the namd2 executable. (assumes namd2 is in the executable search path.)'
    )
    parser.add_argument('--namd2NumCores', dest='namd2NumCores', default=1, type=int,
                        help='Number of CPU cores to be employed for interaction energy calculation '
                             'by NAMD2 executable in a single subprocess. If not specified, it defaults to 1. '
    )
    parser.add_argument('--pairfilterthreshold', dest='pairfilterthreshold', type=float, default=75,
                        help='When given, residue pairs whose contact frequency are over pairfilterthreshold will be taken '
                             'into account in further evaluations. When not given, it defaults to 75%%'
    )
    parser.add_argument('--outfolder', dest='outfolder', type=str, default='outfolder',
                        help='Folder path for storing calculation results. '
    )
    parser.add_argument('--parafolder', dest='parafolder', type=str, default='toppar',
                        help='Folder path for storing calculation results. '
    )
    parser.add_argument('--numCores', dest='numCores', default=multiprocessing.cpu_count(), type=int,
                        help='Number of CPU cores to be employed. '
                             'If not specified, it defaults to the number of cpu cores present '
                             'in your computer.'
    )

    args = parser.parse_args()
    if not os.path.exists(args.outfolder):
        os.makedirs(args.outfolder)

    dis, freq = dis_matrix(args)
    if args.calc:
        pairwise_energy_calc(args,freq)

if __name__ == '__main__':
    main()
