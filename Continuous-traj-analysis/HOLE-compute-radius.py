import MDAnalysis as mda
import numpy as np
import itertools
import math
import pickle
import matplotlib.pyplot as plt
import glob
import os
import mdtraj as md
from MDAnalysis.analysis import hole2
import time

def ComputeHOLE(file):
    '''
    Compute HOLE radius for each .xtc trajectory file
    Args:
        file (str): file name of a trajectory in current directory
    Returns:
        gather (dict) saved as a .pkl file: dictionary output from HOLE MDAnalysis
            containing the relative z-coord, radius at each coordinate, and cen_line_D
            keys: "z-coord", "radius", "cen_line_D"
            values: np.array of shape (num of frames, num of computed radii)
    '''
    name = file.replace('.xtc','')
    print(name)
    dist_file = name + '-HOLE-DICT.pkl'
    a = os.path.isfile(dist_file)
    if a == False:
        #get cvect between SER35-VAL151 (0-indexed)
        topname_parm = glob.glob('*_wat_lip.p*m*')[0]
        print("Topology file is:" + topname_parm)
        traj = md.load(file, top = topname_parm)
        top = traj.topology
        cvect = traj.xyz[0,2284] - traj.xyz[0,581]
        print('cvect calculated')

        #do HOLE
        topname_pdb = glob.glob('*_wat_lip.pdb')[0]
        u = mda.Universe(topname_pdb, file)
        
        localtime = time.asctime(time.localtime(time.time()))
        print('starting HOLE', localtime)

        #replace executable with the path to HOLE2 program
        with hole2.HoleAnalysis(u, executable='/PATH/hole2/exe/hole', cvect=cvect, end_radius=4) as h2:
            h2.run()
            gather = h2.gather()
            h2.delete_temporary_files()
        
        localtime = time.asctime(time.localtime(time.time()))
        print('finish HOLE', localtime)        
        gather_file = open(name + '-HOLE-DICT.pkl', 'wb')
        pickle.dump(gather, gather_file)
        gather_file.close()

if __name__ == '__main__':
    for file in sorted(glob.glob('*MD*_wat_lip.xtc')):
        ComputeHOLE(file)
