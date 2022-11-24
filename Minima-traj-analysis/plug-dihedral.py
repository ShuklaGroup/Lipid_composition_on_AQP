import mdtraj as md
import numpy as np
import itertools
import math
import pickle
import matplotlib.pyplot as plt
import glob
import os

#get the indeces array of seletected atoms
file_phi = glob.glob('plug-phi-atoms.npy')[0]
phi = np.load(file_phi)
assert phi.shape == (1,4)

#get topology file in directory
topname = glob.glob('*_wat_lip.p*m*')[0]
print("Topology file is:" + topname)

#go through all stripped traj (containing protein and water) in directory, calculate dihedral for each
for file in sorted(glob.glob('*MD*_wat_lip.xtc')):
    name = file.replace('.xtc','')
    print(name)
    phi_outfile = name + '-plug-PHI.npy'
    a = os.path.isfile(phi_outfile)
    if a == False:
        traj = md.load(file, top = topname)
        top = traj.topology
        phi_array = md.compute_dihedrals(traj, phi) #shape=(n_frames, num_pairs)
        #making sure final array is of the correct shape
        assert len(phi_array) == len(traj)
        assert len(phi_array[0]) == len(phi)
        np.save(phi_outfile, phi_array)
