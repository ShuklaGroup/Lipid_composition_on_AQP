import numpy as np
import itertools
import math
import pickle
import glob
import os
import mdtraj as md
import time

def ComputeContacts(trajfile, res):
    '''
    Compute residue-residue contacts for an .xtc trajectory file
    Args:
        trajfile (str): .xtc file stripped with CPPTRAJ to contain only protein
        res (np.array): shape (number of res, 2), contains the residue pairs to compute the distance between
    Returns:
        a list of np.array as a saved .pkl file with the distances;
        each array corresponds to a trajectory and is of shape (num of frames, num of residue pairs)
    '''
    name = trajfile.replace('_stripped_wat.xtc','')

    #final file name
    dist_file = name + '-FINAL-FEATURES.pkl'
    a = os.path.isfile(dist_file)
    if a == False:
        #search for parmtop file
        topname_parm = glob.glob('*stripped_wat.p*m*')[0]

        traj = md.load(trajfile, top = topname_parm)
        top = traj.topology
        
        contact_array = md.compute_contacts(traj = traj,
                                            contacts = res,
                                            scheme = 'ca')[0]

        outfile = open(name + '-FINAL-FEATURES.pkl','wb')
        pickle.dump(contact_array, outfile)
        outfile.close()

if __name__ == '__main__':

    #load np.array file that contains the residues
    #could also input residues directly, 
    #as long as res is an np.array of shape (number of pair, 2)
    res = np.load("FINAL-RES-PAIRS-PLPG.npy")
    
    #loop through xtcs containing only proteins 
    for file in sorted(glob.glob('*MD*stripped_wat.xtc')):
        ComputeContacts(file, res)
