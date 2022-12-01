import mdtraj as md
import numpy as np
import pickle
import matplotlib
import matplotlib.pyplot as plt
matplotlib.use("Agg")
import glob
import os
import pyemma
import time
import pandas as pd
from pyemma.coordinates import tica
from multiprocessing import Pool
import itertools

all_contact_files = sorted(glob.glob('*MD*-CONTACT.pkl'))

#by naming convention, extract the bilayer system name
system = all_contact_files[0].replace('_','-').split('-')[0]

#Set up directories
os.makedirs('./RRCS-out-abs-no-Nterm/',exist_ok=True)
outdir = "./RRCS-out-abs-no-Nterm/"

#Read in contact files
for file in all_contact_files:

    out_file = file.replace('-CONTACT.pkl','')

    localtime = time.asctime(time.localtime(time.time()))
    print("Start loading file:", file, localtime)

    feat_file = pickle.load(open(file, 'rb'))

    #Take only the pair-wise distance arrays, excluding res 0-4 (Python index)
    feat = feat_file[0][:,1155:] 

    #Total number of contact features
    num_features = len(feat[0])
    num_frames = len(feat)
    print("Number of features is ", num_features)

    #Load crystal structure data
    open_contact = pickle.load(open('OPEN-STRUCTURE-CONTACT.pkl', 'rb'))[0]
    closed_contact = pickle.load(open('CLOSED-STRUCTURE-CONTACT.pkl', 'rb'))[0]
    
    #excludes res 0-4 (Python index)
    open_contact = open_contact[0][1155:]
    closed_contact = closed_contact[0][1155:]
    print(len(open_contact))
    assert(len(open_contact) == len(closed_contact))
    assert(len(open_contact) == num_features)
    
    print("Feature matrix size:", feat.shape)
    localtime = time.asctime(time.localtime(time.time()))
    print("Start matrix computations:", localtime)

    #compute difference with open and closed structure (delta d open, delta d closed)
    feat_ddo = feat - open_contact.T
    feat_ddc = feat - closed_contact.T

    #compute avg movement in each frame
    feat_ddo_mean = np.mean(feat_ddo, axis=1)
    feat_ddc_mean = np.mean(feat_ddc, axis=1)

    feat_ddo_mean = feat_ddo_mean.reshape(num_frames,1)
    feat_ddc_mean = feat_ddc_mean.reshape(num_frames,1)

    #normalize delta d, across all res-res contacts in one frame
    feat_rrcs_o = (feat_ddo - feat_ddo_mean) / feat_ddo_mean
    feat_rrcs_c = (feat_ddc - feat_ddc_mean) / feat_ddc_mean

    assert(feat_rrcs_o.shape == feat.shape)
    assert(feat_rrcs_c.shape == feat.shape)

    localtime = time.asctime(time.localtime(time.time()))
    print("Finished matrix computations:", localtime)

    #compute sum of normalized distances across all frames:
    traj_rrcs_o_sum = np.sum(np.abs(feat_rrcs_o), axis = 0)
    traj_rrcs_c_sum = np.sum(np.abs(feat_rrcs_c), axis = 0)


    #save files
    open_sum_file = open(outdir + out_file + '-OPEN-ABS-SUM-RRCS.pkl', 'wb')
    pickle.dump(traj_rrcs_o_sum,open_sum_file)
    open_sum_file.close()

    closed_sum_file = open(outdir + out_file + '-CLOSED-ABS-SUM-RRCS.pkl', 'wb')
    pickle.dump(traj_rrcs_c_sum,closed_sum_file)
    closed_sum_file.close()

    #record time
    localtime = time.asctime(time.localtime(time.time()))
    print("Finished file:", file, localtime)
