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
from pyemma.coordinates import tica
from datetime import datetime
date = datetime.today().strftime('%Y%m%d')

#Create numpy file including the mean-first passage time between the open-like and closed-like cluster 

#Read in contact files
file_list = sorted(glob.glob('*-features-2023*.pkl'))
system = file_list[0].replace('_','-').split("-")[0]
tica_obj = pickle.load(open(glob.glob(system + '*dim*cluster*tica_obj.pkl')[0], 'rb'))
tica_trajs = pickle.load(open(glob.glob(system + '*dim*cluster*tica_trajs.pkl')[0], 'rb'))
tica_concatenate = np.concatenate(tica_trajs)
print('length: ',len(tica_concatenate)/10**5)

#Read in cluster files and MSM objects
cluster_obj = pickle.load(open(glob.glob(system + '*dim*cluster*cluster_obj.pkl')[0], 'rb'))
cluster_trajs = cluster_obj.dtrajs
weights = pickle.load(open(system+'-weights.pkl','rb'))
msm = pickle.load(open(system+'-msm_obj.pkl','rb'))

data = []
for file in sorted(glob.glob('*features-2023*.pkl')):
    filedata = pickle.load(open(file,'rb'))
    data.append(filedata)

opencrys_feat = np.load(glob.glob("*2b5f*CC-features-parmtop.npy")[0])
closedcrys_feat = np.load(glob.glob("*1z98*CC-features-parmtop.npy")[0])

open_min = 100
open_traj = 0
open_frame = 0
for traj in range(len(data)):
    for frame in range(len(data[traj])):
        diff = np.linalg.norm(data[traj][frame] - opencrys_feat, 2)
        if diff < open_min:
            open_min = diff
            open_traj = traj
            open_frame = frame
print("open",open_traj,open_frame,open_min)

closed_min = 100
closed_traj = 0
closed_frame = 0
for traj in range(len(data)):
    for frame in range(len(data[traj])):
        diff = np.linalg.norm(data[traj][frame] - closedcrys_feat, 2)
        if diff < closed_min:
            closed_min = diff
            closed_traj = traj
            closed_frame = frame
print("closed",closed_traj,closed_frame,closed_min)

opencc = tica_trajs[open_traj][open_frame]
closedcc = tica_trajs[closed_traj][closed_frame]

#find mfpt between cc
import numpy.linalg as la
opendiff = cluster_obj.clustercenters - opencc
opennorm = []
for i in range(len(opendiff)):
    opennorm.append(la.norm(opendiff[i], 2))
opencluster = np.where(opennorm == min(opennorm))[0][0]

closeddiff = cluster_obj.clustercenters - closedcc
closednorm = []
for i in range(len(closeddiff)):
    closednorm.append(la.norm(closeddiff[i], 2))
closedcluster = np.where(closednorm == min(closednorm))[0][0]

mfptarray = [
    msm.mfpt(closedcluster,opencluster)/(10**4),
    msm.mfpt(opencluster,closedcluster)/(10**4)
    ]

np.save("%s-mfpt-%s.npy"%(system,date),mfptarray)