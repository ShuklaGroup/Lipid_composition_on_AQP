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

#Read in contact files
file_list = sorted(glob.glob('*-features-2023*.pkl'))
system = file_list[0].replace('_','-').split("-")[0]

#Read in msm objects
cluster_obj = pickle.load(open(glob.glob(system + '*dim*cluster*cluster_obj.pkl')[0], 'rb'))
cluster_trajs = cluster_obj.dtrajs
weights = pickle.load(open(system+'-weights.pkl','rb'))
msm = pickle.load(open(system+'-msm_obj.pkl','rb'))

weights_array = weights
raw_cluster_array = cluster_trajs

raw_array = np.concatenate(cluster_trajs)
weighted_array = np.concatenate(weights_array)
num_clusters = max(raw_array) + 1

#Initialize
raw_counts = np.zeros(num_clusters, dtype=float)
msm_counts = np.zeros(num_clusters, dtype=float)

#Iterate through each cluster, count number of states in the cluster for the raw and MSM-reweighted data
for num in range(num_clusters):
    ind = np.where(raw_array == num)[0]
    raw_counts[num] += len(ind)
    msm_counts[num] += np.sum(weighted_array[ind])

raw_prob = raw_counts/len(raw_array)
log_raw_prob = np.log10(raw_prob)
log_msm_prob = np.log10(msm_counts)

np.save(system+'_raw_probabilites.npy',log_raw_prob)
np.save(system+'_msm_probabilites.npy',log_msm_prob)

#get system color/if dont have pkl file, set to any color
color = pickle.load(open('system-colors.pkl','rb'))[system][-2]

#plotting
fig, axs = plt.subplots(1,1,figsize=(6,5))
xmin = min(log_raw_prob)
xmax = max(log_raw_prob)
x = np.linspace(-5,0,100)
plt.plot(x,x,'k-')
axs.scatter(log_raw_prob, log_msm_prob , color=color)
plt.xlabel('Raw Counts Probability ($log_{10}$)', fontsize=18)
plt.ylabel('MSM Counts Probability ($log_{10}$)', fontsize=18)
for axis in ['top','bottom','left','right']:
    axs.spines[axis].set_linewidth(2)
axs.tick_params(width=2)
plt.xticks(np.arange(-5,1,1),fontsize=12,fontweight='bold')
plt.yticks(np.arange(-5,1,1),fontsize=12,fontweight='bold')
plt.xlim([-5,0])
plt.ylim([-5,0])
plt.title(system, fontsize=18,fontweight='bold')
plt.savefig("%s-reweighting-plot-%s.png"%(system,date),dpi=300,bbox_inches='tight')
plt.clf()