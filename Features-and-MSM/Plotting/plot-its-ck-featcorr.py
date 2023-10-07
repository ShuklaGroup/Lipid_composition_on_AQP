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
tica_obj = pickle.load(open(glob.glob(system + '*dim*cluster*tica_obj.pkl')[0], 'rb'))
tica_trajs = pickle.load(open(glob.glob(system + '*dim*cluster*tica_trajs.pkl')[0], 'rb'))
tica_concatenate = np.concatenate(tica_trajs)
print('length: ',len(tica_concatenate)/10**5)

cluster_obj = pickle.load(open(glob.glob(system + '*dim*cluster*cluster_obj.pkl')[0], 'rb'))
its = pickle.load(open(glob.glob(system + '*dim*cluster*its_obj.pkl')[0], 'rb'))
pyemma.plots.plot_implied_timescales(its,units='ns',dt=0.01)
plt.title(system)
plt.savefig("%s-ITS-%s-dpi300.png"%(system, date),dpi=300,bbox_inches='tight')
plt.clf()

cluster_trajs = cluster_obj.dtrajs
weights = pickle.load(open(system+'-weights.pkl','rb'))
msm = pickle.load(open(system+'-msm_obj.pkl','rb'))

#ck
pyemma.plots.plot_cktest(msm.cktest(3),units='ns',dt=0.01)
plt.title(system)
plt.savefig("%s-CK-%s-dpi300.png"%(system, date),dpi=300,bbox_inches='tight')
plt.clf()

#features and correlation with tics
corr = tica_obj.feature_TIC_correlation
feat = pickle.load(open(system+"-features.pkl","rb"))
pdbfile = glob.glob(system + "*wat_lip.pdb")[0]
traj = md.load(pdbfile, top=pdbfile)
topology = traj.topology
pairlist = []
for pair in feat:
    res1str = str(topology.atom(topology.select("resid %s"%pair[0])[0]).residue)[:3]
    res1 = res1str[0] + res1str.lower()[1:3] + str(pair[0]+28)
    
    res2str = str(topology.atom(topology.select("resid %s"%pair[1])[0]).residue)[:3]
    res2 = res2str[0] + res2str.lower()[1:3] + str(pair[1]+28)
    string = res1 + '-' + res2
    pairlist.append(string)
    
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.colors
import seaborn as sns
# colors=[matplotlib.colors.to_rgb('#B4FBB8')]
# cm = LinearSegmentedColormap.from_list(name='a',colors=colors,N=100)
fig, ax = plt.subplots(figsize=(3, 8))
sns.heatmap(tica_obj.feature_TIC_correlation[:,:], vmin=-1, vmax=1, cmap="coolwarm", square=True, linewidth=0.5)
ax.set_yticklabels(pairlist, rotation=0)
ax.set_xticklabels(np.arange(len(tica_obj.feature_TIC_correlation))+1)
ax.set_xlabel("tICA Components")
ax.set_ylabel("Distance Features")
plt.title(system)
plt.savefig("%s-tic-features-%s.png"%(system,date),dpi=300,bbox_inches='tight')
plt.clf()