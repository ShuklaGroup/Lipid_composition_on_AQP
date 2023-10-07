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

#Read in files
file_list = sorted(glob.glob('*-features-2023*.pkl'))
system = file_list[0].replace('_','-').split("-")[0]
tica_obj = pickle.load(open(glob.glob(system + '*dim*cluster*tica_obj.pkl')[0], 'rb'))
tica_trajs = pickle.load(open(glob.glob(system + '*dim*cluster*tica_trajs.pkl')[0], 'rb'))
tica_concatenate = np.concatenate(tica_trajs)
print('length: ',len(tica_concatenate)/10**5)


# cluster_trajs = cluster_obj.dtrajs
weights = pickle.load(open(system+'-weights.pkl','rb'))
# msm = pickle.load(open(system+'-msm_obj.pkl','rb'))

#read in all feature data and crystal structure data
data = []
for file in sorted(glob.glob('*features-2023*.pkl')):
    filedata = pickle.load(open(file,'rb'))
    data.append(filedata)

opencrys_feat = np.load(glob.glob("*2b5f*CC-features-parmtop.npy")[0])
closedcrys_feat = np.load(glob.glob("*1z98*CC-features-parmtop.npy")[0])

#loop through all frames to find the one closest to the cystal structures
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

#save frame corresponding to crystal structures projected on tics
opencc = tica_trajs[open_traj][open_frame] #np.load(system+"-opencc.npy")
closedcc = tica_trajs[closed_traj][closed_frame] #np.load(system+"-closedcc.npy")

#find the most correlated features to the two principal tics
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

tic1name = pairlist[np.where(corr[:,0] == max(corr[:,0]))[0][0]] + ' (' + str(np.round(max(corr[:,0]),2)) + ')'
tic2name = pairlist[np.where(corr[:,1] == max(corr[:,1]))[0][0]] + ' (' + str(np.round(max(corr[:,1]),2)) + ')'

#plot with all information
#%matplotlib inline
from matplotlib import pyplot as plt
from matplotlib.colors import ListedColormap
import tol_colors as tc
from datetime import datetime

date = datetime.today().strftime('%Y%m%d')
purple=['#6247aa','#815ac0','#a06cd5','#b185db','#d2b7e5']
blue=['#2c7da0','#468faf','#61a5c2','#89c2d9','#a9d6e5']
green=['#718355', '#87986a', '#97a97c', '#a3b18a', '#cfe1b9']
orange=["#ffb700","#ffc300","#ffd000","#ffdd00","#ffea00"]
red=['#f25c54', '#f27059', '#f4845f', '#f79d65', '#f7b267']
#cornsilk=['#F1ECCE']
larger=['#f7b267','#f7b267','#f7b267','#f7b267','#f7b267']
colors = purple+blue+green+orange+red+larger

#define parameters
t = 27.6
R = 0.001987
T = 300
fig_wid = 9
fig_hig = 6
#cmap = mpl.cm.jet
Max_energy = 6

#get data
X = tica_concatenate[:,0].T
Y = tica_concatenate[:,1].T

#define limits of bins and plot
x_min = np.min(X)
y_min = np.min(Y)
x_max = np.max(X)
y_max = np.max(Y)

x_hist_lim_low = x_min - 0.5
y_hist_lim_low = y_min - 0.5
x_hist_lim_high = x_max + 0.5
y_hist_lim_high = y_max + 0.5
x_lim_low = (int(np.min(X)/5))*5.0
y_lim_low = (int(np.min(Y)/5))*5.0
x_lim_high = (int(np.max(X)/5) + 1)*5.0
y_lim_high = (int(np.max(Y)/5) + 1)*5.0

#construct 2d histogram from 2 principle tics
hist= np.histogram2d(X,Y,bins=175,
           range = [[x_hist_lim_low,x_hist_lim_high],[y_hist_lim_low,y_hist_lim_high]],
           density = True, weights = np.concatenate(weights))

#define bins from 2d histogram
prob_density = hist[0]
xedge = hist[1]
yedge = hist[2]
x_bin_size = xedge[1]-xedge[0]
y_bin_size = yedge[1]-yedge[0]

#calculate relative free energy from density distribution
free_energy = -R*T*np.log(prob_density*x_bin_size*y_bin_size)
min_free_energy= np.min(free_energy)
delta_free_energy = free_energy - min_free_energy
xx = [(xedge[i]+xedge[i+1])/2 for i in range(len(xedge)-1)]
yy = [(yedge[i]+yedge[i+1])/2 for i in range(len(yedge)-1)]
fig, axs = plt.subplots(1,1,figsize=(fig_wid,fig_hig))
#plt.contour(xx,yy,delta_free_energy.T, linewidths=0.4, levels=[0,1,2,3,4,5], colors='black')

#construct contour plot for free energy
cd = axs.contourf(xx,yy,delta_free_energy.T, np.linspace(0,Max_energy,Max_energy*5+1),
         vmin=0.0, vmax=Max_energy,colors=colors)
cbar = fig.colorbar(cd,ticks=range(Max_energy+1))
cbar.ax.set_yticklabels(range(Max_energy+1), fontsize=22, fontname='Arial', fontweight='bold')
cbar.ax.set_ylabel('Free energy (kcal/mol)', fontsize=22, fontname='Arial', fontweight='bold', rotation=270, labelpad=20)
cbar.outline.set_linewidth(2)
cbar.ax.tick_params(width=2)

#plot two crystal structures frames projected on tic space
plt.plot(opencc[0], opencc[1], color='#EF233C', markeredgecolor='black', marker='o', markersize=10)
plt.plot(closedcc[0], closedcc[1], color='#54F2F2', markeredgecolor='black', marker='o', markersize=10)

#aesthetics
plt.title('SoPIP2:' + system, fontsize=18, fontname='Arial',fontweight='bold') # + ' tICA Decomposition'
plt.xlabel('tIC1:%s'%tic1name, fontsize=22, fontname='Arial', fontweight='bold')
plt.ylabel('tIC2:%s'%tic2name, fontsize=22, fontname='Arial', fontweight='bold')

# change all spines
for axis in ['top','bottom','left','right']:
    axs.spines[axis].set_linewidth(2)
axs.tick_params(width=2)
plt.xticks(fontsize=22, fontweight='bold')
plt.yticks(fontsize=22, fontweight='bold')
plt.tight_layout()

#save figures
plt.savefig(system+'_tic_cor_feat'+date+'_600dpi.png',dpi=600)
plt.savefig(system+'_tic_cor_feat'+date+'_300dpi.png',dpi=300)
plt.clf()
