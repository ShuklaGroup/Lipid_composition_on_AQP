import pathlib
import pickle  # this library is used for saving the analysis objects to file

import numpy as np
import MDAnalysis as mda
from lipyphilic.lib.assign_leaflets import AssignLeaflets, AssignCurvedLeaflets
from lipyphilic.lib.order_parameter import SCC

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

import argparse
import glob
import pickle
from datetime import datetime

#argument parser to run this analysis for a given *.xtc file
ps = argparse.ArgumentParser()
ps.add_argument('filename', type = str, help = '-stripped-wat.xtc frame file input with protein and water')
args = ps.parse_args()
filename = args.filename
system = filename.replace('_','-').split("-")[0]
name = filename.replace(".xtc","")

#get protein state
if "2b5f" in filename:
    protein = "2b5f"
else:
	protein = "1z98"

#initialize MDAnalysis universe with *.gro coordinate file and xtc trajectory
u = mda.Universe(glob.glob("%s_%s_*wat.gro"%(system, protein))[0],
                filename)

#assign leaflet with MDAnalysis detection
leaflets = AssignCurvedLeaflets(
    universe=u,
    lipid_sel='name P O13 O14'
    #document states 3 atoms https://lipyphilic.readthedocs.io/en/latest/reference/lib/areas.html
)
#run leaflet detection
leaflets.run(verbose=True)

#caculate order parameter for each tail
scc_sn1 = SCC(
    universe = u, 
    tail_sel = "name C23 C24 C25 C26 C27 C28 C29 C210 C211 C212 C213 C214 C215 C216 C217 C218"
)
scc_sn1.run(start=None, stop=None, step=None, verbose=True)
scc_sn2 = SCC(
    universe = u, 
    tail_sel = "name C33 C34 C35 C36 C37 C38 C39 C310 C311 C312 C313 C314 C315 C316"
)
scc_sn2.run(start=None, stop=None, step=None, verbose=True)

#average both tails order params
scc = SCC.weighted_average(scc_sn1, scc_sn2)

#average across all trajectories
avg_scc = np.mean(scc.SCC,axis=1)
std_scc = np.std(scc.SCC,axis=1)

pickle.dump(scc_sn1, open(name + "-scc_sn1.pkl","wb"))
pickle.dump(scc_sn2, open(name + "-scc_sn2.pkl","wb"))
pickle.dump(avg_scc, open(name + "-avg_scc.pkl","wb"))
pickle.dump(std_scc, open(name + "-std_scc.pkl","wb"))

#PLOTTING
#protein helices positions for plotting
helices = [
    "resid 10 to 35 and name CA", #TM1
    "resid 46 to 66 and name CA", #TM2
	  "resid 75 to 84 and name CA", #top half
    "resid 90 to 114 and name CA", #TM3
    "resid 135 to 155 and name CA", #TM4
	  "resid 172 to 189 and name CA", #TM5
    "resid 196 to 207 and name CA", #bottom half
    "resid 208 to 234 and name CA", #TM6
    ]
protein_x = []
protein_y = []
for i, helix in enumerate(helices):
	x = []
	y = []
	#trajfile = filename
	universe = u
	helix_sel = universe.select_atoms(helix)
	for ts in universe.trajectory[0:10000:10]:
		x.append(helix_sel.positions[:,0])
		y.append(helix_sel.positions[:,1])
	protein_x.append(np.concatenate(x))
	protein_y.append(np.concatenate(y))
pickle.dump(protein_x,open(name + "-protein-x.pkl", "wb"))
pickle.dump(protein_y,open(name + "-protein-y.pkl", "wb"))

#lipid positions
#all_lipid_ids = u.select_atoms("name P or name O3").resids
scc_lipids_x = u.select_atoms("name P or name O3").positions[:,0]
scc_lipids_y = u.select_atoms("name P or name O3").positions[:,1]

pickle.dump(scc_lipids_x,open(name + "-lipid-x.pkl", "wb"))
pickle.dump(scc_lipids_y,open(name + "-lipid-y.pkl", "wb"))

plt.figure(figsize=(5,6))

#plot protein based on helix color
colors = ["#390099","#73d2de","#9e0059","#ffc2e2","#fc2f00","#ec7d10","#04a777","#ffbd00"]
custom_lines = []
for i in range(len(protein_x)):
	plt.hexbin(protein_x[i],protein_y[i],
	    gridsize=100,bins='log',color=colors[i],mincnt=1,label=str(i),
		alpha=0.1)
	custom_lines.append(Line2D([0],[0],color=colors[i],lw=4))

#scatter each lipid xy position and color based on its average order parameter
plt.scatter(scc_lipids_x,
            scc_lipids_y,
            c=avg_scc,
            #alpha=0.5,
			      s=50,
			      edgecolors='k',
            cmap='coolwarm')

cbar = plt.colorbar(orientation='horizontal',
                    pad=0.1)
#plt.clim(0.05,0.35)
#cbar.set_ticks(np.arange(21,23,0.1))
cbar.set_label("$S_{cd}$")
plt.xlim([-40,40])
plt.ylim([-40,40])
plt.xticks(np.arange(-40,41,10))
plt.yticks(np.arange(-40,41,10))
# custom_lines = [Line2D([0], [0], color='#5E4C5A', lw=4),
#                 Line2D([0], [0], color='#E1F0C4', lw=4),
# 				]
#plt.legend(custom_lines,['TM1','TM2','Half1','TM3','TM4','TM5','Half2','TM6'],
#           fontsize=12, loc=(1.04, 0))
plt.tight_layout()
date = datetime.today().strftime('%Y%m%d')
plt.savefig(name+"-order-per-lip-%s-300dpi.png"%date,dpi=300)
plt.savefig(name+"-order-per-lip-%s-600dpi.png"%date,dpi=600)
