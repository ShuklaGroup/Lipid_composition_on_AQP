#water counts adapted from Gelenter, M.D., et al. Commun Biol 4, 338 (2021). https://doi.org/10.1038/s42003-021-01847-2
import MDAnalysis as mda
import numpy as np
import numpy.linalg as npla
import argparse
import glob
import pickle
import time as timepackage
from datetime import datetime

#take filename argument to process each *.xtc file
ps = argparse.ArgumentParser()
ps.add_argument('filename', type = str, help = '-stripped-state-*.xtc frame file input with protein and water')
args = ps.parse_args()
filename = args.filename
system = filename.replace('_','-').split("-")[0]

#search for macrostate of the trajectory
for state_name in ["c0","c1","c2","i1","i2","i3","i4","o0"]:
     if state_name in filename:
        state = state_name

trajname = filename.replace("_stripped-state-%s.xtc"%state,'')

#define some pore lining residues
pore = "resid 47 or resid 206 or resid 211 or "\
"resid 29 or resid 107 or resid 218 or resid 50 or resid 138 or resid 185 or "\
"resid 21 or resid 58 or resid 222 or "\
"resid 13 or resid 92 or resid 65 or resid 81 or resid 152 or resid 230 "\
"and name CA"

#define normalized channel axis
def channelaxis(dirv,cv):
    '''
    Args: 
    	dirv, numpy array defining direction of channel (intra - extracellular central point)
    	cv, principal z channel axis
    Returns: normalized channel axis
    '''
    r = np.dot(dirv,cv)
    if r<0:
        return -1*cv/npla.norm(cv)
    else:
        return cv/npla.norm(cv)

# minz = -20.0 # Minimum channel coordinate
# maxz = 20.0 # Maximum channel coordinate

# Search for coordinate files
if "2b5f" in filename:
	pdb = glob.glob("*2b5f*stripped.pdb")[0]
if "1z98" in filename:
	pdb = glob.glob("*1z98*stripped.pdb")[0]

# Initialize MDAnalysis traj objects and arrays
trajectory = filename
universe = mda.Universe(pdb, trajectory)
transported = np.zeros((len(universe.trajectory),3)) # Time, Transport to N-term, Transport to C-term
passageTimesO = [] # Empty list to contain the passage times in the N-terminal direction
passageTimesI = [] # Empty list to contain the passage times in the C-terminal direction
transportedresO = []
transportedresI = []

# Select relevant residues

extracellular_region = universe.select_atoms("protein and (resid 29 or resid 107 or resid 218 or resid 50 or resid 138 or resid 185) and name CA")
intracellular_region = universe.select_atoms("protein and (resid 13 or resid 92 or resid 65 or resid 81 or resid 152 or resid 230) and name CA")
pore_region = universe.select_atoms(pore) # All Calphas
waters = universe.select_atoms("resname WAT and name O") # Water oxygens

# Initialize location and time inside array
water_location = np.zeros((len(waters),1),dtype=int)
water_time = np.zeros((len(waters),1),dtype=float)

# Collect data
localtime = timepackage.asctime(timepackage.localtime(timepackage.time()))
print("Start data collection", localtime)

# Iterate through each timestep
for ts in universe.trajectory:
    time = int(ts.time)
    ti = int(time * 0.1)
    inward_direction = intracellular_region.centroid() - extracellular_region.centroid()
    x, y, z = pore_region.principal_axes()
    caxis = channelaxis(inward_direction, z)
    mincor = intracellular_region.centroid() #z of the intracellular centroid
    maxcor = extracellular_region.centroid() #z of the extracellular centroid
    minz = np.dot(mincor - maxcor, caxis)
    
    # Define max/min of x and y coordinates by going through all residues
    x_lims = [np.min(pore_region.positions[:,0]), np.max(pore_region.positions[:,0])]
    y_lims = [np.min(pore_region.positions[:,1]), np.max(pore_region.positions[:,1])]

    # Iterate through each water and check its previous position:
	# 1: water outside channel, in the extracellular side
	# 2: water outside channel, in the intracellular side
	# 3: water inside channel, entered from extracellular side
	# 4: water inside channel, entered from intracellular side
	# then assign its current position
    for wi in range(len(waters)):
        
		#find projection of water onto the channel axis
        water_pos = waters[wi].position - extracellular_region.centroid()
        water_pos_on_caxis = np.dot(water_pos, caxis)
        
		#inside chanel is defined as 0 < water_pos_on_caxis < minz and xy bounded by protein

        if (water_pos_on_caxis<0): # Water on extracellular side
            
			#check its previous position
			#if this water was from the intracellular side and then went inside the channel, count as EXPORT
            if water_location[wi]==4:
                transported[ti-1][1] += 1
                passageTimesO.append(time-water_time[wi]) # Tranport complete, save how many frames it took
                transportedresO.append(waters[wi].resid)
            
			#if previous position not inside channel, assign new position
            water_location[wi] = 1
            
        elif (water_pos_on_caxis>minz): # Water on intracellular side
            
			#check its previous position
			#if this water was from the extracellular side and then went inside the channel, count as IMPORT
            if water_location[wi]==3:
                transported[ti-1][2] += 1 
                passageTimesI.append(time-water_time[wi]) # Tranport complete, save how many frames it took
                transportedresI.append(waters[wi].resid)
                
			#if previous position not inside channel, assign new position
            water_location[wi] = 2
            
        elif (
            0 < water_pos_on_caxis < minz) and (
            x_lims[0] < waters[wi].position[0] < x_lims[1]) and (
            y_lims[0] < waters[wi].position[1] < y_lims[1]): # Water is in the channel
            
			#check previous position

			#if water was from the intracellular side
            if water_location[wi]==1:
                
				#assign its current position and that it comes from the intracellular side
                water_location[wi] = 3
                water_time[wi] = time #save the frame at which this water entered channel
            
			#if water was from the extracellular side
            elif water_location[wi]==2:
                
                #assign its current position and that it comes from the extracellular side
                water_location[wi] = 4
                water_time[wi] = time #save the frame at which this water entered channel
                
localtime = timepackage.asctime(timepackage.localtime(timepackage.time()))
print("Finish data collection", localtime)

date = datetime.today().strftime('%Y%m%d')

#save files for export, import, and total transport
pickle.dump(transported, open("%s-transport-state-%s-%s.pkl"%(trajname,state,date),"wb"))
pickle.dump(passageTimesO, open("%s-passageTimesO-state-%s-%s.pkl"%(trajname,state,date),"wb"))
pickle.dump(passageTimesI, open("%s-passageTimesI-state-%s-%s.pkl"%(trajname,state,date),"wb"))
pickle.dump(transportedresI, open("%s-transportedresI-state-%s-%s.pkl"%(trajname,state,date),"wb"))
pickle.dump(transportedresO, open("%s-transportedresO-state-%s-%s.pkl"%(trajname,state,date),"wb"))
