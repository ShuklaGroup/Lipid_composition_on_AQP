import pickle
import numpy as np

def GetHOLE(filename):
	'''
	Processes dictionary outputs from HOLE MDAnalysis
	Args:
		filename (str): name of pickle file containing the HOLE output of a trajectory
	Returns:
		(np.array): shape (number of frames, number of bins) binned radii for each frame in the directory, all of equal length
		saved as a pickle file
	'''
	#load file
	holefile = open(filename,"rb")
	holedict = pickle.load(holefile)
	holefile.close()
  
	#define edges of bins corresponding to the z-coord
	#(consistent across all trajectories)
	#HOLE z-coordinate outputs are 0.1 apart, so we expect to find one radius value per bin
	edges = np.arange(-20,20,0.1)
	tot_radius = []

	#iterate through each frame
	for frame in range(len(holedict['rxn_coord'])):
    
		#initialize binned radius array for current frame
		radius = np.empty(len(edges)-1)
		radius[:] = np.nan
		frame_z = holedict['rxn_coord'][frame]
		frame_r = holedict['radius'][frame]
    
		#per z-position bin, find the corresponding radius
		for i in range(len(edges)-1):    
			for j in range(len(frame_z)):
				if edges[i] <= frame_z[j] <= edges[i+1]:
					radius[i] = frame_r[j]
		tot_radius.append(radius)
	tot_radius_arr = np.array(tot_radius)
  
	#save binned radius file
	radiusfilename = filename.replace("-HOLE-DICT.pkl","-HOLE-RADIUS.pkl")
	radiusfile = open(radiusfilename, 'wb')
	pickle.dump(tot_radius_arr, radiusfile)
	radiusfile.close()

if __name__ == '__main__':
	import argparse

	ps = argparse.ArgumentParser()
	ps.add_argument('filename', type = str, help = '-HOLE-dict.pkl file input')
	args = ps.parse_args()
	filename = args.filename
	GetHOLE(filename)
