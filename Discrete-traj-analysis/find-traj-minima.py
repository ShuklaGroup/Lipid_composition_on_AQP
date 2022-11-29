import numpy as np
import pandas as pd
import pickle
import mdtraj as md
import itertools
import glob
import os

def GetFrameIndex(tica_trajs, lims):
	"""
	Finding the frames in the selected tICA region
	Args: 
		tica_trajs (list): list of arrays containing outputs from tICA, 
			each array corresponds to a trajectory with shape (number of frames, number of tIC features)
		lims (np.array): shape (4,) array contraining maxs and mins of tIC1 and tIC2
	Returns: 
		frame_ind_all (list): list of arrays containing frame indices that reside in tic1 and tic2 range
			each array corresponds to a trajectory with shape (number of frames, number of tIC features)
			(there could be empty arrays for trajectories that has no frames in the selected region)
	"""
	xmin, xmax, ymin, ymax = lims
	frame_ind_all = []
	for i in range(len(tica_trajs)):
		traj = tica_trajs[i]

		#get boolean array that satisfy input limits for tic1 and tic2
		xbool = np.logical_and(traj[:,0] > xmin, traj[:,0] < xmax)
		ybool = np.logical_and(traj[:,1] > ymin, traj[:,1] < ymax)
		final_bool = np.logical_and(xbool, ybool)

		#extract frame indeces
		frame_ind = np.where(final_bool)[0] + 1
		frame_ind_all.append(frame_ind)
	return frame_ind_all



def GenerateCpptrajInput(state):
	'''
	For current state, save .in files as cpptraj inputs to obtain .xtc files corresponding to found frames
	Args:
		state (int): current state, 1-indexed
	Returns:
		.in files as cpptraj input to generate .xtc files of the frames from GetFrameIndex()
	'''
	#load csv of the x and y limits
	box_df = pd.read_csv(system + '-state%i.csv'%state, names=['x', 'y']) 
	array = np.array([min(box_df['x']), max(box_df['x']), min(box_df['y']), max(box_df['y'])])
	frame_ind_all = GetFrameIndex(tica_trajs, array)
	assert(len(trajnames) == len(frame_ind_all))

	#save selected frames indeces for testing
	pickle.dump(frame_ind_all, open(system + '-frames-state%i.pkl'%state, 'wb'))

	#for each trajectory, write a txt file containing cpptraj inputs
	for j in range(len(frame_ind_all)):
		frames = frame_ind_all[j]
		if len(frames) != 0: 
			traj_name = trajnames[j]
			f = open(traj_name + '-state-%i.txt'%state, 'w')
			f.write('parm *[0-9].p*m*' + '\n')
			for k in range(len(frames)):
				f.write('trajin ' + traj_name + '.xtc' + ' ' + str(frames[k]) + ' ' + str(frames[k]) + ' ' + str(1))
				f.write('\n')
			f.write('trajout ' + traj_name + '-state-%i.xtc'%state + '\n')
			f.write('run' + '\n')
			f.write('quit' + '\n')
			f.close()

if __name__ == '__main__':
	#find the tica output in the directory
	tica_file = glob.glob('*-tica_trajs.pkl')[0]
	tica_trajs = pickle.load(open(tica_file, 'rb'))
	system = tica_file.split('-')[0]

	#get names of trajectories through feature files used for tICA
	trajnames = []
	for file in sorted(glob.glob('*-FINAL-FEATURES.pkl')):
		name = file.replace('-FINAL-FEATURES.pkl','')
		trajnames.append(name)
	
	#iterate through each state
	for i in range(1,5):
		GenerateCpptrajInput(state=i)
