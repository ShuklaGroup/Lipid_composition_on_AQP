#extract frames from tica
import numpy as np
import pandas as pd
import pickle
import mdtraj as md
import glob
import os

def get_frame_index(tica_trajs, lims, sample_size=1000):
	"""
	args: 
		tica_trajs: list of arrays containing tica outputs 
		lims: (4,) ndarray contraining maxs and mins of tic1 and tic2
	return: list of arrays containing frame indices that satisfy specified tic1 and tic2 range
	"""
	xmin, xmax, ymin, ymax = lims
	traj_and_frame = []
	num_traj = len(tica_trajs)

	for i in range(num_traj):
		traj = tica_trajs[i]

		#get boolean array that satisfy input limits for tic1 and tic2
		xbool = np.logical_and(traj[:,0] > xmin, traj[:,0] < xmax)
		ybool = np.logical_and(traj[:,1] > ymin, traj[:,1] < ymax)
		final_bool = np.logical_and(xbool, ybool)

		#extract frame indeces
		frame_ind = np.where(final_bool)[0] + 1
		#frame_ind_all.append(frame_ind)

		#save traj and frame index for random selection
		if len(frame_ind) != 0:
			for frame in frame_ind:                                                                                                                                                                                      
				traj_and_frame.append(np.array([i, frame]))
	
	traj_and_frame = np.array(traj_and_frame)

	#randomly select for sample_size number of frames
	sample_ind = np.random.choice(len(traj_and_frame), sample_size, replace=False)
	sample = traj_and_frame[sorted(sample_ind)]

	#save selected sample in [n_trajs,n_frame] structure
	frame_ind_all = []
	for i in range(num_traj):
		frame_ind = []
		for traj_frame_pair in sample:
			if traj_frame_pair[0] == i:
				frame_ind.append(traj_frame_pair[1])
		frame_ind_all.append(frame_ind)

	return frame_ind_all

#take argument for this?
tica_file = glob.glob('*-tica_trajs.pkl')[0]
tica_trajs = pickle.load(open(tica_file, 'rb'))
system = tica_file.split('-')[0]

#get names of trajectories
trajnames = []
for file in sorted(glob.glob('*-features-2023*.pkl')):
    name = file[:-22]
    trajnames.append(name)

for lim_file in sorted(glob.glob("*-state-*.csv")): #each state
	box_df = pd.read_csv(lim_file, names=['x', 'y']) #load csv of the x and y limits
	array = np.array([min(box_df['x']), max(box_df['x']), min(box_df['y']), max(box_df['y'])])
	frame_ind_all = get_frame_index(tica_trajs, array)

	assert(len(trajnames) == len(frame_ind_all))
	
	state_name = lim_file.replace(".csv", '').replace(system + "-state-",'')

	#save selected frames indeces for testing
	pickle.dump(frame_ind_all, open(system + '-frames-state%s.pkl'%state_name, 'wb'))

	#for each trajectory, write a txt file containing cpptraj inputs
	for j in range(len(frame_ind_all)):
		frames = frame_ind_all[j]
		if len(frames) != 0: 
			traj_name = trajnames[j]
			f = open(traj_name + '-state-%s.txt'%state_name, 'w')
			f.write('parm *[0-9].p*m*' + '\n')
			for k in range(len(frames)):
				f.write('trajin ' + traj_name + '.xtc' + ' ' + str(frames[k]) + ' ' + str(frames[k]) + ' ' + str(1))
				f.write('\n')
			f.write('trajout ' + traj_name + '-state-%s.xtc'%state_name + '\n')
			f.write('run' + '\n')
			f.write('quit' + '\n')
			f.close()
