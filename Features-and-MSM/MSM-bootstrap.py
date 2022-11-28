import glob
import pickle
import pyemma
import numpy as np
import os

def BootStrap(tica_file, cluster_file, msmlag, edges, nbins=175, nsamples=200):
	'''
	Perform bootstrapping on tICA output to evaluate relative errors

	Args: 
		tica_file (str): .pkl filename of the tICA output from pyEMMA tICA object
			contains a list of np.array, each array corresponds to a trajectory 
				and has shape (num of frames, num of tICA dimension)
		cluster_file (str): .pkl filename of the cluster object from pyEMMA
		msmlag (int): chosen lagtime of MSM, 
			typically the lagtime at which the implied timescale plot has converged
		edges (tuple): shape (4), containing the total data edges
		nbins (int): number of bins to assign data to
		nsamples (int): number of subsets of data to randomly select

	Returns:
		max_prob (np.array): shape (200, ), contains the maximum probability of each subset
		prob_arr (np.array): shape (200, nbins, nbins)
	'''
  	#get the tica outputs and cluster object from MSM
	tica_trajs = pickle.load(open(tica_file, "rb"))
	tica_concatenate = np.concatenate(tica_trajs)

	cluster_obj = pickle.load(open(cluster_file, "rb"))
	cluster_traj = cluster_obj.dtrajs

	#system
	system_name = tica_file.split('-')[0]
	trajnum = len(tica_trajs)
  
	#create output directory
	os.makedirs('./bootstrap-test/',exist_ok=True)
	outdir = "./bootstrap-test/"
	x_bins = nbins
	y_bins = nbins

	#Initialize arrays
	prob_arr = np.zeros([nsamples,x_bins,y_bins])
	max_prob = np.zeros(nsamples)
	
	for i in range(nsamples):
		selected_cluster = []
		selected_tica_x = []
		selected_tica_y = []
		#randomly sample trajectory indeces for 80% of data
		selected_ind = np.random.choice(np.arange(0,trajnum,1), size=int(trajnum*0.8), replace=False)
		for ind in selected_ind:
			selected_cluster.append(cluster_traj[ind])
			selected_tica_x.append(tica_trajs[ind][:,0])
			selected_tica_y.append(tica_trajs[ind][:,1])
		#make msm from data subset
		msm = pyemma.msm.estimate_markov_model(selected_cluster, lag = msmlag)	
		msm_file = open(outdir + system_name + '-msm_obj-bt-%i'%i + '.pkl','wb')
		pickle.dump(msm, msm_file)
		msm_file.close()
		#get weights from msm
		weights = msm.trajectory_weights()
		weight_file = open(outdir + system_name + '-weights-bt-%i'%i + '.pkl','wb')
		pickle.dump(weights, weight_file)
		weight_file.close()
		weights_conc = np.concatenate(weights)
		
		#concatenate the 80% data subset and total dataset
		x_80 = np.concatenate(selected_tica_x)
		y_80 = np.concatenate(selected_tica_y)
		
		#Define the bin edges with the maximum and minimum from the total data
		x_lower_bound, y_lower_bound, x_upper_bound, y_upper_bound = edges

		hist= np.histogram2d(x_80,y_80,bins=175,
				range = [[x_lower_bound,x_upper_bound],[y_lower_bound,y_upper_bound]],
				density= True, weights = weights_conc)
		
		#collect the probability density
		prob_density = hist[0]
		assert(prob_density.shape == (x_bins, y_bins))
		max_prob[i] = np.max(prob_density)
		prob_arr[i,:,:] = prob_density

	#save the probability density and max probability density
	prob_file = open(outdir + system_name + '-prob-bt.pkl','wb')
	pickle.dump(prob_arr, prob_file)
	prob_file.close()

	maxprob_file = open(outdir + system_name + '-maxprob-bt.pkl','wb')
	pickle.dump(max_prob, maxprob_file)
	maxprob_file.close()

	return max_prob, prob_arr

def GetErr(tica_file, max_prob, prob_arr):
	'''
	Calculate the free energy from the binned probability density and plot it
	
	Args:
		tica_file (str): .pkl filename of the tICA output from pyEMMA tICA object
			contains a list of np.array, each array corresponds to a trajectory 
				and has shape (num of frames, num of tICA dimension)
		max_prob (np.array): shape (nsamples, ) containing the max probability acquired
			from each randomly selected sample
		prob_arr (np.array): shape (nsamples, nbins, nbins) containing the probability
			of each bins computed from bootstrapping
	
	Returns:
		err_free_energy (np.array): shape (nbins, nbins) 
			as a saved .pkl file with the standard error in free energy computed from
			max_prob and prob_arr

	'''
	#calculation parameters
	R = 0.001987
	T = 300

	mean_prob_density = np.mean(prob_arr,axis=0)
	err_prob_density = np.std(prob_arr,axis=0)

	err_free_energy = R*T*np.abs(np.divide(err_prob_density,mean_prob_density) - np.std(max_prob)/np.mean(max_prob))

	#save figure
	system_name = tica_file.split('-')[0]
	errfile = open('FINAL-'+system_name+'_tic-err.pkl','wb')
	pickle.dump(err_free_energy, errfile)
	errfile.close()

if __name__ == '__main__':

	import argparse
	#input the MSM lagtime as a script argument
	ps = argparse.ArgumentParser()
	ps.add_argument('msmlag', type = int, help = 'MSM lagtime')
	args = ps.parse_args()
	msmlag = args.msmlag

	#find tica output, cluster object, and weights list files
	tica_file = glob.glob("*tica_trajs.pkl")[0]
	cluster_file = glob.glob("*cluster_obj.pkl")[0]
	weights_file = glob.glob("*weights-final.pkl")[0]

	tica_trajs = pickle.load(open(tica_file, "rb"))
	tica_concatenate = np.concatenate(tica_trajs)

	#from first 2 tICA component, define bin edges
	x_data = tica_concatenate[:,0].T
	y_data = tica_concatenate[:,1].T

	edges = (np.min(x_data) - 0.5, np.min(y_data) - 0.5,
			np.max(x_data) + 0.5 , np.max(y_data) + 0.5)

	#perform boostrapping and get probability of each frame
	max_prob, prob_arr = BootStrap(tica_file, cluster_file, int(msmlag), edges=edges,
									nbins=175, nsamples=200)
	
	#calculate the error in free energy
	GetErr(tica_file, max_prob, prob_arr)
