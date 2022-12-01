import numpy as np
import pickle
import glob
import os
import pyemma
import time
import matplotlib.pyplot as plt

def GridSearch(system, dim, cluster):
	'''
	Perform one step of grid searching on a set of specified hyperparameters
	Args:
		system (str): name of bilayer system, extracted from feature file names
		dim (int): number of tIC
		cluster (int): number of clusters
	Returns:
		system_name (str): name of bilayer, along with tested hyperparameters
		vamp_avg (float): VAMP-2 score of current hyperparameter set
	'''
	
	#set up system name for file saving
	system_name = system + '-' + str(dim) + 'dim-' + str(cluster) + 'cluster-'

	#tICA
	tica_obj = pyemma.coordinates.tica(data = feat_all, dim = dim) #lag should be 10-30% of msm lag time
	tica_trajs = tica_obj.get_output()
	tica_file = open(outdir + system_name + 'tica_obj' + '.pkl','wb')
	pickle.dump(tica_obj, tica_file)
	tica_file.close()
	tica_traj_file = open(outdir + system_name + 'tica_trajs' + '.pkl','wb')
	pickle.dump(tica_trajs, tica_traj_file)
	tica_traj_file.close()

	#Clustering
	localtime = time.asctime(time.localtime(time.time()))
	print("Start clustering", system_name, localtime)

	cluster_obj = pyemma.coordinates.cluster_mini_batch_kmeans(tica_obj, k = cluster, max_iter = 200, stride = 5)
	#For large systems, we recommend to pass the tica object itself into the subsequent stages, 
	#e.g., clustering, in order to avoid loading all transformed data into memory.
	cluster_trajs = cluster_obj.dtrajs
	cluster_file = open(outdir + system_name + 'cluster_obj' + '.pkl','wb')
	pickle.dump(cluster_obj, cluster_file)
	cluster_file.close()

	#Implied timescale plot

	its = pyemma.msm.its(cluster_obj.dtrajs, lags = [50, 100, 200, 400, 500, 700, 800], nits = 5)
	its_file = open(outdir + system_name + 'its_obj' + '.pkl','wb')
	pickle.dump(its, its_file)
	its_file.close()
	pyemma.plots.plot_implied_timescales(its)
	plt.savefig(outdir + system_name + 'its.png')
	plt.clf()

	#MSM construction and validation
	msm = pyemma.msm.estimate_markov_model(cluster_trajs, lag = 200)
	
	#VAMP-2
	vamp_avg = np.mean(msm.score_cv(cluster_trajs))

	msm_file = open(outdir + system_name + 'msm_obj' + '.pkl','wb')
	pickle.dump(msm, msm_file)
	msm_file.close()

	
	return system_name, vamp_avg


if __name__ == '__main__':
	

	os.makedirs('./GRID/',exist_ok=True)
	outdir = "./GRID/"

	#Read in feature files
	feat_all = []

	localtime = time.asctime(time.localtime(time.time()))
	print("Start loading feature files", localtime)

	for file in sorted(glob.glob('*-FINAL-FEATURES.pkl')):
		feat = pickle.load(open(file, 'rb'))
		feat_all.append(feat)
	
	#get system name
	system = file.replace('_','-').split('-')[0]
	
	#total number of contact features
	num_features = len(feat_all[0][0])

	#iterate through tICA dims, clusters
	#change depends on system
	dims = np.linspace(4,8,5)
	clusters = np.linspace(300,1000,8)
	vamp_scores = {}

	for i, dim in enumerate(dims):
		for j, cluster in enumerate(clusters):

			localtime = time.asctime(time.localtime(time.time()))
			print("Start Grid for %i dim, %i cluster" %(dim, cluster), localtime)
			system_name, vamp_avg = GridSearch(system, int(dim), int(cluster))
			vamp_scores[system_name] = vamp_avg
	
	#save VAMP score dictionary
	with open(outdir + system + '_vamp_dict.pkl', 'wb') as f:
		pickle.dump(vamp_scores, f)
		f.close()
