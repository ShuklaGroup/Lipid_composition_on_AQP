import mdtraj as md
import numpy as np
import pickle
import matplotlib
import matplotlib.pyplot as plt
matplotlib.use("Agg")
import glob
import pyemma
import time
from pyemma.coordinates import tica

#Function that performs oASIS, returns column index
def run_oasis(input_feature_data, num_features,
	      lag = 4, max_columns = 50, spectraldir = "."):
	'''
	Args: input_feature_data: numpy array of all feature data
		lag: (int) tica lag time
		max_columns: (int) number of desired final features
		num_features: (int) number of total features in input_feature_data
	Returns: col_ind (array of int): column indices of selected features
	'''
	t_obj = pyemma.coordinates.tica_nystroem(data = input_feature_data,
						lag = lag,
						max_columns = max_columns,
						initial_columns = np.random.choice(num_features,1,replace=False))

	t_file = open("%s/nystroem-tica-obj-%d.pkl"%(spectraldir, max_columns),'wb')
	pickle.dump(t_obj, t_file)
	t_file.close()

	#get final column indeces of the chosen features, sorted
	col_ind = sorted(t_obj.column_indices)

	col_file = open("%s/nystroem-col-obj-%d.pkl"%(spectraldir, max_columns),'wb')
	pickle.dump(col_ind, col_file)
	col_file.close()
    
	return col_ind

if __name__ == "__main__":
	#combine all post-RRCS features
	feat_all = []

	for file in sorted(glob.glob('*-RRCS-FEATURES.pkl')):
		feat = pickle.load(open(file, 'rb'))
		feat_all.append(feat) #Take only the pair-wise distance arrays

	#get rrcs residue pairs
	rrcs_res_pair = np.load('post-RRCS-res-index.npy')

	#Total number of contact features
	num_features = len(feat_all[0][0])

	#Set up parameters for optimization
	#lag = np.array([4,5])
	cols = np.arange(10,105,5)

	os.makedirs('./oASIS-out/',exist_ok=True)
	outdir = "./oASIS-out/"
	vampfile = open('oASIS-vamp.txt','w')
	vampfile.close()
	
	for i in range(len(cols)):
		localtime = time.asctime(time.localtime(time.time()))
		print("Start oASIS round %i"%i, localtime)
		col_ind = run_oasis(feat_all,
					#lag = lagtime,#
					max_columns = cols[i],
					num_features = num_features)
		localtime = time.asctime(time.localtime(time.time()))
		print("Done oASIS round %i"%i, localtime)

		#visualize full tica with given features and lagtime
		feat = []
		for traj in feat_all:
			feat.append(traj[:,col_ind])
		tica_obj = tica(feat)
		tica_trajs = tica_obj.get_output()
		tica_concatenated = np.concatenate(tica_trajs)
		plt.hexbin(tica_concatenated[:,0].T, tica_concatenated[:,1].T, 
					gridsize=300,bins='log',cmap='jet',mincnt=1)
		plt.savefig("tica-visualization-full-%i-feature-post-oASIS.png"%(cols[i]))
		plt.clf()

		#get residue pairs that passed oASIS
		oasis_res_pair = []
		oasis_res_pair = rrcs_res_pair[col_ind]
		oasis_res_file = open("oasis-res-pairs-%i-feature.pkl"%(cols[i]), "wb")
		pickle.dump(oasis_res_pair, oasis_res_file)
		oasis_res_file.close()
		
		#VAMP-score optimization
		cluster_obj = pyemma.coordinates.cluster_mini_batch_kmeans(tica_obj, k = 500, max_iter = 200, stride = 5)
		cluster_trajs = cluster_obj.dtrajs
		msm = pyemma.msm.estimate_markov_model(cluster_trajs, lag = 200)
		vamp_avg = np.mean(msm.score_cv(cluster_trajs))
		vamp_scores[cols[i]] = vamp_avg
		
		vampfile = open('oASIS-vamp.txt','a')
		vampfile.write(str(cols[i]) + ': ' + str(vamp_avg) + '\n')
		vampfile.close()
