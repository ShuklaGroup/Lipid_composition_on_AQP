import numpy as np
import pickle
import matplotlib.pyplot as plt
import glob
import pyemma

###Obtain cluster object and weights array corresponding to each frame, count the probability of each cluster, then plot. 
###This analysis ensures reweighting is not too diverged from original probabilities 

cluster_file = glob.glob("*cluster_obj.pkl")[0]
cluster_obj = pickle.load(open(cluster_file, 'rb'))
cluster_trajs = cluster_obj.dtrajs

weights_file = glob.glob("*weights-final*pkl")[0]
weights_array = pickle.load(open(weights_file, 'rb'))
raw_cluster_array = cluster_trajs

raw_array = np.concatenate(raw_cluster_array)
weighted_array = np.concatenate(weights_array)
num_clusters = max(raw_array) + 1

raw_counts = np.zeros(num_clusters, dtype=float)
msm_counts = np.zeros(num_clusters, dtype=float)
for num in range(num_clusters):
    ind = np.where(raw_array == num)[0]
    raw_counts[num] += len(ind)
    msm_counts[num] += np.sum(weighted_array[ind])

raw_prob = raw_counts/len(raw_array)
log_raw_prob = np.log10(raw_prob)
log_msm_prob = np.log10(msm_counts)

system_name = cluster_file.split('-')[0]
np.save(system_name+'_raw_probabilites.npy',log_raw_prob)
np.save(system_name+'_msm_probabilites.npy',log_msm_prob)

xmin = min(log_raw_prob)
xmax = max(log_raw_prob)
x = np.linspace(xmin-1,xmax+1,100)
plt.plot(x,x,'k-')
plt.scatter(log_raw_prob, log_msm_prob)
plt.xlabel('Raw Counts Probability ($log_{10}$)')
plt.ylabel('MSM Counts Probability ($log_{10}$)')
plt.title('SoPIP2:%s MSM Reweighting' %system_name)
plt.savefig(system_name + '-reweighting-plot.png',dpi=600)
plt.clf()
