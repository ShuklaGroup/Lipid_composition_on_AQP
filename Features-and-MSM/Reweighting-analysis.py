import numpy as np
import pickle
import matplotlib.pyplot as plt
import glob
import pyemma

###Obtain cluster object and weights array corresponding to each frame, count the probability of each cluster, then plot. 
###This analysis ensures reweighting is not too diverged from original probabilities 

def GetRawAndReweightProb(cluster_trajs, weights_array, cluster_file):
    '''
    Gets the log10 probability of each cluster, from the cluster object and reweighting list
    Args:
        cluster_trajs (list): list of np.array, 
			each array corresponds to a trajectory with shape (num of frames, )
			and contains the assigned cluster of each frame in the trajectory
        weights_array (list): list of np.array,
			each array corresponds to a trajectory with shape (num of frames, )
			and contains the weights each frame in the trajectory assigned by PyEMMA
		cluster_file (str): .pkl filename of the cluster object from PyEMMA for saving outputs
    Returns:
        log_raw_prob (np.array): shape (num of cluster,), log of the raw probability of each cluster in the cluster object
        log_msm_prob (np.array): shape (num of cluster,), log of the raw probability of each cluster in the reweights
    '''

    
    raw_cluster_conc = np.concatenate(cluster_trajs)
    weights_conc = np.concatenate(weights_array)
    num_clusters = max(raw_cluster_conc) + 1

    raw_counts = np.zeros(num_clusters, dtype=float)
    msm_counts = np.zeros(num_clusters, dtype=float)
    for num in range(num_clusters):
        ind = np.where(raw_cluster_conc == num)[0]
        raw_counts[num] += len(ind)
        msm_counts[num] += np.sum(weights_conc[ind])

    raw_prob = raw_counts/len(raw_cluster_conc)
    log_raw_prob = np.log10(raw_prob)
    log_msm_prob = np.log10(msm_counts)

    system_name = cluster_file.split('-')[0]
    np.save(system_name+'_raw_probabilites.npy',log_raw_prob)
    np.save(system_name+'_msm_probabilites.npy',log_msm_prob)
    return log_raw_prob, log_msm_prob

if __name__ == '__main__':
    
    #find cluster object and weights list in directory
    cluster_file = glob.glob("*cluster_obj.pkl")[0]
    cluster_obj = pickle.load(open(cluster_file, 'rb'))
    cluster_trajs = cluster_obj.dtrajs
    	
    weights_file = glob.glob("*weights-final*pkl")[0]
    weights_array = pickle.load(open(weights_file, 'rb'))

    system_name = cluster_file.split('-')[0]
    log_raw_prob, log_msm_prob = GetRawAndReweightProb(cluster_file, weights_file)
    plt.scatter(log_raw_prob, log_msm_prob)
    
    #plot x=y line
    xmin = min(log_raw_prob)
    xmax = max(log_raw_prob)
    x = np.linspace(xmin-1,xmax+1,100)
    plt.plot(x,x,'k-')
    plt.xlabel('Raw Counts Probability ($log_{10}$)')
    plt.ylabel('MSM Counts Probability ($log_{10}$)')
    plt.title('SoPIP2:%s MSM Reweighting' %system_name)
    plt.savefig(system_name + '-reweighting-plot.png',dpi=600)
    plt.clf()
