import numpy as np
import pickle
import matplotlib.pyplot as plt
import glob
import pyemma

###Obtain cluster object and weights array corresponding to each frame, count the probability of each cluster, then plot. 
###This analysis ensures reweighting is not too diverged from original probabilities 

def GetRawAndReweightProb(cluster_file, weights_file):
    '''
    Gets the log10 probability of each cluster, from the cluster object and reweighting list
    Args:
        cluster_file (str): filename of the cluster object outputted from PyEMMA.
        weights_file (str): filename of the reweighting obtained from the MSM built with PyEMMA
    Returns:
        log_raw_prob (np.array): log of the raw probability of each cluster in the cluster object, shape (num of cluster,)
        log_msm_prob (np.array): log of the raw probability of each cluster in the reweights, shape (num of cluster,)
    '''
    cluster_obj = pickle.load(open(cluster_file, 'rb'))
    cluster_trajs = cluster_obj.dtrajs
    raw_cluster_array = cluster_trajs
    
    weights_array = pickle.load(open(weights_file, 'rb'))
    
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
    return log_raw_prob, log_msm_prob

if __name__ == '__main__':
    
    #find cluster object and weights object in directory
    cluster_file = glob.glob("*cluster_obj.pkl")[0]
    weights_file = glob.glob("*weights-final*pkl")[0]
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
