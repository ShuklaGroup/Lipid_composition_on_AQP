import numpy as np
import pickle
import matplotlib
import matplotlib.pyplot as plt
matplotlib.use("Agg")
import glob
# import seaborn as sns
# sns.set_style("white")
from scipy import stats

feat_all_o = []
for file in sorted(glob.glob('*MD*-OPEN-ABS-SUM-RRCS.pkl')):

	#out_file = file.replace('-CONTACT.pkl','')

	feat_file_o = pickle.load(open(file, 'rb'))
	feat_all_o.append(feat_file_o)

feat_all_o = np.array(feat_all_o)
feat_all_c = []

for file in sorted(glob.glob('*MD*-CLOSED-ABS-SUM-RRCS.pkl')):

	#out_file = file.replace('-CONTACT.pkl','')

	feat_file_c = pickle.load(open(file, 'rb'))
	feat_all_c.append(feat_file_c)

feat_all_c = np.array(feat_all_c)

open_rrcs = np.sum(feat_all_o,axis=0)
closed_rrcs = np.sum(feat_all_c,axis=0)

open_rrcs_log = np.log(open_rrcs/10**6)
closed_rrcs_log = np.log(closed_rrcs/10**6)

zscores = stats.zscore(open_rrcs_log)
selected_dist_ind_o = []
for i in range(len(zscores)):
    if zscores[i] >= 1.5 or zscores[i] <= -1.5:
        selected_dist_ind_o.append(i)

minnumo = np.min(open_rrcs_log)
maxnumo = np.max(open_rrcs_log)
#getting the max/min values of -1.5<z-score<1.5 region
for i in selected_dist_ind_o:
    if open_rrcs_log[i] > 0:
        if open_rrcs_log[i] < maxnumo:
            maxnumo = open_rrcs_log[i]
    if open_rrcs_log[i] < 0:
        if open_rrcs_log[i] > minnumo:
            minnumo = open_rrcs_log[i]

zscores = stats.zscore(closed_rrcs_log)
selected_dist_ind_c = []
for i in range(len(zscores)):
    if zscores[i] >= 1.5 or zscores[i] <= -1.5:
        selected_dist_ind_c.append(i)
minnumc = np.min(closed_rrcs_log)
maxnumc = np.max(closed_rrcs_log)
for i in selected_dist_ind_c:
    if closed_rrcs_log[i] > 0:
        if open_rrcs_log[i] < maxnumc:
            maxnumc = closed_rrcs_log[i]
    if closed_rrcs_log[i] < 0:
        if open_rrcs_log[i] > minnumc:
            minnumc = closed_rrcs_log[i]

final_ind = np.union1d(sorted(selected_dist_ind_c),sorted(selected_dist_ind_o))
np.save('post-RRCS-features-ind.npy', final_ind)
