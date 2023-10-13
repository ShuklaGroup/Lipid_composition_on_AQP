import glob
import numpy as np
import pylipid
from pylipid.api import LipidInteraction

def GetLipResTime(filename):
    '''
    Compute the lipid residence time at each residue for continuous input trajectory
    Args:
        filename (str): file name of the .xtc file with only protein and lipid
    Returns:
        residence_time (np.array) of shape (236, ) containing the lipid residence
            time at each residue
        pbd files of the protein in the output directory "pylipid-cont" with the residence time
            as b-factor values
    '''
    trajfile_list = [filename]
    if "2b5f" in filename:
        topfile_list = [glob.glob("*2b5f*_wat.gro")[0]]
    if "1z98" in filename:
        topfile_list = [glob.glob("*1z98*_wat.gro")[0]]
    if "complex" in filename:
        lipid = "complex_A"
    else:
        lipid = filename[:4]

    #define pylipid parameters and output directory for the pdbs
    cutoffs = [0.4, 0.8]
    nprot = 1
    timeunit = "ns"
    save_dir = "pylipid-cont"

    li = LipidInteraction(trajfile_list, topfile_list=topfile_list, cutoffs=cutoffs, lipid=lipid, timeunit=timeunit,
                        nprot=nprot, save_dir=save_dir)

    #perform necessary calculations to acquire lipid res time & binding sites
    li.collect_residue_contacts()
    koffs, res_times = li.compute_residue_koff(fig_close=True)
    node_list, modularity = li.compute_binding_nodes(threshold=4, print_data=True)

    #li.dataset is a dataframe including calculated parameters for each residue,
    #fetch the residence time as an array from this dataframe
    residence_time = list(li.dataset["Residence Time"])
    outname = filename.replace('_wat.xtc','')[-40:]
    np.save(outname + "-pylipid-restime.npy",residence_time)
    li.save_coordinate(item="Residence Time",fn_coord=outname+'_restime.pdb')

    #get df for binding sites as well
    num_binding_sites = len(node_list)

    for i in range(num_binding_sites):
        df = li.dataset[li.dataset["Binding Site ID"] == i]
        df.to_csv(outname + '-binding-site-index-' + str(i) + '.csv')

if __name__ == '__main__':

	import argparse
    	#fetch user input argument as the .xtc filename
	ps = argparse.ArgumentParser()
	ps.add_argument('filename', type = str, help = '_wat.xtc file input with only protein and lipid')
	args = ps.parse_args()
	filename = args.filename
	GetLipResTime(filename)
