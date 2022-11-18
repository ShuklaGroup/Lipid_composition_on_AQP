import MDAnalysis
import matplotlib.pyplot as plt
import numpy as np
import numpy.linalg as la
import time
import pickle
import glob

def GetWatsInFrame(filename, slice, slicenum):
	'''
	Args: 
		a xtc file containing protein and water, type 'str'
		a slice indication in the pore, type 'str'
	Return: 
		list of np.arrays, each corresponds to all wat resid identified in the slice in a frame
	'''
	if '2b5f' in filename:
		pdbtraj = glob.glob("*2b5f*stripped.pdb")[0]
	if '1z98' in filename:
		pdbtraj = glob.glob("*1z98*stripped.pdb")[0]
	trajectorytraj = filename
	universetraj = MDAnalysis.Universe(pdbtraj, trajectorytraj)
	wat = universetraj.select_atoms("resname WAT and name O and cyzone 8 4 -4 %s" %slice, updating=True)
	allwat_allframes = []
	allwat_allframes.append(wat.resids)
	for i in range(len(universetraj.trajectory)-1):
		universetraj.trajectory.next()
		allwat_allframes.append(wat.resids)
	outfile = filename.replace('.xtc', '-WAT-PER-FRAME-SLICE-%i.pkl'%(slicenum+1))
	file = open(outfile, "wb")
	pickle.dump(allwat_allframes, file)
	file.close()
	return allwat_allframes

def GetFramesofWat(filename, allwat_allframes, slicenum):
	'''
	From the list of np.arrays with all water resids occupying slice at each frame,
	get another list of np.arrays, each array corresponds to a water resid, containing
	the frames in which it occupies the slice
	Args: 'str' filename (to save files), list of np.arrays
	Return: list of np.arrays of frames per wat id, list of all wats
	'''
	#get all water ids that have been in slice
	allwatset = list(sorted(set(np.concatenate(allwat_allframes))))
	frames_of_wat_arr = []
	#for each water id, record which frame it's in the slice
	for i in range(len(allwatset)):
		frames = []
		for j in range(len(allwat_allframes)):
			if allwatset[i] in allwat_allframes[j]:
				frames.append(j)
		frames_of_wat_arr.append(np.array(frames))
	#save
	outfile = filename.replace('.xtc', '-FRAMES-PER-WAT-SLICE-%i.pkl'%(slicenum+1))
	file = open(outfile, "wb")
	pickle.dump(frames_of_wat_arr, file)
	file.close()

	allwatoutfile = filename.replace('.xtc', '-ALL-WAT-SLICE-%i.pkl'%(slicenum+1))
	allwatfile = open(allwatoutfile, "wb")
	pickle.dump(allwatset, allwatfile)
	allwatfile.close()

	return frames_of_wat_arr

def GetResTime(filename, frames_of_wat_arr, slicenum):
	'''
	From the list of np.arrays with each water and the frames in which it's in the slice,
	calculate how many frames this water is in the slice continuously
	Args: 'str' filename (to save files), list of np.arrays, list of ints
	Return: list of np.arrays
	'''
	all_wat_count_frame_num = []
	for wat in range(len(frames_of_wat_arr)):
		count = 0
		count_frame_num = []
		prev = frames_of_wat_arr[wat][0]
		for i in range(len(frames_of_wat_arr[wat])):
			curr = frames_of_wat_arr[wat][i]
			if curr == prev + 1:
				count += 1
			else:
				count_frame_num.append(count)
				count = 0
			prev = curr
		all_wat_count_frame_num.append(np.array(count_frame_num))

	avg = []
	for wat in all_wat_count_frame_num:
		avg.append(np.mean(wat))
	
	std = []
	for wat in all_wat_count_frame_num:
		std.append(np.std(wat))

	outfile = filename.replace('.xtc', '-AVG-RES-TIME-SLICE-%i.pkl'%(slicenum+1))
	file = open(outfile, "wb")
	pickle.dump(avg, file)
	file.close()

	stdoutfile = filename.replace('.xtc', '-STD-RES-TIME-SLICE-%i.pkl'%(slicenum+1))
	stdfile = open(stdoutfile, "wb")
	pickle.dump(std, stdfile)
	stdfile.close()

	watoutfile = filename.replace('.xtc', '-ALL-WAT-COUNT-SLICE-%i.pkl'%(slicenum+1))
	watfile = open(watoutfile, "wb")
	pickle.dump(all_wat_count_frame_num, watfile)
	watfile.close()

if __name__ == '__main__':
	import argparse

	ps = argparse.ArgumentParser()
	ps.add_argument('filename', type = str, help = '-stripped-state-*.xtc 5000 frame file input with protein and water')
	args = ps.parse_args()
	filename = args.filename
	slices = ["resid 47 or resid 206 or resid 211", 
				"resid 32 or resid 137 or resid 187 or resid 202",
				"resid 29 or resid 107 or resid 218",
				"resid 54 or resid 141 or resid 183 or resid 199",
				"resid 21 or resid 58 or resid 222",
				"resid 17 or resid 96 or resid 175 or resid 226",
				"resid 13 or resid 81 or resid 152 or resid 230",
				"resid 10 or resid 90 or resid 233"]
	for i in range(8):
		start_time = time.time()
		allwat_allframes = GetWatsInFrame(filename, slices[i], i)
		frames_of_wat_arr = GetFramesofWat(filename, allwat_allframes, i)
		GetResTime(filename, frames_of_wat_arr, i)
		print("--- %s seconds ---" % (time.time() - start_time))
