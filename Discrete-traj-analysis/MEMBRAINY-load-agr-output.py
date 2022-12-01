import glob
import pickle

def GetThickness(agr):
	'''
	Gather thickness information from .agr files (Membrainy ouput) into
	a dictionary
	Args:
		agr (str): the name of the .agr file to be extracted information
	Returns:
		thickness (dict) saved as .pkl file: {keys: frame, values: thickness} 
			with length equals the number of frames 
	'''
	thickness = {}
	f = open(agr,'r')
	name = agr.replace('.agr','')

	#loop through each line of file
	for line in f:

		#all lines not containing values of interest starts with '@'
		if ('@' not in line):
			time, val = line.split()
			thickness[float(time)] = float(val)
	
	f.close()
	#save dictionary
	outname = agr.replace('.agr','.pkl')
	thickness_file = open(outname,'wb')
	pickle.dump(thickness, thickness_file)
	thickness_file.close()

if __name__ == '__main__':
  
  #process file name as argument
	import argparse
	ps = argparse.ArgumentParser()
	ps.add_argument('filename', type = str, help = '.agr file input')
	args = ps.parse_args()
	filename = args.filename
	GetThickness(filename)
