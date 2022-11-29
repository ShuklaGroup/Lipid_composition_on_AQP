import glob
import os

#Generate bash scripts for Membrainy to extract .gro files
#containing the annular shell of each frame
for xtcfile in glob.glob('*MD*[0-9]_wat.xtc'):
	name = xtcfile.replace('.xtc','')
	outname = 'shell-trajectory.agr'
	print('name', name)
	b = os.path.isfile(outname)
	if b == False:
		os.makedirs('./output/',exist_ok=True)
		outdir = "./output/"

		filename = 'membrainy-shell-' + name + '.sh'
		print(filename)
		f = open(filename, 'w')
		if ('2b5f' in name):
			top = glob.glob('*2b5f*wat.gro')[0]
		if '1z98' in name:
			top = glob.glob('*1z98*wat.gro')[0]
		f.write('/PATH/jdk-11.0.15/bin/java -jar /PATH/Membrainy-2021.2.jar')
		f.write(' -f ' + xtcfile + ' -s ' + top + ' -ldm vector' + ' -l /PATH/library_tutorial.mbr' + ' -ff CHARMM')
		f.write(' -nt 16' + ' -shell -traj -analysis' + '\n')
		f.write('mv shell-trajectory.gro %s'%outdir + name + '-shell-trajectory.gro' + '\n')
