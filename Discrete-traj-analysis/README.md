Scripts used for the analysis of randomly selected frames at the energy minima from the total simulation. This includes HOLE radius, LEU197 'pore plug' dihedral, and hydrophobic mistmatch (using Membrainy)

- [*get-cpptraj-script-minima.py*](get-cpptraj-script-minima.py) is used to generate CPPTRAJ scripts to extract frames selected in a minima of the tICA free energy landscape. Each minima is defined by a box with max/min of width/height in \*.csv files on Box (path: SoPIP2-lipid-GitHub/02-discrete-frames-and-analysis/minima-box/). Outputs of this script are \*.txt files, each file corresponds to a trajectory with frames inside indicated box.

- HOLE radius analysis: 
  - To obtain the HOLE radius of the protein in each frame of each trajectory, [*HOLE-compute-radius.py*](HOLE-compute-radius.py) is used to obtain a dictionary output that follows the MDAnalysis HOLE package (data structure is explained in details in the code). 
  - To acquire the average HOLE radius at each z position over all trajectories in one macrostate, radii values had to be binned by z coordinate values, as each frame outputs different number of values and relative z positions. This is done with [*HOLE-bin-radius.py*](HOLE-bin-radius.py)

- Hydrophobic mismatch:
  - With [Membrainy](http://www.membrainy.net/) and dependencies (Java, jdk) installed, generate bash scripts per trajectory to extract annular shell as a .gro file (or thickness, entropy by changing an argument) using [*MEMBRAINY-make-shell-scripts.py*](MEMBRAINY-make-shell-scripts.py)
  - Membrainy outputs .agr files containing calculated thickness. Each of these files are saved as .pkl dictionary with [*MEMBRAINY-load-agr-output.py*](MEMBRAINY-load-agr-output.py)
