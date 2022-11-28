Scripts used for the analysis of randomly selected frames at the energy minima from the total simulation. This includes HOLE radius, LEU197 'pore plug' dihedral, and hydrophobic mistmatch (using Membrainy)

- HOLE radius analysis: 
  - Each trajectory are calculated with *HOLE-compute-radius.py* to obtain a dictionary output that follows the MDAnalysis HOLE package. 
  - To acquire the average HOLE radius at each z position over all trajectories in one macrostate, radii values had to be binned by z coordinate values, as each frame outputs different number of values and relative z positions. This is done with *HOLE-bin-radius.py*
