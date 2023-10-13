Scripts for the time-dependent analysis of continuous trajectories occupying the minima of the energy landscape, with water count for each 100-ns macrostate trajectory, water residence time per lateral pore slice, and lipid order parameter. Python 3.7.12 environment used is provided in [watlip-env.yml](watlip-env.yml). Dependencies include MDAnalysis 2.1.0 and LiPyphilic 0.10.0.

- [*count-water.py*](count-water.py) takes an argument of the \*.xtc file in the directory, save the import and export water information for that trajectory (Figure 5 of manuscript)
 
- [*lipid-ord-param-lipyphilic.py*](lipid-ord-param-lipyphilic.py) also takes an argument of the \*.xtc file in the trajectory corresponding to only homogeneous lipid systems, uses LiPyphilic to get an average lipid order parameter through the trajectory of each lipid's acyl tails. The acquired data is saved and also plotted (Figure 7 of manuscript)

- [*water-residence-time.py*](water-residence-time.py) acquires residence time of each water in each pore lateral slice. (Figure S20)

- HOLE radius analysis: 
  - To obtain the HOLE radius of the protein in each frame of each trajectory, [*HOLE-compute-radius.py*](HOLE-compute-radius.py) is used to obtain a dictionary output that follows the MDAnalysis HOLE package (data structure is explained in details in the code). 
  - To acquire the average HOLE radius at each z position over all trajectories in one macrostate, radii values had to be binned by z coordinate values, as each frame outputs different number of values and relative z positions. This is done with [*HOLE-bin-radius.py*](HOLE-bin-radius.py)

Lipid binding analyses were also performed on continuous trajectories as part of the peer review process. These analyses were done using PyLipID. The installation for PyLipID is a bit tricky. We provide [anh_pylipid.yml](anh_pylipid.yml) to help with installation. We recommend installing each of dependencies in [anh_pylipid_no_pip.yml](anh_pylipid_no_pip.yml) line-by-line first using `conda`, and then installing the remaining dependencies lline-by-line using `pip`. PyLipID also constitutes a Python 3.7.12 environment. Once a pylipid environment has been created, the following script can be used:

- [*lipid-residence-time.py*](lipid-residence-time.py) takes an argument of the \*.xtc file in the directory. This file return a .npy array of lipid residence time per residue. $lipid = any lipid name i.e. "POPG" for homogeneous system. If running for a complex bilayer, change just the `lipid` field value within the `LipidInteraction` function to be the particular lipid residue name you are interested in studying. (Supplementary heatmaps from revisions of manuscript)
