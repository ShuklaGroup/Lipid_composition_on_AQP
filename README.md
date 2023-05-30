# MD Data Analysis

Scripts for processing and analyzing AMBER .xtc trajectory files

## Features-and-MSM
Contains scripts to perform distance feature selection from trajectory information, then use the final feature set to construct and validate Markov State Models.

## Discrete-traj-analysis
Contains scripts to select 1300-2000 frames from minima of each macrostate in the energy landscape and scripts to perform analyses on them (HOLE radius, pore plug diheral, hydrophobic mismatch).

## Continuous-traj-analysis
Contains scripts for time-dependent analyses, such as water and lipid residence time on a 5000 frame (50 ns) continous trajectory in a energy landscape minima.

# Data Availability

All parameter files and trajectories are organized by protein-structure/lipid and are being uploaded on a Box folder. Smaller files related to each analysis section above are on [Box](https://uofi.box.com/s/uc33gid1jhyuc0oru8tr30to3x9kyj8z) with the following tree

```
.
└── SoPIP2-lipid-GitHub/
    ├── 00-starting-coordinates-and-params/
    │   ├── open-structure/
    │   │   ├── *.rst/ncrst
    │   │   └── *.parm7/prmtop
    │   └── closed-structure
    ├── 01-features-and-MSM/
    │   ├── features-per-system
    │   └── MSM-related-objs/
    │       ├── POPC/
    │       │   ├── *features.pkl: indices of protein residue pairs
    │       │   ├── *tica_obj.pkl: PyEMMA tICA object
    │       │   ├── *tica_trajs.pkl: PyEMMA tICs for each traj and frame
    │       │   ├── *cluster_obj.pkl: PyEMMA k-means clustering object
    │       │   ├── *msm_obj.pkl: PyEMMA MSM object
    │       │   ├── *ITS-error.pkl: PyEMMA implied timescale object
    │       │   ├── *weights.pkl: PyEMMA stationary distribution of MSM
    │       │   └── *probabilities.npy: msm and raw probability of each cluster
    │       ├── POPE
    │       ├── POPG
    │       └── ...
    ├── 02-discrete-frames-and-analysis/
    │   ├── minima-box (*.csv files)
    │   ├── trajectories
    │   ├── dihedral
    │   └── mismatch
    └── 03-continuous-trajs-and-analysis/
        ├── cont-traj-data.csv (all cont trajs information and calculated results)
        ├── trajectories
        ├── water
        └── lipid-order
