# MD Data Analysis

Scripts for processing and analyzing AMBER .xtc trajectory files

## Features-and-MSM
Contains scripts to perform distance feature selection from trajectory information, then use the final feature set to construct and validate Markov State Models.

## Discrete-traj-analysis
Contains scripts to select 1000 frames from minima of each macrostate in the energy landscape and scripts to perform analyses on them (HOLE radius, pore plug diheral, hydrophobic mismatch).

## Continuous-traj-analysis
Contains scripts for time-dependent analyses, such as water count, water residence time, and lipid order parameter on a 10000 frame (100 ns) continous trajectory in a energy landscape minima.

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
    │       └── ... (other systems)
    ├── 02-discrete-frames-and-analysis/
    │   ├── minima-box (*.csv files)
    │   ├── trajectories (*.xtc files with water stripped along with *gro parameters)/
    │   │   ├── POPC-frames-stripped-wat.tar.bz2
    │   │   └── ... (other systems)
    │   ├── dihedral/
    │   │   ├── POPC-dihedral.tar.bz2 (all *npy containing dihedral data of POPC system)
    │   │   ├── POPE-dihedral.tar.bz2
    │   │   ├── POPG-dihedral.tar.bz2
    │   │   └── ... (other systems)
    │   └── mismatch/
    │       ├── prot-bulk-mismatch.tar.bz2 (all *pkl containing prot-bulk data)
    │       └── shell-bulk-mismatch.tar.bz2 (all *pkl containing shell-bulk data)
    └── 03-continuous-trajs-and-analysis/
        ├── cont-traj-data.csv (all cont trajs information and calculated results)
        ├── trajectories-stripped-lipid/
        │   ├── POPC-traj-stripped-lipid.tar.bz2/
        │   │   ├── *stripped*xtc of 3 trajs (lipids stripped) per macrostate
        │   │   └── *pdb of coordinate files for analyses
        │   ├── POPE-traj-stripped-lipid.tar.bz2
        │   ├── POPG-traj-stripped-lipid.tar.bz2
        │   └── ... (other systems)
        ├── trajectories-stripped-water/
        │   ├── POPC-traj-stripped-water.tar.bz2/
        │   │   ├── *wat.xtc of 3 trajs (waters stripped) per macrostate
        │   │   └── *pdb of coordinate files for analyses
        │   ├── POPE-traj-stripped-water.tar.bz2
        │   ├── POPG-traj-stripped-water.tar.bz2
        │   └── ... (other systems)
        ├── water/
        │   ├── POPC-wat-transport.tar.bz2/
        │   │   ├── *passagetime*: time it takes for each water to transport
        │   │   ├── *transportedres*: residue ID of the waters that was transported
        │   │   ├── *transport*: 3d array recording whether a water was transported at each frame
        │   │   └── *wat-in-pore*: number of water occupying the pore at each frame
        │   ├── POPC-wat-restime.tar.bz2/
        │   │   ├── *AVG-RES-TIME*: average time water continuously spends in slice, per water
        │   │   ├── *STD-RES-TIME*: standard deviation of average residence time
        │   │   └── *FRAMES-PER-WAT*: per water, record frames during which water occupy slice
        │   │   ├── *ALL-WAT-COUNT*: number of water continuously in slice, per water
        │   │   ├── *ALL-WAT*: record all water in slice per frame
        │   └── ... (other systems)
        └── lipid-order/
            ├── POPC-lipid-ord.tar.bz2/
            │   ├── *protein-x.pkl: protein x positions across traj
            │   ├── *protein-y.pkl: protein y positions across traj
            │   ├── *lipid-x.pkl: lipid x positions of each lipid
            │   ├── *lipid-y.pkl: lipid y positions of each lipid
            │   ├── *scc_sn1.pkl: order param of each lipid's tail sn1
            │   ├── *scc_sn2.pkl: order param of each lipid's tail sn2
            │   ├── *avg_scc.pkl: average order param of each lipid
            │   └── *std_scc.pkl: standard deviation of order param
            └── ... (other systems)
