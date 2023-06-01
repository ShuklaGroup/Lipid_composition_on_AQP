Scripts for feature selections (computing residue-residue distances, RRCS calculations, oASIS calculations), Markov State Models building from these features (grid searching, implied timescale plots), and validation of built MSM (reweighting plots and bootstrapping). Dependencies include PyEMMMA 2.5.6. with Python 3.6. Conda environment used is included within [msm-env.yml](msm-env.yml) (NOTE: this environment is REQUIRED for opening and analyzing PyEMMA objects from this project, which are included on [Box](https://uofi.box.com/s/uc33gid1jhyuc0oru8tr30to3x9kyj8z))

- To select features:
  - [*compute-distances.py*](compute-distances.py) is used to obtain residue-residue distances between specified pairs of residues or all residues.
  - Residue-residue contact scoring (RRCS) are calculated for all residue pairs, excluding the N-terminus residues 28-33 (0-5 in Python index) with [*RRCS.py*](RRCS.py)
  - Process RRCS outputs with [*post-RRCS.py*](post-RRCS.py) to get indices of features selected from RRCS
  - Perform oASIS optimization from RRCS outputs with [*do-oASIS.py*](do-oASIS.py), then visualize final features on VMD for rational design
  
- MSM scripts are used for the Markov State Model building workflow:
  - To determine hyperparameters, [*MSM-grid-search.py*](MSM-grid-search.py) is used to perform a grid search. Outputs are .pkl tICA and clustering objects, an implied timescale plot from the resulting hyperparameter sets' MSM, and a dict containing VAMP-2 scores. Five parameter sets with the highest VAMP-2 scores are examined. MSM lagtime is selected based on the implied timescale plot
  - With built MSMs, validations are done with [*MSM-bootstrap.py*](MSM-bootstrap.py) to estimate relative errors through bootstrapping, as well as [*MSM-reweighting-analysis.py*](MSM-reweighting-analysis.py) to evaluate the stationary distribution from the MSM and its deviation from the original distribution

- Final MSM-related objects are included in [Box](https://uofi.box.com/s/uc33gid1jhyuc0oru8tr30to3x9kyj8z) with the following tree

```
.
└── SoPIP2-lipid-GitHub/
    └── 01-features-and-MSM/
        ├── features-per-system
        └── MSM-related-objs/
            ├── POPC/
            │   ├── *features.pkl: indices of protein residue pairs
            │   ├── *tica_obj.pkl: PyEMMA tICA object
            │   ├── *tica_trajs.pkl: PyEMMA tICs for each traj and frame
            │   ├── *cluster_obj.pkl: PyEMMA k-means clustering object
            │   ├── *msm_obj.pkl: PyEMMA MSM object
            │   ├── *ITS-error.pkl: PyEMMA implied timescale object
            │   ├── *weights.pkl: PyEMMA stationary distribution of MSM
            │   └── *probabilities.npy: msm and raw probability of each cluster
            ├── POPE
            ├── POPG
            └── ...
