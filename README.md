# Cross-Nearest Neighbor Distance Analysis (Multi-Order) with CSR Simulation

Compares cluster-center coordinates (from single-molecule localization microscopy (SMLM) experiments) from two related HDF5 files (dataset A and B)
by computing cross-nearest neighbor (NN) distances: for each point (cluster center) in file A, the distance to the 1st, 2nd, … N-th nearest point (cluster center) in file B. 
A complete spatial randomness (CSR) simulation is run for each file and the same cross-NN analysis is repeated on the simulated coordinates.


## Overview

1. Pair HDF5 files from two input folders (INPUT_DIR_1, INPUT_DIR_2) that share
   the same base name after removing FILE_SUFFIX_1 / FILE_SUFFIX_2.
2. Load cluster centers from the ``locs`` dataset (fields ``x``, ``y``).
3. Apply a region of interest (ROI) per file:
   - Recommended: polygon from a matching Picasso YAML file
     (<base_name> + ROI_YAML_SUFFIX + ".yaml", e.g. ``_ROI_picks.yaml`` with
     key ``Vertices``).
   - Alternatively: global bounding box, mask image, or YAML path via ROI.
4. Measured cross-NN (data):
   - Uses only points (cluster centers) inside the ROI in each file.
   - For each point in A (inside ROI), finds the N nearest neighbors in B (inside ROI).
5. Simulation:
   - Places the same number of points as in the original HDF5, uniformly at random
     inside the YAML polygon ROI (complete spatial randomness, CSR).
   - Saves ``<stem>_simulation.hdf5`` to OUTPUT_FOLDER_1 / OUTPUT_FOLDER_2.
6. Simulated cross-NN: same cross-NN analysis on the simulated coordinates.
7. Distances are converted from pixels to nanometers using PIXEL_SIZE_NM and saved to a CSV file (one input CSV file per input file with one column per NN order (data + simulation))
8. Optional neighbor counting (NEIGHBOR_RADIUS_NM): per point (cluster center), counts
   within the radius for cross (A↔B) and self (A/A, B/B) neighbors, for measured
   (ROI-filtered) and simulated data; saved to a separate CSV with explicit
   column names.

Parts of this script were developed with assistance from OpenAI's ChatGPT and Cursor AI. The resulting code was reviewed, modified, and validated by the author.

## Requirements

* Python **3.11**
* PyYAML **6.0.1**
* pandas **2.1.4**
* h5py **3.9.0**
* numpy **1.26.4**
* scipy **1.11.4**

## Input Data

### HDF5 files

HDF5 files (.hdf5) containing a list of cluster center coordinates in the 'locs' dataset. The script uses the 'x' and 'y' fields of this dataset to calculate nearest neighbor distances.
The input coordinates are assumed to be in pixels. Nearest neighbor distances are converted from pixels to nanometers using PIXEL_SIZE_NM. 
The script was tested on HDF5 files with cluster centers generated from localization data via the clustering algorithm DBSCAN with the 
[Picasso Software](https://github.com/jungmannlab/picasso) version v0.7.3 and with the [PicassoBatchProcess](https://github.com/HeilemannLab/PicassoBatchProcess) Software.  

HDF5 pairing (dataset A vs dataset B):
The script recursively finds HDF5 files ending with FILE_SUFFIX_1 in INPUT_DIR_1 and
FILE_SUFFIX_2 in INPUT_DIR_2. A pair is formed when both files share the same base name
after those suffixes are removed.

For example, with:
FILE_SUFFIX_1 = "_ROI_dbscan_centers.hdf5"
FILE_SUFFIX_2 = "_ROI_dbscan_centers.hdf5"

these two files are paired (base name ``cell1``):
INPUT_DIR_1 / cell1_ROI_dbscan_centers.hdf5
INPUT_DIR_2 / cell1_ROI_dbscan_centers.hdf5

Each base name must be unique within each input directory (after suffix removal).

The script can optionally use a ROI to restrict the analysis to a specific area. If no ROI is specified, all cluster centers in the HDF5 file are used. A global ROI can be defined and applied 
to all input files. Alternatively, individual YAML polygon ROIs can be automatically matched to each HDF5 file.

### YAML files (with ROI vertices)

YAML files (.yaml) contain polygon vertex coordinates under the `Vertices` field. The script was tested on YAML files (polygonal ROIs) 
generated with the [Picasso Software](https://github.com/jungmannlab/picasso) version v0.7.3 (modul "Render": -> Tools -> Tool Settings -> Shape: Polygon, -> Tools -> Pick).

Per-file YAML ROI matching:
If ROI_YAML_SUFFIX is set, each HDF5 file is matched independently to a YAML ROI.
The HDF5 suffix is stripped from the file stem to determine the base name, then the script looks for:

<base_name> + ROI_YAML_SUFFIX + ".yaml"

This is done in both INPUT_DIR_1 and INPUT_DIR_2.

For example, in INPUT_DIR_1 with:
FILE_SUFFIX_1 = "_ROI_dbscan_centers.hdf5"
ROI_YAML_SUFFIX = "_ROI_picks"

these files are associated (base name ``cell1``):
cell1_ROI_dbscan_centers.hdf5
cell1_ROI_picks.yaml

The common base name is:
cell1 

For example:
```
input_data/ 

├── protein1/ 

│ ├── cell1_ROI_picks.yaml 

│ └── cell1_ROI_dbscan_centers.hdf5 

│ 
└── protein2/ 

  ├── cell1_ROI_picks.yaml 
  
  └── cell1_ROI_dbscan_centers.hdf5
```
  

## Configuration

Before running the script, edit the following variables in the **CONFIGURATION** section in `cross_nearest_neighbor_distances.py`:
```python
# Dataset A
INPUT_DIR_1 = r"C:\cross_nearest_neighbor_distances\example_data\input_data\protein1"
FILE_SUFFIX_1 = "_ROI_dbscan_centers.hdf5"
OUTPUT_FOLDER_1 = r"C:\cross_nearest_neighbor_distances\example_data\output_data\protein1"
# Dataset B
INPUT_DIR_2 = r"C:\cross_nearest_neighbor_distances\example_data\input_data\protein2"
FILE_SUFFIX_2 = "_ROI_dbscan_centers.hdf5"
OUTPUT_FOLDER_2 = r"C:\cross_nearest_neighbor_distances\example_data\output_data\protein2"
N_NEAREST_NEIGHBORS = 2  
# Number of nearest neighbors to compute (e.g. 3 → 1st, 2nd, 3rd NN)
NEIGHBOR_RADIUS_NM = 70.0
# If set (e.g. 500.0): per-point neighbor counts within this radius (nm):
# cross (A↔B), self (A/A and B/B, excluding self), measured + simulated.
# None disables this step.
PIXEL_SIZE_NM = 157.0  
# Conversion factor: 1 pixel = X nanometers (used for NN distances)
ROI = None  
# Global ROI definition (used if no per-file YAML is applied)
#
# Options:
#   None → use all data points
#   (x_min, y_min, x_max, y_max) → rectangular ROI
#   "path/to/mask.png" → mask image (255 = valid region)
#   "path/to/file.yaml" → polygon ROI from YAML (same ROI for all files)
ROI_YAML_SUFFIX = "_ROI_picks"  
# Enables per-file ROI from YAML, input e.g. "_ROI_picks"
# For each HDF5 file, the script searches for:
#   <base_name> + ROI_YAML_SUFFIX + ".yaml"
# Supported YAML formats (Picasso / cluster_density):
#   - Polygon ROI (recommended): key "Vertices"
#   - Legacy circle picks: keys "Centers" + "Diameter"
# Example:
#   cell1_ROI_dbscan_centers.hdf5 → cell1_ROI_picks.yaml
# Behavior:
#   If matching YAML is found → overrides ROI
#   If not found:
#       - If ROI is None → use full dataset
#       - Else → fallback to global ROI
# Set to None to disable per-file ROI completely
RANDOM_SEED = None  
# Seed for random number generator (simulation reproducibility)
# None → different random simulation each run
# Integer (e.g. 42) → reproducible simulation results
```

## Installation
```PowerShell
conda create --name cross_nearest_neighbor_distances python=3.11
conda activate cross_nearest_neighbor_distances
cd filepath\cross_nearest_neighbor_distances
conda install --file requirements.txt
```

## Usage

1. Place the YAML and HDF5 input files in the specified folders.
2. Open `cross_nearest_neighbor_distances.py`.
3. Set variables in CONFIGURATION section.
4. Open environment:
```PowerShell
conda activate cross_nearest_neighbor_distances
```
5. Navigate to the file path, where nearest_neighbor_distances.py is stored:
```PowerShell
cd filepath\cross_nearest_neighbor_distances
```
6. Run:
```PowerShell
python cross_nearest_neighbor_distances.py
```

## Output

- CSV with nearest neighbor distances is saved in OUTPUT_FOLDER_1: ``<base_name>_cross_nearest_neighbors.csv``
  Columns: cross_nn_1_nm … cross_nn_N_nm (measured data),
           sim_cross_nn_1_nm … sim_cross_nn_N_nm (simulation).
- Simulation HDF5 files are saved in OUTPUT_FOLDER_1 and OUTPUT_FOLDER_2.
- Optional: neighbor-count CSV is saved in OUTPUT_FOLDER_1: (if NEIGHBOR_RADIUS_NM is set)
  ``<base_name>_neighbor_counts_<radius>nm.csv`` — measured and simulated
  cross- and self-neighbor counts with dataset names in each column header.


