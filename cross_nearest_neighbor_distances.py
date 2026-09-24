"""
Cross-Nearest Neighbor Distance Analysis (Multi-Order) with CSR Simulation

Compares cluster-center coordinates (from single-molecule localization microscopy (SMLM) experiments) from two related HDF5 files (dataset A and B)
by computing cross-nearest neighbor (NN) distances: for each point (cluster center) in file A, the distance to the 1st, 2nd, … N-th nearest point (cluster center) in file B. 
A complete spatial randomness (CSR) simulation is run for each file and the same cross-NN analysis is repeated on the simulated coordinates.

Author: Tanja Menche
Affiliation: Research group of Mike Heilemann, Goethe University Frankfurt am Main, Germany
Version: v1.0.0
Date: 2026-09-22

About this script:
------------------
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

User Inputs (edit the variables below):
---------------------------------------
- INPUT_DIR_1: The root directory for dataset A containing the HDF5 files and optionally the YAML files with the polygonal ROI coordinates.
- FILE_SUFFIX_1: The suffix that the HDF5 files must end with for dataset A. For example, '_protein1_ROI_dbscan_centers.hdf5'. Use '.hdf5' to process all HDF5 files.
- OUTPUT_FOLDER_1: The folder where the output CSV files and simulated HDF5 files are saved for dataset A.
- INPUT_DIR_2: The root directory for dataset B containing the HDF5 files and optionally the YAML files with the polygonal ROI coordinates.
- FILE_SUFFIX_2: The suffix that the HDF5 files must end with for dataset B. For example, '_protein1_ROI_dbscan_centers.hdf5'. Use '.hdf5' to process all HDF5 files.
- OUTPUT_FOLDER_2: The folder where the output CSV files and simulated HDF5 files are saved for dataset B.
- N_NEAREST_NEIGHBORS: The number of nearest neighbors to calculate. For example, 4 calculates the 1st, 2nd, 3rd, and 4th nearest neighbors.
- NEIGHBOR_RADIUS_NM: per-point counts of neighboring cluster centers within this radius in nm: cross (A↔B), self (A/A and B/B, excluding self), measured + simulated. None disables this step.
- PIXEL_SIZE_NM: The pixel size in nanometers. This value is used to convert the calculated nearest neighbor distances from pixels to nanometers.
- ROI: The global Region of Interest to use for filtering the data. Set to None to use all data points, provide a bounding box as (x_min, y_min, x_max, y_max), or provide a path to a YAML polygon or image mask.
- ROI_YAML_SUFFIX: The suffix used to identify per-file YAML ROI polygons. For example, '_ROI_picks'. Set to None to disable automatic per-file YAML ROI matching.
- RANDOM_SEED: The random seed used for the CSR simulation. Set to None to generate different random simulations on each run. Specify an integer, such as 42, to make the simulation reproducible.

Information about the input data:
--------------------------------
The HDF5 files contain cluster center coordinates in the 'locs' dataset. The script uses the 'x' and 'y' fields of this dataset to calculate nearest neighbor distances.
The input coordinates are assumed to be in pixels. Nearest neighbor distances are converted from pixels to nanometers using PIXEL_SIZE_NM.
The script can optionally use a ROI to restrict the analysis to a specific area. If no ROI is specified, all cluster centers in the HDF5 file are used. A global ROI can be defined and applied  to all input files (bounding box, mask image, or a single YAML path). Alternatively, a YAML ROI can be matched automatically to each HDF5 file (polygon Vertices or Picasso circle picks).
Corresponding HDF5 files and YAML files can be generated with the Picasso Software
(https://github.com/jungmannlab/picasso) version v0.7.3 from SMLM data (see README.md).

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

The YAML file must contain a 'Vertices' entry defining the polygon coordinates of the ROI. The script uses these vertices to determine which cluster centers are located within the polygon and to simulate a CSR distribution of points within the selected ROI.

How to use:
-----------
1. Edit the variables in the CONFIGURATION section below to match your data and analysis requirements. Optionally, set RANDOM_SEED to an integer if reproducible simulations are required.
2. Run the script.
   The script will recursively process all matching HDF5 files.
3. For each processed file, the script saves a CSV file containing the experimental and simulated nearest neighbor distances. The script also saves the simulated cluster centers as an HDF5 file 
accompanied by a YAML file of the same filename for compatibility with the Picasso Software.

Outputs (per matched file pair)
-------------------------------
- CSV in OUTPUT_FOLDER_1: ``<base_name>_cross_nearest_neighbors.csv``
  Columns: cross_nn_1_nm … cross_nn_N_nm (measured data),
           sim_cross_nn_1_nm … sim_cross_nn_N_nm (simulation).
- Simulation HDF5 files in OUTPUT_FOLDER_1 and OUTPUT_FOLDER_2.
- Optional neighbor-count CSV (if NEIGHBOR_RADIUS_NM is set):
  ``<base_name>_neighbor_counts_<radius>nm.csv`` — measured and simulated
  cross- and self-neighbor counts with dataset names in each column header.
  
"""

import os
import argparse
from pathlib import Path
from typing import Optional, Tuple, Union, List
import h5py
import numpy as np
import pandas as pd
import yaml
from scipy import spatial

# =============================================================================
# CONFIGURATION
# =============================================================================

# --- Input / Output -----------------------------------------------------------

# Dataset A
INPUT_DIR_1 = r"C:\cross_nearest_neighbor_distances\example_data\input_data\protein1"
FILE_SUFFIX_1 = "_ROI_dbscan_centers.hdf5"
OUTPUT_FOLDER_1 = r"C:\cross_nearest_neighbor_distances\example_data\output_data\protein1"

# Dataset B
INPUT_DIR_2 = r"C:\cross_nearest_neighbor_distances\example_data\input_data\protein2"
FILE_SUFFIX_2 = "_ROI_dbscan_centers.hdf5"
OUTPUT_FOLDER_2 = r"C:\cross_nearest_neighbor_distances\example_data\output_data\protein2"


# --- Nearest Neighbor Analysis -----------------------------------------------

N_NEAREST_NEIGHBORS = 2  
# Number of nearest neighbors to compute (e.g. 3 → 1st, 2nd, 3rd NN)

NEIGHBOR_RADIUS_NM = 70.0
# If set (e.g. 500.0): per-point neighbor counts within this radius (nm):
# cross (A↔B), self (A/A and B/B, excluding self), measured + simulated.
# None disables this step.


# --- Units -------------------------------------------------------------------

PIXEL_SIZE_NM = 157.0  
# Conversion factor: 1 pixel = X nanometers (used for NN distances)


# --- ROI (Region of Interest) ------------------------------------------------

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
#
# For each HDF5 file, the script searches for:
#   <base_name> + ROI_YAML_SUFFIX + ".yaml"
#
# Supported YAML formats (Picasso / cluster_density):
#   - Polygon ROI (recommended): key "Vertices"
#   - Legacy circle picks: keys "Centers" + "Diameter"
#
# Example:
#   cell1_ROI_dbscan_centers.hdf5 → cell1_ROI_picks.yaml
#
# Behavior:
#   If matching YAML is found → overrides ROI
#   If not found:
#       - If ROI is None → use full dataset
#       - Else → fallback to global ROI
#
# Set to None to disable per-file ROI completely


# --- Simulation --------------------------------------------------------------

RANDOM_SEED = None  
# Seed for random number generator (simulation reproducibility)
#
# None → different random simulation each run
# Integer (e.g. 42) → reproducible simulation results

# =============================================================================



ROI_COORDS_IN_NM = False  
# If True: YAML ROI coordinates are scaled to nm using PIXEL_SIZE_NM
# If False: YAML ROI coordinates are assumed to already match data units (usually pixels)

ROI_PICKS_CONVEX_HULL_RATIO = None
# For Picasso Centers picks: if Diameter < ratio * pick spread, use convex hull
# instead of tiny circles (Diameter 1 px with spread ~150 px → hull).

# Optional: use utils if available in same directory
try:
    from utils import load_hdf5, save_hdf5, load_info, load_mask
except ImportError:
    # Fallback implementations
    def _load_info_fallback(file_path):
        path_base = str(file_path).rsplit(".", 1)[0]
        filename = path_base + ".yaml"
        try:
            with open(filename, "r") as info_file:
                return list(yaml.load_all(info_file, Loader=yaml.UnsafeLoader))
        except FileNotFoundError:
            return []

    def load_hdf5(file_path):
        with h5py.File(file_path, 'r') as hdf5_file:
            locs_data = hdf5_file['locs'][:]
            info = _load_info_fallback(file_path)
            return locs_data, info

    def save_hdf5(file_path, locs, info):
        with h5py.File(file_path, "w") as locs_file:
            locs_file.create_dataset("locs", data=locs)
        if info:
            yaml_path = str(file_path).replace('.hdf5', '.yaml')
            with open(yaml_path, "w") as f:
                yaml.dump_all(info, f, default_flow_style=False)

    def load_info(file_path):
        return _load_info_fallback(file_path)

    def load_mask(mask_path):
        try:
            import cv2
            mask = cv2.imread(str(mask_path), cv2.IMREAD_GRAYSCALE)
            return mask
        except (ImportError, Exception):
            return None


# ROI type: None, bbox tuple, mask path, YAML path, or dict with 'polygon' key
ROI_TYPE = Union[
    None,
    Tuple[float, float, float, float],
    str,
    Path,
    dict,  # polygon: {'polygon': [...], 'bounds': ...}
           # circles:  {'circles': [(cx, cy, r), ...], 'bounds': ...}
]


def _point_in_polygon(x: np.ndarray, y: np.ndarray, vertices: List[Tuple[float, float]]) -> np.ndarray:
    """Point-in-polygon test. Returns boolean array."""
    try:
        from matplotlib.path import Path
        poly_path = Path(vertices)
        points = np.column_stack([x, y])
        return poly_path.contains_points(points)
    except ImportError:
        # Fallback: ray-casting
        n = len(vertices)
        if n < 3:
            return np.zeros(len(x), dtype=bool)
        inside = np.zeros(len(x), dtype=bool)
        for i in range(n):
            j = (i + 1) % n
            xi, yi = vertices[i]
            xj, yj = vertices[j]
            if abs(yj - yi) < 1e-10:
                continue
            mask = ((yi > y) != (yj > y)) & (x < (xj - xi) * (y - yi) / (yj - yi) + xi)
            inside ^= mask
        return inside


def _point_in_circles(
    x: np.ndarray,
    y: np.ndarray,
    circles: List[Tuple[float, float, float]],
) -> np.ndarray:
    """True where point lies inside any circle (cx, cy, radius)."""
    if not circles:
        return np.zeros(len(x), dtype=bool)
    inside = np.zeros(len(x), dtype=bool)
    for cx, cy, radius in circles:
        inside |= (x - cx) ** 2 + (y - cy) ** 2 <= radius ** 2
    return inside


def _points_in_roi_dict(x: np.ndarray, y: np.ndarray, roi: dict) -> np.ndarray:
    """Point-in-ROI test for YAML-derived polygon or circle-pick ROIs."""
    if 'polygon' in roi:
        return _point_in_polygon(x, y, roi['polygon'])
    if 'circles' in roi:
        return _point_in_circles(x, y, roi['circles'])
    return np.zeros(len(x), dtype=bool)


def _read_roi_from_yaml(
    yaml_path: Path,
    pixel_size_nm: float = 1.0,
    scale_to_nm: bool = False
) -> Optional[dict]:
    """
    Read ROI from YAML.

    Supports:
    - Polygon picks (cluster_density / Picasso): key ``Vertices``
    - Circular picks (Picasso): keys ``Centers`` + ``Diameter`` (pixels) or ``Diameter (nm)``

    scale_to_nm: if True, multiply coordinates by pixel_size_nm (for locs in nm).
    """
    try:
        with open(yaml_path, 'r') as f:
            data = yaml.safe_load(f)
        if not data:
            return None

        scale = pixel_size_nm if scale_to_nm else 1.0

        if 'Vertices' in data:
            vertices = [
                (coord[0] * scale, coord[1] * scale)
                for coord_group in data['Vertices']
                for coord in coord_group
            ]
            if len(vertices) < 3:
                return None
            xs = [v[0] for v in vertices]
            ys = [v[1] for v in vertices]
            bounds = (min(xs), min(ys), max(xs), max(ys))
            return {'polygon': vertices, 'bounds': bounds, 'roi_kind': 'polygon'}

        if 'Centers' in data:
            if 'Diameter (nm)' in data:
                diameter_px = float(data['Diameter (nm)']) / pixel_size_nm
            else:
                diameter_px = float(data.get('Diameter', 1.0))

            pts = np.array(
                [[float(c[0]) * scale, float(c[1]) * scale] for c in data['Centers']],
                dtype=np.float64,
            )
            if len(pts) == 0:
                return None

            spread = max(float(np.ptp(pts[:, 0])), float(np.ptp(pts[:, 1])), 1e-9)
            use_hull = (
                len(pts) >= 3
                and diameter_px < ROI_PICKS_CONVEX_HULL_RATIO * spread
            )

            if use_hull:
                try:
                    hull = spatial.ConvexHull(pts)
                    vertices = [tuple(pts[i]) for i in hull.vertices]
                    xs = pts[:, 0]
                    ys = pts[:, 1]
                    bounds = (float(xs.min()), float(ys.min()), float(xs.max()), float(ys.max()))
                    return {'polygon': vertices, 'bounds': bounds, 'roi_kind': 'convex_hull'}
                except Exception:
                    pass

            radius_px = diameter_px / 2.0
            circles = [(pts[i, 0], pts[i, 1], radius_px) for i in range(len(pts))]
            x_mins = [c[0] - c[2] for c in circles]
            y_mins = [c[1] - c[2] for c in circles]
            x_maxs = [c[0] + c[2] for c in circles]
            y_maxs = [c[1] + c[2] for c in circles]
            bounds = (min(x_mins), min(y_mins), max(x_maxs), max(y_maxs))
            return {'circles': circles, 'bounds': bounds, 'roi_kind': 'circles'}

        return None
    except Exception:
        return None


def _read_roi_polygon_from_yaml(
    yaml_path: Path,
    pixel_size_nm: float = 1.0,
    scale_to_nm: bool = False
) -> Optional[dict]:
    """Backward-compatible alias for :func:`_read_roi_from_yaml`."""
    return _read_roi_from_yaml(yaml_path, pixel_size_nm, scale_to_nm)


def _get_roi_for_hdf5(
    hdf5_path: Path,
    input_dir: Path,
    roi_yaml_suffix: Optional[str],
    hdf5_file_suffix: Optional[str] = None
) -> Optional[dict]:
    """Find matching YAML ROI file for an HDF5 (cluster_density_v1 style)."""
    if not roi_yaml_suffix:
        return None
    stem = hdf5_path.stem
    # Remove HDF5 suffix to get base name (like cluster_density get_base_name)
    if hdf5_file_suffix:
        suffix_no_ext = hdf5_file_suffix.replace(".hdf5", "").replace(".h5", "")
        if stem.endswith(suffix_no_ext):
            stem = stem[: -len(suffix_no_ext)]
    else:
        for suffix in ['_ROI_fl_dbscan_centers_filter', '_filtered_in', '_centers', '_dbscan']:
            if stem.endswith(suffix):
                stem = stem[: -len(suffix)]
                break
    yaml_name = stem + roi_yaml_suffix + ".yaml"
    # Look in same directory as HDF5 first
    yaml_path = hdf5_path.parent / yaml_name
    if yaml_path.exists():
        return _read_roi_from_yaml(yaml_path, PIXEL_SIZE_NM, scale_to_nm=ROI_COORDS_IN_NM)
    # Search recursively in input_dir
    for p in Path(input_dir).rglob(yaml_name):
        if p.exists():
            return _read_roi_from_yaml(p, PIXEL_SIZE_NM, scale_to_nm=ROI_COORDS_IN_NM)
    return None


def filter_locs_by_roi(
    locs: np.ndarray,
    roi: Optional[ROI_TYPE]
) -> np.ndarray:
    """
    Filter localization data to points within ROI.

    Parameters
    ----------
    locs : np.ndarray
        Structured array with 'x' and 'y' fields
    roi : None, tuple, dict, or str/Path
        - None: use full data (no filtering)
        - (x_min, y_min, x_max, y_max): bounding box
        - dict with 'polygon': YAML polygon ROI
        - str/Path to .yaml: polygon from YAML
        - str/Path to image: mask (255 = valid)

    Returns
    -------
    np.ndarray : Filtered locs
    """
    x = locs['x']
    y = locs['y']

    if roi is None:
        return locs

    if isinstance(roi, dict) and ('polygon' in roi or 'circles' in roi):
        valid = _points_in_roi_dict(x, y, roi)
        return locs[valid]

    if isinstance(roi, (tuple, list)) and len(roi) == 4:
        x_min, y_min, x_max, y_max = roi
        valid = (x >= x_min) & (x <= x_max) & (y >= y_min) & (y <= y_max)
        return locs[valid]

    # Path to YAML or mask
    roi_path = Path(roi) if isinstance(roi, (str, Path)) else None
    if roi_path and roi_path.suffix.lower() in ('.yaml', '.yml'):
        poly_dict = _read_roi_from_yaml(roi_path, PIXEL_SIZE_NM, scale_to_nm=ROI_COORDS_IN_NM)
        if poly_dict:
            valid = _point_in_polygon(x, y, poly_dict['polygon'])
            return locs[valid]
        print(f"Warning: Could not read ROI from YAML {roi}, using full data")
        return locs

    # ROI is mask path
    mask = load_mask(roi)
    if mask is None:
        print(f"Warning: Could not load ROI mask from {roi}, using full data")
        return locs

    x_idx = np.clip(x.astype(int), 0, mask.shape[1] - 1)
    y_idx = np.clip(y.astype(int), 0, mask.shape[0] - 1)
    valid = mask[y_idx, x_idx] == 255
    return locs[valid]


def get_roi_bounds_for_simulation(
    locs: np.ndarray,
    roi: Optional[ROI_TYPE]
) -> Tuple[float, float, float, float]:
    """Get (x_min, y_min, x_max, y_max) for simulation area."""
    if roi is None:
        return (
            float(locs['x'].min()),
            float(locs['y'].min()),
            float(locs['x'].max()),
            float(locs['y'].max())
        )
    if isinstance(roi, dict) and 'bounds' in roi:
        return roi['bounds']
    if isinstance(roi, (tuple, list)) and len(roi) == 4:
        return tuple(roi)
    # For mask: use mask bounding box of valid pixels
    if isinstance(roi, (str, Path)):
        mask = load_mask(roi)
        if mask is not None:
            coords = np.where(mask == 255)
            if len(coords[0]) > 0:
                y_min, y_max = coords[0].min(), coords[0].max()
                x_min, x_max = coords[1].min(), coords[1].max()
                return (float(x_min), float(y_min), float(x_max), float(y_max))
    return (
        float(locs['x'].min()),
        float(locs['y'].min()),
        float(locs['x'].max()),
        float(locs['y'].max())
    )


def get_roi_area_for_simulation(
    roi: Optional[Union[Tuple[float, float, float, float], str, Path]],
    bounds: Tuple[float, float, float, float]
) -> float:
    """Get area of ROI for density calculation."""
    if isinstance(roi, (tuple, list)) and len(roi) == 4:
        x_min, y_min, x_max, y_max = roi
        return (x_max - x_min) * (y_max - y_min)
    if roi is None:
        x_min, y_min, x_max, y_max = bounds
        return (x_max - x_min) * (y_max - y_min)
    # Mask: count valid pixels
    mask = load_mask(roi)
    if mask is None:
        x_min, y_min, x_max, y_max = bounds
        return (x_max - x_min) * (y_max - y_min)
    return float(np.sum(mask == 255))


def cross_nearest_neighbor_distances(
    coords_source: np.ndarray,
    coords_target: np.ndarray,
    k_neighbors: int
) -> np.ndarray:
    """
    Compute nearest-neighbor distances from coords_source
    to coords_target.

    Example:
        source = file A
        target = file B

    Returns:
        distances from every point in source
        to nearest neighbors in target.
    """

    if len(coords_source) == 0 or len(coords_target) == 0:
        return np.empty((0, k_neighbors), dtype=np.float64)

    tree = spatial.cKDTree(coords_target.astype(np.float64))

    k_use = min(k_neighbors, len(coords_target))

    d, _ = tree.query(
        coords_source.astype(np.float64),
        k=k_use
    )

    if k_use < k_neighbors:
        pad = np.full(
            (len(coords_source), k_neighbors - k_use),
            np.nan
        )
        d = np.hstack([d, pad])

    if k_neighbors == 1:
        d = d[:, np.newaxis]

    return d


def count_neighbors_within_radius(
    coords_source: np.ndarray,
    coords_target: np.ndarray,
    radius_px: float,
    exclude_self: bool = False,
) -> np.ndarray:
    """
    For each point in coords_source, count points in coords_target within radius_px.

    Coordinates must be in the same units as radius_px (pixels).
    If exclude_self is True and source is the same as target, the query point itself
    is not counted.
    """
    if len(coords_source) == 0:
        return np.array([], dtype=np.int64)
    if len(coords_target) == 0:
        return np.zeros(len(coords_source), dtype=np.int64)

    tree = spatial.cKDTree(coords_target.astype(np.float64))
    neighbors = tree.query_ball_point(
        coords_source.astype(np.float64),
        r=float(radius_px),
    )
    counts = np.array([len(nbrs) for nbrs in neighbors], dtype=np.int64)
    if exclude_self and len(coords_source) == len(coords_target):
        counts = np.maximum(counts - 1, 0)
    return counts


def compute_neighbor_count_sets(
    coords_a: np.ndarray,
    coords_b: np.ndarray,
    radius_px: float,
) -> dict:
    """Cross- and self-neighbor counts for dataset A and B."""
    return {
        "a_cross_b": count_neighbors_within_radius(coords_a, coords_b, radius_px),
        "b_cross_a": count_neighbors_within_radius(coords_b, coords_a, radius_px),
        "a_self": count_neighbors_within_radius(
            coords_a, coords_a, radius_px, exclude_self=True
        ),
        "b_self": count_neighbors_within_radius(
            coords_b, coords_b, radius_px, exclude_self=True
        ),
    }


def build_neighbor_counts_dataframe(
    label_a: str,
    label_b: str,
    radius_nm: float,
    data_counts: dict,
    sim_counts: dict,
) -> pd.DataFrame:
    """
    Build a wide CSV table for measured and simulated cross- and self-neighbor
    counts. Each count column reads as: how many points from dataset X are
    within the radius of this point in dataset Y (e.g.
    measured_n_neighbors_B_close_to_this_point_A = B points near one A point).
    """
    r = f"{radius_nm:g}"

    n_a = max(
        len(data_counts["a_cross_b"]),
        len(data_counts["a_self"]),
        len(sim_counts["a_cross_b"]),
        len(sim_counts["a_self"]),
        0,
    )
    n_b = max(
        len(data_counts["b_cross_a"]),
        len(data_counts["b_self"]),
        len(sim_counts["b_cross_a"]),
        len(sim_counts["b_self"]),
        0,
    )
    n_rows = max(n_a, n_b)

    return pd.DataFrame({
        "dataset_A_file": [label_a] * n_rows,
        "dataset_B_file": [label_b] * n_rows,
        "point_index_A": _pad_1d(np.arange(n_a, dtype=np.float64), n_rows),
        f"measured_n_neighbors_B_close_to_this_point_A_{r}nm": _pad_1d(
            data_counts["a_cross_b"].astype(np.float64), n_rows
        ),
        f"measured_n_neighbors_A_close_to_this_point_A_{r}nm_excl_self": _pad_1d(
            data_counts["a_self"].astype(np.float64), n_rows
        ),
        f"simulated_n_neighbors_B_close_to_this_point_A_{r}nm": _pad_1d(
            sim_counts["a_cross_b"].astype(np.float64), n_rows
        ),
        f"simulated_n_neighbors_A_close_to_this_point_A_{r}nm_excl_self": _pad_1d(
            sim_counts["a_self"].astype(np.float64), n_rows
        ),
        "point_index_B": _pad_1d(np.arange(n_b, dtype=np.float64), n_rows),
        f"measured_n_neighbors_A_close_to_this_point_B_{r}nm": _pad_1d(
            data_counts["b_cross_a"].astype(np.float64), n_rows
        ),
        f"measured_n_neighbors_B_close_to_this_point_B_{r}nm_excl_self": _pad_1d(
            data_counts["b_self"].astype(np.float64), n_rows
        ),
        f"simulated_n_neighbors_A_close_to_this_point_B_{r}nm": _pad_1d(
            sim_counts["b_cross_a"].astype(np.float64), n_rows
        ),
        f"simulated_n_neighbors_B_close_to_this_point_B_{r}nm_excl_self": _pad_1d(
            sim_counts["b_self"].astype(np.float64), n_rows
        ),
    })


def simulate_clusters_same_density(
    locs: np.ndarray,
    roi: Optional[ROI_TYPE],
    random_seed: Optional[int] = None,
    n_points: Optional[int] = None,
) -> np.ndarray:
    """
    Simulate cluster centers with Complete Spatial Randomness (CSR) inside ROI.

    Parameters
    ----------
    locs : np.ndarray
        Template locs array (dtype / metadata); may be filtered or full.
    roi : ROI specification (None, bbox, polygon dict, mask path)
    random_seed : int, optional
    n_points : int, optional
        Number of points to simulate. Default: len(locs). Use the original
        HDF5 count when analysis uses ROI-filtered data but simulation should not.

    Returns
    -------
    np.ndarray : Simulated locs with same structure, new x,y coordinates
    """
    if random_seed is not None:
        rng = np.random.default_rng(random_seed)
    else:
        rng = np.random.default_rng()

    n_points = len(locs) if n_points is None else int(n_points)
    if n_points == 0:
        return np.array([], dtype=locs.dtype)

    bounds = get_roi_bounds_for_simulation(locs, roi)
    x_min, y_min, x_max, y_max = bounds

    # YAML polygon or Picasso circle picks: rejection sampling inside ROI
    if isinstance(roi, dict) and ('polygon' in roi or 'circles' in roi):
        x_sim_list, y_sim_list = [], []
        for _ in range(500):  # Safety limit
            if len(x_sim_list) >= n_points:
                break
            batch = max((n_points - len(x_sim_list)) * 4, 100)
            x_cand = rng.uniform(x_min, x_max, batch)
            y_cand = rng.uniform(y_min, y_max, batch)
            inside = _points_in_roi_dict(x_cand, y_cand, roi)
            x_sim_list.extend(x_cand[inside].tolist())
            y_sim_list.extend(y_cand[inside].tolist())
        n_got = len(x_sim_list)
        x_sim = np.empty(n_points, dtype=np.float64)
        y_sim = np.empty(n_points, dtype=np.float64)
        if n_got > 0:
            x_sim[:n_got] = x_sim_list[:n_points]
            y_sim[:n_got] = y_sim_list[:n_points]
        if n_got < n_points:
            print(
                f"  Warning: placed {n_got}/{n_points} simulated points inside YAML ROI; "
                "check Diameter / coordinate units."
            )
            if n_got == 0:
                x_sim[:] = rng.uniform(x_min, x_max, n_points)
                y_sim[:] = rng.uniform(y_min, y_max, n_points)
            else:
                x_sim[n_got:] = rng.uniform(x_min, x_max, n_points - n_got)
                y_sim[n_got:] = rng.uniform(y_min, y_max, n_points - n_got)
    # Rectangular ROI
    elif roi is None or (isinstance(roi, (tuple, list)) and len(roi) == 4):
        x_sim = rng.uniform(x_min, x_max, n_points)
        y_sim = rng.uniform(y_min, y_max, n_points)
    # Mask: sample within valid pixels
    elif isinstance(roi, (str, Path)):
        mask = load_mask(roi)
        if mask is None:
            x_sim = rng.uniform(x_min, x_max, n_points)
            y_sim = rng.uniform(y_min, y_max, n_points)
        else:
            valid_indices = np.argwhere(mask == 255)
            if len(valid_indices) == 0:
                x_sim = rng.uniform(x_min, x_max, n_points)
                y_sim = rng.uniform(y_min, y_max, n_points)
            else:
                chosen = rng.choice(len(valid_indices), size=n_points, replace=True)
                selected = valid_indices[chosen]
                jitter = rng.random(size=(n_points, 2))
                coords = selected + jitter
                y_sim = coords[:, 0]
                x_sim = coords[:, 1]
    else:
        x_sim = rng.uniform(x_min, x_max, n_points)
        y_sim = rng.uniform(y_min, y_max, n_points)

    sim_locs = np.empty(n_points, dtype=locs.dtype)
    if len(locs) > 0:
        template = locs[0]
        for name in locs.dtype.names:
            if name not in ('x', 'y'):
                sim_locs[name] = template[name]
    sim_locs['x'] = x_sim
    sim_locs['y'] = y_sim
    return sim_locs


def process_file_pair(
    file_1: Path,
    file_2: Path,
    roi_1,
    roi_2,
    random_seed: Optional[int],
    output_dir_1: Optional[Path] = None,
    output_dir_2: Optional[Path] = None,
    neighbor_radius_nm: Optional[float] = None,
):
    """
    Cross-nearest-neighbor analysis between two files.

    Computes:
        file1 -> file2 NN distances

    Simulations:
        simulated_file1 -> simulated_file2 NN distances
    """

    # Load data
    locs_1, info_1 = load_hdf5(str(file_1))
    locs_2, info_2 = load_hdf5(str(file_2))

    n_orig_1, n_orig_2 = len(locs_1), len(locs_2)

    # ROI filtering (for measured cross-NN analysis only)
    locs_1_roi = filter_locs_by_roi(locs_1, roi_1)
    locs_2_roi = filter_locs_by_roi(locs_2, roi_2)

    print(
        f"  Dataset A: {n_orig_1} points in HDF5, {len(locs_1_roi)} inside ROI"
    )
    print(
        f"  Dataset B: {n_orig_2} points in HDF5, {len(locs_2_roi)} inside ROI"
    )

    if len(locs_1_roi) == 0 or len(locs_2_roi) == 0:
        print("  Warning: One dataset empty after ROI filtering")
        return np.array([]), np.array([]), None

    coords_1 = np.column_stack([locs_1_roi['x'], locs_1_roi['y']])
    coords_2 = np.column_stack([locs_2_roi['x'], locs_2_roi['y']])

    neighbor_counts = None

    # CROSS NN ANALYSIS
    nn_distances = cross_nearest_neighbor_distances(
        coords_1,
        coords_2,
        N_NEAREST_NEIGHBORS
    ) * PIXEL_SIZE_NM

    # SIMULATIONS: same count as original HDF5, uniform random inside YAML ROI
    sim_locs_1 = simulate_clusters_same_density(
        locs_1,
        roi_1,
        random_seed,
        n_points=n_orig_1,
    )

    sim_locs_2 = simulate_clusters_same_density(
        locs_2,
        roi_2,
        random_seed,
        n_points=n_orig_2,
    )

    print(
        f"  Simulated {len(sim_locs_1)} + {len(sim_locs_2)} points "
        f"(matched to original HDF5 counts)"
    )

    coords_sim_1 = np.column_stack([
        sim_locs_1['x'],
        sim_locs_1['y']
    ])

    coords_sim_2 = np.column_stack([
        sim_locs_2['x'],
        sim_locs_2['y']
    ])

    nn_distances_sim = cross_nearest_neighbor_distances(
        coords_sim_1,
        coords_sim_2,
        N_NEAREST_NEIGHBORS
    ) * PIXEL_SIZE_NM

    if neighbor_radius_nm is not None and neighbor_radius_nm > 0:
        radius_px = neighbor_radius_nm / PIXEL_SIZE_NM
        data_counts = compute_neighbor_count_sets(coords_1, coords_2, radius_px)
        sim_counts = compute_neighbor_count_sets(coords_sim_1, coords_sim_2, radius_px)
        neighbor_counts = {
            "label_a": file_1.stem,
            "label_b": file_2.stem,
            "radius_nm": neighbor_radius_nm,
            "data": data_counts,
            "sim": sim_counts,
        }
        print(
            f"  Neighbor counts within {neighbor_radius_nm:g} nm "
            f"({radius_px:.4f} px): cross A↔B, self A/A, self B/B "
            f"(measured + simulated)"
        )

    # Save simulations
    if output_dir_1 is not None:
        output_dir_1.mkdir(parents=True, exist_ok=True)
        sim_path_1 = output_dir_1 / f"{file_1.stem}_simulation.hdf5"
        save_hdf5(str(sim_path_1), sim_locs_1, info_1)

    if output_dir_2 is not None:
        output_dir_2.mkdir(parents=True, exist_ok=True)
        sim_path_2 = output_dir_2 / f"{file_2.stem}_simulation.hdf5"
        save_hdf5(str(sim_path_2), sim_locs_2, info_2)

    return nn_distances, nn_distances_sim, neighbor_counts


def find_hdf5_files(root_dir: Path, file_suffix: str) -> List[Path]:
    """Recursively find all HDF5 files matching the suffix in subfolders."""
    root_dir = Path(root_dir)
    if not root_dir.exists():
        return []

    files = []
    for path in root_dir.rglob("*"):
        if path.is_file() and path.suffix.lower() in ('.hdf5', '.h5'):
            if file_suffix and not path.name.endswith(file_suffix):
                continue
            if "_simulation" in path.name:
                continue  # Skip our own output files
            files.append(path)
    return sorted(files)

def pair_hdf5_files(
    files_1: List[Path],
    files_2: List[Path],
    suffix_1: str,
    suffix_2: str
):
    """
    Pair files from dataset 1 and dataset 2 using matching base names.
    """

    def remove_suffix(name: str, suffix: str):
        if name.endswith(suffix):
            return name[:-len(suffix)]
        return Path(name).stem

    dict_1 = {
        remove_suffix(f.name, suffix_1): f
        for f in files_1
    }

    dict_2 = {
        remove_suffix(f.name, suffix_2): f
        for f in files_2
    }

    common_keys = sorted(set(dict_1.keys()) & set(dict_2.keys()))

    pairs = []

    for key in common_keys:
        pairs.append((key, dict_1[key], dict_2[key]))

    return pairs

def _pad_1d(values: np.ndarray, length: int, fill: float = np.nan) -> np.ndarray:
    """Pad or truncate a 1D array to a fixed length (for aligned CSV columns)."""
    out = np.full(length, fill, dtype=np.float64)
    n = min(len(values), length)
    if n > 0:
        out[:n] = np.asarray(values[:n], dtype=np.float64)
    return out


def _sanitize_column_name(name: str) -> str:
    """Create a valid CSV column name from file path."""
    s = str(name).replace("\\", "_").replace("/", "_").replace(".", "_")
    # Remove or replace characters that might cause issues
    s = "".join(c if c.isalnum() or c == "_" else "_" for c in s)
    return s[:80]  # Limit length


def run_batch_analysis(
    input_dir_1,
    file_suffix_1,
    input_dir_2,
    file_suffix_2,
    roi=None,
    roi_yaml_suffix=None,
    output_folder_1=None,
    output_folder_2=None,
    random_seed=None,
    neighbor_radius_nm: Optional[float] = None,
):

    input_dir_1 = Path(input_dir_1)
    input_dir_2 = Path(input_dir_2)
    
    output_dir_1 = Path(output_folder_1)
    output_dir_2 = Path(output_folder_2)

    output_dir_1.mkdir(parents=True, exist_ok=True)
    output_dir_2.mkdir(parents=True, exist_ok=True)

    files_1 = find_hdf5_files(input_dir_1, file_suffix_1)
    files_2 = find_hdf5_files(input_dir_2, file_suffix_2)

    print(f"Found {len(files_1)} file(s) in dataset 1 (suffix: {file_suffix_1!r})")
    print(f"Found {len(files_2)} file(s) in dataset 2 (suffix: {file_suffix_2!r})")

    pairs = pair_hdf5_files(
        files_1,
        files_2,
        file_suffix_1,
        file_suffix_2
    )

    print(f"Found {len(pairs)} matching file pairs")

    if len(pairs) == 0 and (files_1 or files_2):

        def _base_key(path: Path, suffix: str) -> str:
            name = path.name
            if suffix and name.endswith(suffix):
                return name[: -len(suffix)]
            return path.stem

        keys_1 = {_base_key(f, file_suffix_1) for f in files_1}
        keys_2 = {_base_key(f, file_suffix_2) for f in files_2}
        only_1 = sorted(keys_1 - keys_2)
        only_2 = sorted(keys_2 - keys_1)
        if only_1:
            print("  Base names only in dataset 1:", ", ".join(only_1))
        if only_2:
            print("  Base names only in dataset 2:", ", ".join(only_2))
    elif len(pairs) == 0:
        print("  Hint: check INPUT_DIR paths and FILE_SUFFIX values match your filenames.")
    
    for base_name, fp1, fp2 in pairs:

        print(f"\nProcessing pair:")
        print(f"  A: {fp1.name}")
        print(f"  B: {fp2.name}")

        # ROI handling
        roi_1 = roi
        roi_2 = roi

        if roi_yaml_suffix:

            poly_roi_1 = _get_roi_for_hdf5(
                fp1,
                input_dir_1,
                roi_yaml_suffix,
                hdf5_file_suffix=file_suffix_1
            )

            poly_roi_2 = _get_roi_for_hdf5(
                fp2,
                input_dir_2,
                roi_yaml_suffix,
                hdf5_file_suffix=file_suffix_2
            )
            
            if poly_roi_1 is not None:
                roi_1 = poly_roi_1
                roi_type = poly_roi_1.get("roi_kind") or (
                    "polygon" if "polygon" in poly_roi_1 else "circle picks"
                )
                print(f"  ROI A: YAML ({roi_type}) from {roi_yaml_suffix}")
            else:
                expected = fp1.stem
                if file_suffix_1:
                    suffix_no_ext = file_suffix_1.replace(".hdf5", "").replace(".h5", "")
                    if expected.endswith(suffix_no_ext):
                        expected = expected[: -len(suffix_no_ext)]
                print(
                    f"  ROI A: no YAML found ({expected}{roi_yaml_suffix}.yaml); "
                    "using data bounding box"
                )

            if poly_roi_2 is not None:
                roi_2 = poly_roi_2
                roi_type = poly_roi_2.get("roi_kind") or (
                    "polygon" if "polygon" in poly_roi_2 else "circle picks"
                )
                print(f"  ROI B: YAML ({roi_type}) from {roi_yaml_suffix}")
            else:
                expected = fp2.stem
                if file_suffix_2:
                    suffix_no_ext = file_suffix_2.replace(".hdf5", "").replace(".h5", "")
                    if expected.endswith(suffix_no_ext):
                        expected = expected[: -len(suffix_no_ext)]
                print(
                    f"  ROI B: no YAML found ({expected}{roi_yaml_suffix}.yaml); "
                    "using data bounding box"
                )

        nn_data, nn_sim, neighbor_counts = process_file_pair(
            fp1,
            fp2,
            roi_1,
            roi_2,
            random_seed,
            output_dir_1,
            output_dir_2,
            neighbor_radius_nm=neighbor_radius_nm,
        )

        if neighbor_counts is not None:
            radius_label = f"{neighbor_counts['radius_nm']:g}".replace(".", "p")
            df_counts = build_neighbor_counts_dataframe(
                neighbor_counts["label_a"],
                neighbor_counts["label_b"],
                neighbor_counts["radius_nm"],
                neighbor_counts["data"],
                neighbor_counts["sim"],
            )
            counts_csv = (
                output_dir_1 / f"{base_name}_neighbor_counts_{radius_label}nm.csv"
            )
            df_counts.to_csv(counts_csv, index=False)
            n_a = len(neighbor_counts["data"]["a_cross_b"])
            n_b = len(neighbor_counts["data"]["b_cross_a"])
            print(
                f"  Saved neighbor counts: {counts_csv.name} "
                f"({n_a} points in {neighbor_counts['label_a']}, "
                f"{n_b} points in {neighbor_counts['label_b']})"
            )

        if nn_data.size == 0:
            continue

        # Data NN: one row per ROI-filtered point in file A.
        # Sim NN: one row per simulated point in file A (full HDF5 count).
        n_rows = max(nn_data.shape[0], nn_sim.shape[0])
        data_dict = {}

        for i in range(N_NEAREST_NEIGHBORS):
            data_dict[f"cross_nn_{i+1}_nm"] = _pad_1d(nn_data[:, i], n_rows)
            data_dict[f"sim_cross_nn_{i+1}_nm"] = _pad_1d(nn_sim[:, i], n_rows)

        df = pd.DataFrame(data_dict)

        csv_path = output_dir_1 / f"{base_name}_cross_nearest_neighbors.csv"

        df.to_csv(csv_path, index=False)

        print(f"  Saved CSV: {csv_path.name} ({n_rows} rows)")
        if nn_data.shape[0] != nn_sim.shape[0]:
            print(
                f"    cross_nn_* columns: {nn_data.shape[0]} rows (ROI-filtered data); "
                f"sim_cross_nn_* columns: {nn_sim.shape[0]} rows (simulation); "
                "extra rows padded with NaN"
            )



def main():
    parser = argparse.ArgumentParser(
        description="Nearest neighbor analysis and simulation for cluster centers in HDF5 files"
    )
    parser.add_argument(
        "--input", "-i",
        default=INPUT_DIR_1,
        help="Root directory to search for HDF5 files"
    )
    parser.add_argument(
        "--suffix", "-s",
        default=FILE_SUFFIX_1,
        help="File suffix to match (e.g., _filtered_in.hdf5). Use .hdf5 for all."
    )
    parser.add_argument(
        "--roi",
        default=None,
        help="ROI: path to mask image, or 'x_min,y_min,x_max,y_max' for bounding box"
    )
    parser.add_argument(
        "--roi-yaml-suffix",
        default=None,
        help="YAML ROI suffix (e.g. _ROI_picks) - match per-file ROI like cluster_density_v1"
    )
    parser.add_argument(
        "--output", "-o",
        default=None,
        help="Output folder for CSV and simulation HDF5 files. If not set, uses OUTPUT_FOLDER from config"
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=RANDOM_SEED,
        help="Random seed for simulation"
    )
    parser.add_argument(
        "--neighbor-radius-nm",
        type=float,
        default=None,
        help=(
            "Count partner-dataset points within this radius (nm) around each "
            "ROI-filtered point (A→B and B→A). Overrides NEIGHBOR_RADIUS_NM."
        ),
    )
    args = parser.parse_args()

    roi = ROI if args.roi is None else args.roi
    if isinstance(roi, str) and "," in roi:
        parts = [float(p.strip()) for p in roi.split(",")]
        if len(parts) == 4:
            roi = tuple(parts)

    roi_yaml = args.roi_yaml_suffix if args.roi_yaml_suffix is not None else ROI_YAML_SUFFIX
    neighbor_radius_nm = (
        args.neighbor_radius_nm
        if args.neighbor_radius_nm is not None
        else NEIGHBOR_RADIUS_NM
    )

    run_batch_analysis(
        input_dir_1=INPUT_DIR_1,
        file_suffix_1=FILE_SUFFIX_1,
        input_dir_2=INPUT_DIR_2,
        file_suffix_2=FILE_SUFFIX_2,
        roi=roi,
        roi_yaml_suffix=roi_yaml,
        output_folder_1=OUTPUT_FOLDER_1,
        output_folder_2=OUTPUT_FOLDER_2,
        random_seed=args.seed,
        neighbor_radius_nm=neighbor_radius_nm,
    )


if __name__ == "__main__":
    main()
