"""
22/09/25
---------
Plotting false discovery rate corrected P-values for significance stippling for Model tail width trends (alpha=0.1).

To run: python3 Model_FDR_correction.py

"""

import sys
import numpy as np
import os
import iris
import glob
from statsmodels.stats.multitest import fdrcorrection
from Lat_Land_mask_function import lat_land_mask

# Load p-value data
datadir = '/data/users/samuel.sancar/Heat_extremes/Trends/Trend_significance/P-vals/*_p_1958.nc'
Models = sorted(glob.glob(datadir, recursive=True))
lat_minimum = -60
lat_maximum =int(sys.argv[1])  # e.g., 90

# Output directories
outdir_all = "/data/users/samuel.sancar/Heat_extremes/Trends/Trend_significance/FDR_corrections_All/"
outdir_pos = "/data/users/samuel.sancar/Heat_extremes/Trends/Trend_significance/FDR_corrections_PositiveOnly/"
os.makedirs(outdir_all, exist_ok=True)
os.makedirs(outdir_pos, exist_ok=True)

for f in Models:
    print(f"Processing file: {f}")
    base_name = os.path.splitext(os.path.basename(f))[0]

    # Load p-value cube
    p_cube = iris.load_cube(f)
    pvals = p_cube.data

    # Load corresponding Kendall tau cube
    
    tau_file = f"/data/users/samuel.sancar/Heat_extremes/Trends/Trend_significance/P-vals/{base_name.replace('_p_1958', '_mktau_1958-2022')}.nc"
    tau_cube = iris.load_cube(tau_file)
    tau_vals = tau_cube.data


    # Apply lat_land_mask
    masked_cube, no_of_gridpoints = lat_land_mask(pvals, lat_minimum, lat_maximum)
    masked_data = masked_cube.data

    # Flatten valid values and apply FDR correction
    valid_mask = np.isfinite(masked_data)
    pvals_flat = masked_data[valid_mask]
    rejected, corrected_pvals = fdrcorrection(pvals_flat, alpha=0.1)

    # Create full mask with rejected values
    rejected_mask_full = np.full(masked_data.shape, False)
    rejected_mask_full[valid_mask] = rejected

    # Save full FDR mask (all trends)
    np.save(os.path.join(outdir_all, f"{base_name}_fdr_significant_all.npy"), rejected_mask_full)
    print(f"Saved all-trend FDR mask: {np.sum(rejected)} significant points")

    # Save positive trend FDR mask
    tau_masked = tau_vals[valid_mask]
    pos_mask = np.full(masked_data.shape, False)
    pos_mask[valid_mask] = np.logical_and(rejected, tau_masked > 0)
    np.save(os.path.join(outdir_pos, f"{base_name}_fdr_significant_positive.npy"), pos_mask)
    print(f"Saved positive-trend FDR mask: {np.sum(pos_mask)} significant points")

