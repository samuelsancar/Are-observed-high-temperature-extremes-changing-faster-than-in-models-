"""
22/09/25
---------
Plotting false discovery rate corrected P-values for significance stippling for CANARI tail width trends (alpha=0.1).

To run: python3 Model_FDR_correction.py

"""

import sys
import xarray as xr
import numpy as np
import os
import iris
import glob
from statsmodels.stats.multitest import fdrcorrection
from Lat_Land_mask_function import lat_land_mask

datadir = '/data/users/samuel.sancar/Heat_extremes/Trends/Trend_significance/CANARI/P-vals/*_p_1958.nc'
Members = sorted(glob.glob(datadir, recursive=True))
lat_minimum = -60
lat_maximum = 90

# Output directories
outdir_all = "/data/users/samuel.sancar/Heat_extremes/Trends/Trend_significance/CANARI/FDR_corrections_All/"
outdir_pos = "/data/users/samuel.sancar/Heat_extremes/Trends/Trend_significance/CANARI/FDR_corrections_PositiveOnly/"
os.makedirs(outdir_all, exist_ok=True)
os.makedirs(outdir_pos, exist_ok=True)


def convert_to_iris(rejected_mask_full):
    
    fdr_cube = iris.cube.Cube(rejected_mask_full.astype(np.int8), long_name='fdr corrected trends',
                                var_name='fdr_corrected_trends',units='1',attributes={})

    # Add lat/lon coords from dataset
    latitudes = Stat_sig['latitude'].values
    longitudes = Stat_sig['longitude'].values

    fdr_cube.add_dim_coord(iris.coords.DimCoord(latitudes, standard_name='latitude', units='degrees'),0)
    fdr_cube.add_dim_coord(iris.coords.DimCoord(longitudes, standard_name='longitude', units='degrees'),1)
        
    return fdr_cube

for f in Members:
    
    print(f"Processing file: {f}")
    base_name = os.path.splitext(os.path.basename(f))[0]

    # Load p-value cube
    Stat_sig = xr.open_dataset(f)
    pvals = Stat_sig['Tailwidth_trend_p-values'].values
    
    # Number of valid gridpoints
    _ , no_of_gridpoints = lat_land_mask(pvals, lat_minimum, lat_maximum)

    # Flatten valid values and apply FDR correction
    valid_mask = np.isfinite(pvals)
    pvals_flat = pvals[valid_mask]
    rejected, corrected_pvals = fdrcorrection(pvals_flat, alpha=0.1)

    # Create full mask with rejected values
    rejected_mask_full = np.full(pvals.shape, False)
    rejected_mask_full[valid_mask] = rejected
    number_of_rejected = np.sum(rejected)

    print(f"Number of significant grid points after FDR correction: {number_of_rejected} out of {no_of_gridpoints} valid grid points")

    # Extract stippling coordinates for plot
    latitudes = Stat_sig['latitude'].values
    longitudes = Stat_sig['longitude'].values
    lon2d, lat2d = np.meshgrid(longitudes, latitudes)
    stipple_lats, stipple_lons = lat2d[rejected_mask_full], lon2d[rejected_mask_full]

    # Save all FDR significant trends
    fdr_cube = convert_to_iris(rejected_mask_full)
    iris.save(fdr_cube, os.path.join(outdir_all, f"{base_name}_fdr_significant_all.nc"))
    print(f"Saved all-trend FDR mask: {np.sum(rejected)} significant points")
    
    # Save FDR significant positive trends
    tau_data = xr.open_dataset(f"/data/users/samuel.sancar/Heat_extremes/Trends/Trend_significance/CANARI/P-vals/{base_name.replace('_p_1958', '_mktau_1958-2022')}.nc")
    tau_vals = tau_data["Tailwidth_mktau"].values
    tau_masked = tau_vals[valid_mask]

    pos_mask = np.full(pvals.shape, False)
    pos_mask[valid_mask] = np.logical_and(rejected, tau_masked > 0)

    pos_fdr_cube = convert_to_iris(pos_mask)
    iris.save(pos_fdr_cube, os.path.join(outdir_pos, f"{base_name}_tailwidth_fdr_positive.nc"))
    print(f"Saved positive-trend FDR mask: {np.sum(pos_mask)} significant points")