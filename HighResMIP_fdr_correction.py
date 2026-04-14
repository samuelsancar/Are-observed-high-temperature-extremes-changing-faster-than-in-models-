"""
22/09/25
---------
Correcting for significant trends for false discovery rate  for HighResMIP members(alpha=0.1).

To run: python3 HighResMIP_FDR_correction.py

"""
import sys
import numpy as np
import os
import iris
import glob
from statsmodels.stats.multitest import fdrcorrection
from Lat_Land_mask_function import lat_land_mask
from iris import Constraint

t1 = sys.argv[1]

lat_minimum = -60
lat_maximum = 90

datadir = f'/data/users/samuel.sancar/Heat_extremes/Trends/Trend_significance/P-vals/*_p_{t1}.nc'
Models = sorted(glob.glob(datadir, recursive=True))

# Output directories
outdir_all = "/data/users/samuel.sancar/Heat_extremes/Trends/Trend_significance/FDR_corrections_All/"
outdir_pos = "/data/users/samuel.sancar/Heat_extremes/Trends/Trend_significance/FDR_corrections_PositiveOnly/"

os.makedirs(outdir_all, exist_ok=True)
os.makedirs(outdir_pos, exist_ok=True)

def _var_or_name(name):
    return Constraint(cube_func=lambda c: (getattr(c, 'var_name', None) == name) or (c.name() == name))

for f in Models:
    
    print(f"Processing file: {f}")
    base_name = os.path.splitext(os.path.basename(f))[0]
    save_basename = "_".join(f.split('/')[-1].split('_')[0:3])

    # Load p-value cube
    Stat_sig = iris.load_cube(f, _var_or_name('Tailwidth_trend_p-values'))
    pvals = Stat_sig.data
    
    # Number of valid gridpoints
    _, no_of_gridpoints = lat_land_mask(pvals, lat_minimum, lat_maximum)

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
    latitudes = Stat_sig.coord('latitude').points
    longitudes = Stat_sig.coord('longitude').points
    lon2d, lat2d = np.meshgrid(longitudes, latitudes)
    stipple_lats, stipple_lons = lat2d[rejected_mask_full], lon2d[rejected_mask_full]

    # Save all FDR significant trends
    fdr_cube = Stat_sig.copy()
    fdr_cube.data = rejected_mask_full.astype(np.int8)
    fdr_cube.long_name = 'fdr corrected trends'
    fdr_cube.var_name = 'fdr_corrected_trends'
    fdr_cube.units = '1'
    fdr_cube.attributes = {}
    iris.save(fdr_cube, os.path.join(outdir_all, f"{save_basename}_fdr_significant_all_{t1}.nc"))
    print(f"Saved all-trend FDR mask: {np.sum(rejected)} significant points")
    
    # Save FDR significant positive trends
    tau_path = f"/data/users/samuel.sancar/Heat_extremes/Trends/Trend_significance/P-vals/{base_name.replace(f'_p_{t1}', f'_mktau_{t1}-2022')}.nc"
    tau_cube = iris.load_cube(tau_path, _var_or_name("Tailwidth_mktau"))
    tau_vals = tau_cube.data
    tau_masked = tau_vals[valid_mask]

    pos_mask = np.full(pvals.shape, False)
    pos_mask[valid_mask] = np.logical_and(rejected, tau_masked > 0)

    pos_fdr_cube = Stat_sig.copy()
    pos_fdr_cube.data = pos_mask.astype(np.int8)
    pos_fdr_cube.long_name = 'fdr corrected positive trends'
    pos_fdr_cube.var_name = 'fdr_corrected_positive_trends'
    pos_fdr_cube.units = '1'
    pos_fdr_cube.attributes = {}
    iris.save(pos_fdr_cube, os.path.join(outdir_pos, f"{save_basename}_tailwidth_fdr_positive_{t1}.nc"))
    print(f"Saved positive-trend FDR mask: {np.sum(pos_mask)} significant points")