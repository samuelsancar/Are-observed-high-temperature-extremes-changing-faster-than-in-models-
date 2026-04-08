"""
24/07/25
---------
Plotting false discovery rate corrected P-values for significance stippling for ERA5 tail width trends (alpha=0.1).
To run: python ERA5_MKstatsig_plot.py [1958] [2023]
"""

import iris
import matplotlib.cm as mpl_cm
import numpy as np
import os
import sys
from statsmodels.stats.multitest import fdrcorrection
from Lat_Land_mask_function import lat_land_mask
from iris import Constraint

cmap_rb = mpl_cm.get_cmap('RdBu_r', 15)

t1 = sys.argv[1]
lat_minimum = -60
lat_maximum = 90
alpha = 0.1

reanalysis = 'era5'
variable = 'daily_maximum_2m_temperature'
dir = '/data/users/samuel.sancar/Heat_extremes/Trends/Reanalysis_trends/'
Stat_sig_name = f"{dir}{reanalysis}_{variable}_tailwidth_p_{t1}-2023.nc"
Mk_tau_name = f"{dir}{reanalysis}_{variable}_tailwidth_mktau_{t1}-2023.nc"

def _var_or_name(name):
    return Constraint(cube_func=lambda c: (getattr(c, 'var_name', None) == name) or (c.name() == name))

# Load p-value cube
Stat_sig = iris.load_cube(Stat_sig_name, _var_or_name('Tailwidth_trend_p-values'))
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
filename = f"{reanalysis}_{variable}_tailwidth_fdr_all_{t1}-2023.nc"
iris.save(fdr_cube, os.path.join(dir, filename))
print(f"Saved all-trend FDR mask: {np.sum(rejected)} significant points")

# Save FDR significant positive trends
tau_cube = iris.load_cube(Mk_tau_name)
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
filename_pos = f"{reanalysis}_{variable}_tailwidth_fdr_positive_{t1}-2023.nc"
iris.save(pos_fdr_cube, os.path.join(dir, filename_pos))
print(f"Saved positive-trend FDR mask: {np.sum(pos_mask)} significant points")

