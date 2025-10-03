"""
24/07/25
---------
Plotting false discovery rate corrected P-values for significance stippling for ERA5 tail width trends (alpha=0.1).
To run: python ERA5_MKstatsig_plot.py [1958] [2023]
"""

import sys
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.cm as mpl_cm
import xarray as xr
import numpy as np
import os
from statsmodels.stats.multitest import fdrcorrection
from Lat_Land_mask_function import lat_land_mask

# Colour map
cmap_rb = mpl_cm.get_cmap('RdBu_r', 15)

# Command-line arguments
reanalysis = 'era5'
variable = 'daily_maximum_2m_temperature'
t1 = 1958
t2 = 2023
lat_minimum = -60
lat_maximum = int(sys.argv[1])
alpha = float(sys.argv[2])

datadir = '/data/users/samuel.sancar/Heat_extremes/Trends/Reanalysis_trends/'
Stat_sig = xr.open_dataset(f"{datadir}{reanalysis}_{variable}_tailwidth_p_{t1}-{t2}.nc")
pvals = Stat_sig["P-values"].values

# Apply lat_land_mask before FDR correction
masked_cube, no_of_gridpoints = lat_land_mask(pvals, lat_minimum, lat_maximum)
masked_data = masked_cube.data

# Flatten valid values and apply FDR correction
valid_mask = np.isfinite(masked_data)
pvals_flat = masked_data[valid_mask]
rejected, corrected_pvals = fdrcorrection(pvals_flat, alpha)

# Create full mask with rejected values
rejected_mask_full = np.full(masked_data.shape, False)
rejected_mask_full[valid_mask] = rejected
number_of_rejected = np.sum(rejected)

print(f"Number of significant grid points after FDR correction: {number_of_rejected} out of {no_of_gridpoints} valid grid points")

# Extract stippling coordinates
latitudes = masked_cube.coord('latitude').points
longitudes = masked_cube.coord('longitude').points
lon2d, lat2d = np.meshgrid(longitudes, latitudes)
stipple_lats, stipple_lons = lat2d[rejected_mask_full], lon2d[rejected_mask_full]

# Save significance mask to file
outpath = "/data/users/samuel.sancar/Heat_extremes/Trends/Reanalysis_trends"
outdir_pos = "/data/users/samuel.sancar/Heat_extremes/Trends/Reanalysis_trends/"

filename = f"{reanalysis}_{variable}_tailwidth_fdr_{t1}-{t2}.npy"
np.save(os.path.join(outpath, filename), rejected_mask_full)

# Save FDR significant positive trends
tau_data = xr.open_dataset(f"{datadir}{reanalysis}_{variable}_tailwidth_mktau_1958-2023.nc")
tau_vals = tau_data["P-values"].values
tau_masked = tau_vals[valid_mask]
pos_mask = np.full(masked_data.shape, False)
pos_mask[valid_mask] = np.logical_and(rejected, tau_masked > 0)
filename_pos = f"{reanalysis}_{variable}_tailwidth_fdr_positive_{t1}-{t2}.npy"
np.save(os.path.join(outdir_pos, filename_pos), pos_mask)
print(f"Saved positive-trend FDR mask: {np.sum(pos_mask)} significant points")

fig = plt.figure(figsize=(12, 6))
ax = plt.axes(projection=ccrs.PlateCarree())
ax.set_extent([-180, 180, -60, 90], crs=ccrs.PlateCarree())

ax.scatter(stipple_lons, stipple_lats, color='green', s=0.5, alpha=0.5, transform=ccrs.PlateCarree(), label='Stippling')
ax.coastlines(resolution='50m', linewidth=1)
ax.add_feature(cfeature.LAND, facecolor='white')
ax.add_feature(cfeature.OCEAN, facecolor='white')

gl = ax.gridlines(draw_labels=True, linewidth=0.5, color='gray', alpha=0.5, linestyle='--')
gl.top_labels = False
gl.right_labels = False

plt.title(f"Significant tailwidth trends after FDR correction ($alpha$ < 0.1) for ERA5 reanalysis (1958-2022)", fontsize=14, weight='bold')
plt.tight_layout()
plt.show()

