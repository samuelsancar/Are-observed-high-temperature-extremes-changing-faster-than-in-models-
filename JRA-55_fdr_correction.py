"""

Code to calculate FDR correction for trends in tailwidths for JRA-55 reanalysis data.

To run: python JRA-55_fdr_correction.py

"""
import sys
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.cm as mpl_cm
import xarray as xr
import numpy as np
import os
import iris
import iris.coords
import iris.cube
from statsmodels.stats.multitest import fdrcorrection
from Lat_Land_mask_function import lat_land_mask

# Colour map
cmap_rb = mpl_cm.get_cmap('RdBu_r', 15)

# Command-line arguments
reanalysis = 'jra55'
t1 = int(sys.argv[1])
t2 = int(sys.argv[2])
variable = '_minmax_surf.015_tmax.reg_tl319.'

# Load p-value data
datadir = '/data/users/samuel.sancar/Heat_extremes/Trends/Reanalysis_trends/'
Stat_sig = xr.open_dataset(f"{datadir}{reanalysis}_{variable}_tailwidth_p_{t1}-{t2}.nc")
pvals = Stat_sig["__xarray_dataarray_variable__"].values

# Apply lat_land_mask before FDR correction
masked_cube, no_of_gridpoints = lat_land_mask(pvals, -60, 90)
masked_data = masked_cube.data

# Flatten valid values and apply FDR correction
valid_mask = np.isfinite(masked_data)
pvals_flat = masked_data[valid_mask]
rejected, corrected_pvals = fdrcorrection(pvals_flat, alpha=0.1)

# Create full mask with rejected values
rejected_mask_full = np.full(masked_data.shape, False)
rejected_mask_full[valid_mask] = rejected

# Extract stippling coordinates
latitudes = masked_cube.coord('latitude').points
longitudes = masked_cube.coord('longitude').points
lon2d, lat2d = np.meshgrid(longitudes, latitudes)
stipple_lats, stipple_lons = lat2d[rejected_mask_full], lon2d[rejected_mask_full]

# Save significance mask to file
outpath = "/data/users/samuel.sancar/Heat_extremes/Trends/Reanalysis_trends"
filename = f"{reanalysis}_{variable}_tailwidth_fdr_{t1}-{t2}.npy"
np.save(os.path.join(outpath, filename), rejected_mask_full)

# Plotting
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

plt.title(f"Significant tailwidth trend after FDR correction ($alpha$ < 0.1) for \n{reanalysis} {variable.replace('_', ' ')} ({t1}-{t2})", fontsize=14, weight='bold')
plt.tight_layout()
plt.show()