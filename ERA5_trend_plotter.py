"""
17/07/25
---------
Plot trend in tail width calculated by ERA5_tailwidth.py
To run: python plot_tailwidth_trend.py [1958] [2023] [fdr|p]
"""

import iris
import sys
import os
import matplotlib.pyplot as plt
import iris.plot as iplt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.cm as mpl_cm
import numpy as np
import xarray as xr
from Lat_Land_mask_function import lat_land_mask

# Define colour map
cmap_rb = mpl_cm.get_cmap('coolwarm', 15)

# Read command-line arguments
reanalysis = 'era5'
t1 = int(sys.argv[1])
t2 = int(sys.argv[2])
sig = sys.argv[3]
variable = 'daily_maximum_2m_temperature'

# Load trend data
datadir = '/data/users/samuel.sancar/Heat_extremes/Trends/Reanalysis_trends/'
trend = iris.load_cube(f"{datadir}{reanalysis}_{variable}_tailwidth_trends_{t1}-{t2}.nc")
template = iris.load_cube('template_cube_1deg.nc')

# Apply latitude and land mask using function
model_mask, _ = lat_land_mask(trend, -60, 90)

# Extract grid point locations for stippling
lon = template.coord('longitude').points
lat = template.coord('latitude').points
lon2d, lat2d = np.meshgrid(lon, lat)

# Load significance mask (already masked for land and latitude)
if sig == 'fdr':
    statsig = xr.DataArray(
        data=np.load(f"{datadir}{reanalysis}_{variable}_tailwidth_fdr_{t1}-{t2}.npy"),
        dims=["latitude", "longitude"],
        coords={"latitude": lat, "longitude": lon},
        name="significance_mask"
    )
    stipple_lats, stipple_lons = lat2d[statsig], lon2d[statsig]
else:
    statsig_ds = xr.open_dataset(f"{datadir}{reanalysis}_{variable}_tailwidth_p_{t1}-{t2}.nc")
    statsig = list(statsig_ds.data_vars.values())[0]
    sig_mask_np = statsig.values < 0.05
    stipple_lats, stipple_lons = lat2d[sig_mask_np], lon2d[sig_mask_np]

# Create figure and axis
fig = plt.figure(figsize=(12, 6))
ax = plt.axes(projection=ccrs.PlateCarree())
ax.set_extent([-180, 180, -60, 90], crs=ccrs.PlateCarree())

# Plot using Iris
iplt.pcolormesh(model_mask, cmap=cmap_rb, vmin=-0.3, vmax=0.3)

# Add map features
ax.coastlines(resolution='50m', linewidth=0.5)
ax.add_feature(cfeature.NaturalEarthFeature(category='cultural', name='admin_0_boundary_lines_land', facecolor='none', scale='50m', edgecolor='0', lw=.05))
ax.add_feature(cfeature.LAKES.with_scale('110m'), edgecolor='0', facecolor='none', lw=.05)
ax.add_feature(cfeature.STATES.with_scale('110m'), edgecolor='0', facecolor='none', lw=.2)
ax.add_feature(cfeature.RIVERS.with_scale('110m'), edgecolor='tab:blue', facecolor='none', lw=.4)
ax.add_feature(cfeature.LAND, facecolor='lightgray')
ax.add_feature(cfeature.OCEAN, facecolor='white')

# Plot stippling
ax.scatter(stipple_lons, stipple_lats, color='black', s=0.5, alpha=0.5, transform=ccrs.PlateCarree(), label='Stippling')

# Add gridlines
gl = ax.gridlines(draw_labels=True, linewidth=0.5, color='gray', alpha=0.5, linestyle='--')
gl.top_labels = False
gl.right_labels = False

# Add colourbar
cbar = plt.colorbar(ax.collections[0], orientation='vertical', pad=0.05, aspect=50)
cbar.set_label('Trend in Tail Width (°C per decade)', fontsize=12)

# Add title
plt.title(f"Trend in Tail Width\n{reanalysis} {variable.replace('_', ' ')} ({t1}-{t2})", fontsize=14, weight='bold')

# Layout
plt.tight_layout()
plt.savefig("/home/users/samuel.sancar/Documents/Heat_Extremes/CSSP_Images/ERA5_fdr_corrected.png", dpi = 500)
plt.show()




