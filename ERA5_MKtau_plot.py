"""
24/07/25
---------
Plot trend in MannKendall (KendallTau) results for ERA5 tail width trends.
To run: python plot_tailwidth_trend.py [1958] [2023]

"""

import iris
import sys
import matplotlib.pyplot as plt
import iris.quickplot as qplt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.cm as mpl_cm
import xarray as xr
import numpy as np

# Define color map
cmap_rb = mpl_cm.get_cmap('coolwarm', 15)

# Read command-line arguments
reanalysis = 'era5'
t1 = int(sys.argv[1])
t2 = int(sys.argv[2])
variable = 'daily_maximum_2m_temperature'

# Load trend data
datadir = '/data/users/samuel.sancar/Heat_extremes/trends/Reanalysis_trends/'
trend = xr.open_dataset(f"{datadir}{reanalysis}_{variable}_tailwidth_mktau_{t1}-{t2}.nc")

# Download Significance stippling
statsig = np.load('/data/users/samuel.sancar/Heat_extremes/trends/Reanalysis_trends/era5_daily_maximum_2m_temperature_tailwidth_p_mask_1958-2023.npy')

# Assuming the trend variable is the first data variable in the dataset
data_var = list(trend.data_vars)[0]

lon = trend['longitude'].values
lat = trend['latitude'].values
lon2d, lat2d = np.meshgrid(lon, lat)

stipple_lats, stipple_lons = lat2d[statsig], lon2d[statsig]

# Create figure and axis
fig = plt.figure(figsize=(12, 6))
ax = plt.axes(projection=ccrs.PlateCarree())

# Plot using xarray
mesh = trend[data_var].plot.pcolormesh(ax=ax,cmap=cmap_rb,vmin=-0.3,vmax=0.3,add_colorbar=False,transform=ccrs.PlateCarree())

# Add map features
ax.coastlines(resolution='50m', linewidth=0.5)
ax.add_feature(cfeature.NaturalEarthFeature(category='cultural', name='admin_0_boundary_lines_land', facecolor='none', scale='50m', edgecolor='0', lw=.05))
ax.add_feature(cfeature.LAKES.with_scale('110m'), edgecolor='0', facecolor='none', lw=.05)
ax.add_feature(cfeature.STATES.with_scale('110m'), edgecolor='0', facecolor='none', lw=.2)
ax.add_feature(cfeature.RIVERS.with_scale('110m'), edgecolor='tab:blue', facecolor='none', lw=.4)
ax.add_feature(cfeature.LAND, facecolor='lightgray')
ax.add_feature(cfeature.OCEAN, facecolor='white')
ax.scatter(stipple_lons, stipple_lats,color='black', s=0.5, alpha=0.5, transform=ccrs.PlateCarree(), label='Stippling')

# Add gridlines
gl = ax.gridlines(draw_labels=True, linewidth=0.5, color='gray', alpha=0.5, linestyle='--')
gl.top_labels = False
gl.right_labels = False

# Add colorbar
cbar = plt.colorbar(mesh, orientation='vertical', pad=0.05, aspect=50)
cbar.set_label('Trend in Tail Width (°C per decade)', fontsize=12)

# Add title
plt.title(f"KendallTau in Tail Width\n{reanalysis} {variable.replace('_', ' ')} ({t1}-{t2})", fontsize=14, weight='bold')

# Layout
plt.tight_layout()
plt.show()



