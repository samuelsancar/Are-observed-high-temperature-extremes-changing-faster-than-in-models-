""""
30/07/2025

Code to calculate difference in daily tasmax tailwidth trends between ERA5 and model, for years 1958-2022.

To run: python Obs-Coupled_diff.py


"""

import iris
import numpy as np
import sys
import iris.analysis.stats
import matplotlib.pyplot as plt
import iris.plot as iplt
import iris.quickplot as qplt
import xarray as xr
import glob
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.cm as mpl_cm


datadir1 = '/data/users/samuel.sancar/Heat_extremes/trends/Analysis_trends/avg_trend.nc'
datadir2 = '/data/users/samuel.sancar/Heat_extremes/trends/Reanalysis_trends/era5_daily_maximum_2m_temperature_tailwidth_trends_1958-2023.nc'


# Define color map
cmap_rb = mpl_cm.get_cmap('RdBu_r', 15)

avg_trend = xr.open_dataset(datadir1)
trends = xr.open_dataset(datadir2) > 0 

difference = trends - avg_trend

# Download sea mask
lsm = xr.open_dataset('lsm_1279l4_0.1x0.1.grb_v4_unpack.nc')['lsm'].squeeze()
lsm_interp = lsm.interp(latitude=difference.latitude, longitude=difference.longitude, method="nearest")

# Apply mask: keep land (1), mask sea (0)

masked_trends = difference.where(lsm_interp > 0.5)

fig = plt.figure(figsize=(12, 6))
ax = plt.axes(projection=ccrs.PlateCarree())

data_var = list(masked_trends.data_vars)[0]

# Plot using xarray
mesh = masked_trends[data_var].plot.pcolormesh(ax=ax,cmap=cmap_rb,vmin=-0.4,vmax=0.4,add_colorbar=False,transform=ccrs.PlateCarree())

ax.coastlines(resolution='50m', linewidth=1)
ax.add_feature(cfeature.BORDERS, linestyle=':')
ax.add_feature(cfeature.LAND, facecolor='lightgray')
ax.add_feature(cfeature.OCEAN, facecolor='white')

gl = ax.gridlines(draw_labels=True, linewidth=0.5, color='gray', alpha=0.5, linestyle='--')
gl.top_labels = False
gl.right_labels = False

cbar = plt.colorbar(mesh, orientation='horizontal', pad=0.05, aspect=50)
cbar.set_label('Trend in Tail Width (°C per decade)', fontsize=12)

plt.title(f"Difference between Average coupled model trend and ERA5 trend in Tail Widths (1958-2022)", fontsize=14, weight='bold')

plt.tight_layout()
plt.show()