""""
In progress

Example code to calculate average trend in daily tasmax over all models, for 
years 1958-2022 (combines historical and future model runs).


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


datadir = '/data/users/samuel.sancar/Heat_extremes/trends/Coupled_model_trends/*'
outdir = '/data/users/samuel.sancar/Heat_extremes/trends/Analysis_trends'

# Define color map
cmap_rb = mpl_cm.get_cmap('RdBu_r', 15)

trends = xr.open_mfdataset(datadir, combine='nested',compat='override')

print(trends)
avg_trend = trends.mean(dim=['latitude','longitude'])

# Clear encoding for all coordinates to avoid xarray datetime bug
for coord in avg_trend.coords:
    avg_trend[coord].encoding.clear()
avg_trend.to_netcdf(outdir + '/avg_trend.nc')

# sea mask
lsm = xr.open_dataset('lsm_1279l4_0.1x0.1.grb_v4_unpack.nc')['lsm'].squeeze()
lsm_interp = lsm.interp(latitude=trends.latitude, longitude=trends.longitude, method="nearest")

# Apply mask:land (>0.5), sea (0)
masked_trends = trends.where(lsm_interp > 0.5)

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

plt.title(f"Average trend in Tail Widths over all coupled models (1958-2022)", fontsize=14, weight='bold')

plt.tight_layout()
plt.show()