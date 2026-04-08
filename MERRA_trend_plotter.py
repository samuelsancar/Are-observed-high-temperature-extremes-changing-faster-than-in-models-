"""
Code to plot MERRA2 tailwidth trends
"""

import iris
import sys
import matplotlib.pyplot as plt
import iris.plot as iplt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.cm as mpl_cm
import numpy as np
import xarray as xr
from Lat_Land_mask_function import lat_land_mask
from Save_figure_function import save_figure
from ERA5_MKtest import significance_plot
from Map_features import map_features

cmap_rb = mpl_cm.get_cmap('coolwarm', 12)

t1 = 1980
t2 = 2023
sig = sys.argv[1] # either 'fdr' or 'pvals'
alpha = 0.05 # For MK_test
lat_minimum = -60
lat_maximum = 90

reanalysis = 'MERRA'
variable = 'daily_maximum_2m_temperature'
datadir = '/data/users/samuel.sancar/Heat_extremes/Trends/Reanalysis_trends/'
trend = iris.load_cube(f"{datadir}{reanalysis}_{variable}_tailwidth_trends_{t1}-{t2}.nc")
p87p5_trend = iris.load_cube(f"{datadir}{reanalysis}_{variable}_p87p5_trends_{t1}-{t2}.nc")
template = iris.load_cube('template_cube_1deg.nc')

# Apply latitude and land mask using function
model_mask, _ = lat_land_mask(trend, lat_minimum, lat_maximum)
p87p5_masked, _ = lat_land_mask(p87p5_trend, lat_minimum, lat_maximum)

# Extract grid point locations for stippling
lon = template.coord('longitude').points
lat = template.coord('latitude').points
lon2d, lat2d = np.meshgrid(lon, lat)

# Load significance mask (already masked for land and latitude)
if sig == 'fdr':
    statsig = iris.load_cube(f"{datadir}{reanalysis}_{variable}_tailwidth_fdr_all_{t1}-{t2}.nc")
    statsig_bool =  statsig.data.astype(bool)
    stipple_lats, stipple_lons = lat2d[statsig_bool], lon2d[statsig_bool]
    
else:
    statsig_ds = xr.open_dataset(f"{datadir}{reanalysis}_{variable}_tailwidth_p_{t1}-{t2}.nc")
    statsig = list(statsig_ds.data_vars.values())[0]
    sig_mask_np = statsig < alpha
    stipple_lats, stipple_lons = lat2d[sig_mask_np], lon2d[sig_mask_np]

neg_indices = np.where(p87p5_masked.data < 0)
neg_lats = lat2d[neg_indices]
neg_lons = lon2d[neg_indices]

# Plot trends
fig = plt.figure(figsize=(14, 8))
ax = plt.axes(projection=ccrs.Robinson(central_longitude=0))

map = iplt.pcolormesh(model_mask, cmap=cmap_rb, vmin=-0.4, vmax=0.4)

# Plot significant trends
ax.scatter(stipple_lons, stipple_lats, color='black', s=0.5, alpha=0.6, transform=ccrs.PlateCarree(central_longitude=0), label='Stippling')

# Overlay white dots where p87p5_trend is negative
ax.scatter(neg_lons, neg_lats, color='white', s=0.9, alpha=0.6, transform=ccrs.PlateCarree(), label='p87p5 trend < 0')

map_features('Rob', ax, -60, 90)

cbar = plt.colorbar(map, ax=ax, orientation='horizontal', pad=0.05, aspect=50)
cbar.ax.tick_params(labelsize=14)
cbar.set_label('Trend in Tailwidth (°C per decade)', fontsize=15)

plt.title(f"MERRA tail width trends ({t1}-{t2-1})", fontsize=18, weight='bold')

plt.tight_layout()
save_figure(plt, f'MERRA_taiwidth_trends_({t1}-{t2-1}({sig}))', 'Maps')
plt.show()