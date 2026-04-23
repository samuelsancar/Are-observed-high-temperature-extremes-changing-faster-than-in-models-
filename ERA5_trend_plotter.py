"""
17/07/25
---------
Plot trend in tail width calculated by ERA5_trend_calc.py

Used to produce ERA5 plot for figure 1.

Note: ERA5 significant trends are calculated in ERA5_MKtest.py and saved as netCDF files. These must be calculated first for code to work. 
These are loaded here to add stippling to the trend map.

To run: 'python ERA5_trend_plotter.py [fdr|p]' where fdr or p determines whether to use the 
FDR-corrected significance mask or the raw p-value mask for stippling.
"""

import iris
import sys
import matplotlib.pyplot as plt
import iris.plot as iplt
import cartopy.crs as ccrs
import matplotlib.cm as mpl_cm
import numpy as np

from Lat_Land_mask_function import lat_land_mask
from Save_figure_function import save_figure
from Map_features import map_features

# Define colormap
cmap_rb = mpl_cm.get_cmap('coolwarm', 12)

# Beginning & end year
t1 = 1958
t2 = 2023

sig = sys.argv[1] # either 'fdr' or 'pvals'
alpha = 0.05 # Test level for significance testing 
  
# Latitude limits for plotting
lat_minimum = -60
lat_maximum = 90

reanalysis = 'era5'
variable = 'daily_maximum_2m_temperature'
datadir = '/data/users/samuel.sancar/Heat_extremes/Trends/Reanalysis_trends/'

template = iris.load_cube('template_cube_1deg.nc')
trend = iris.load_cube(f"{datadir}{reanalysis}_{variable}_tailwidth_trends_{t1}-{t2-1}.nc")
p87p5_trend = iris.load_cube(f"{datadir}{reanalysis}_{variable}_p87p5_trends_{t1}-{t2-1}.nc")

# Guess bounds
lon_coord = template.coord('longitude')
if not lon_coord.has_bounds():
    lon_coord.guess_bounds()

lat_coord = template.coord('latitude')
if not lat_coord.has_bounds():
    lat_coord.guess_bounds()

# Apply latitude and land mask using function
model_mask, _ = lat_land_mask(trend, template, lat_minimum, lat_maximum)
p87p5_masked, _ = lat_land_mask(p87p5_trend, template, lat_minimum, lat_maximum)

# Extract grid point locations for stippling
lon = template.coord('longitude').points
lat = template.coord('latitude').points
lon2d, lat2d = np.meshgrid(lon, lat)

# Load significance mask (already masked for land and latitude)
if sig == 'fdr':
    statsig = iris.load_cube(f"{datadir}{reanalysis}_{variable}_fdr_corrected_p_{t1}-{t2-1}.nc")
    data = statsig.data  # masked array
    stipple_mask = (~data.mask) & data.astype(bool)
    stipple_lats, stipple_lons = lat2d[stipple_mask], lon2d[stipple_mask]
else:
    statsig = iris.load_cube(f"{datadir}{reanalysis}_{variable}_tailwidth_p_{t1}-{t2-1}.nc")
    statsig_masked, _ = lat_land_mask(statsig, template, lat_minimum, lat_maximum)
    sig_mask_np = np.ma.filled(statsig_masked.data < alpha, False)
    stipple_lats, stipple_lons = lat2d[sig_mask_np], lon2d[sig_mask_np]
    
# Extract grid point locations where p87p5 trend is negative for overlaying white stippling
neg_mask = np.ma.filled(p87p5_masked.data < 0, False)
neg_indices = np.where(neg_mask)

neg_lats = lat2d[neg_indices]
neg_lons = lon2d[neg_indices]

# Plot trends
fig = plt.figure(figsize=(12, 6))
ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=0))

map_features('Detail', ax, lat_minimum, lat_maximum)

map = iplt.pcolormesh(model_mask, cmap=cmap_rb, vmin=-0.3, vmax=0.3)

# Plot significant trends
ax.scatter(stipple_lons, stipple_lats, color='black', s=0.3, alpha=0.6, transform=ccrs.PlateCarree(), label='Stippling')

# Pot white stippling where p87p5_trend is negative
ax.scatter(neg_lons, neg_lats, color='white', s=0.9, alpha=0.6, transform=ccrs.PlateCarree(), label='p87p5 trend < 0')

cbar = plt.colorbar(map, ax=ax, orientation='horizontal', pad=0.05, aspect=50)
cbar.ax.tick_params(labelsize=16)
cbar.set_label('Trend in Tail width (°C per decade)', fontsize=17)

plt.title(f"ERA5", fontsize=17, weight='bold')
plt.tight_layout()
save_figure(plt, f'ERA5_taiwidth_trends_({t1}-{t2-1}({sig}))', 'Maps')
plt.show()

