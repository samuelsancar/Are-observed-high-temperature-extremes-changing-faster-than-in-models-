"""
1/04/26
--------

Plot HighResMIP multi-model mean tail width trend and inter member std.

Plots MMM tail width trend calculated, by average multi-model ensemble each year for years 1958-2022, 
then calculating the trend across this time period at each grid point.
Produces std plot for figure S1a in supplementary material.

To run: python HighResMIP_mmm_trend_plotter_2.0.py [fdr|p]
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

# Define colour map
cmap_rb = mpl_cm.get_cmap('coolwarm', 15)

t1 = 1958
t2= 2023

alpha = 0.05

lat_minimum = -60
lat_maximum = 90

# Command-line arguments
sig = sys.argv[1]  # 'fdr' or 'p'
datadir = f'/data/users/samuel.sancar/Heat_extremes/Trends/'

template = iris.load_cube('template_cube_1deg.nc')
filename = f"HighResMIP_MMM_tailwidth_trend_1958-2022.nc"
cube = iris.load_cube(datadir + filename)
p87p5_trend = iris.load_cube(datadir + filename.replace('tailwidth','p87p5'))

# Guess bounds
lon_coord = template.coord('longitude')
if not lon_coord.has_bounds():
    lon_coord.guess_bounds()

lat_coord = template.coord('latitude')
if not lat_coord.has_bounds():
    lat_coord.guess_bounds()

# Apply land mask
masked_cube, _ = lat_land_mask(cube, template, lat_minimum, lat_maximum)
p87p5_masked, _ = lat_land_mask(p87p5_trend, template, lat_minimum, lat_maximum)

# Extract grid point locations for stippling
lon = cube.coord('longitude').points
lat = cube.coord('latitude').points
lon2d, lat2d = np.meshgrid(lon, lat)

# Load significance mask (already masked for land and latitude)
if sig == 'fdr':
    sig_datadir = f'/data/users/samuel.sancar/Heat_extremes/Trends/Trend_significance/FDR_corrections/'
    statsig = iris.load_cube(f"{sig_datadir}HighResMIP_MMM_fdr_corrected_p_{t1}-{t2-1}.nc")
    data = statsig.data  # masked array
    stipple_mask = (~data.mask) & data.astype(bool)
    stipple_lats, stipple_lons = lat2d[stipple_mask], lon2d[stipple_mask]
else:
    sig_datadir = f'/data/users/samuel.sancar/Heat_extremes/Trends/Trend_significance/P-vals/'
    statsig = iris.load_cube(f"{sig_datadir}HighResMIP_MMM_p_{t1}-{t2-1}.nc")
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

map = iplt.pcolormesh(masked_cube, cmap=cmap_rb, vmin=-0.15, vmax=0.15)

# Plot significant trends
#ax.scatter(stipple_lons, stipple_lats, color='black', s=0.6, alpha=0.6, transform=ccrs.PlateCarree(), label='Stippling')

# Pot grey stippling where p87p5_trend is negative
#ax.scatter(neg_lons, neg_lats, color="grey", s=0.4, alpha=0.6, marker = 's', transform=ccrs.PlateCarree(), label='p87p5 trend < 0')

cbar = plt.colorbar(map, ax = ax, orientation='horizontal', pad=0.05, aspect=50)
cbar.ax.tick_params(labelsize=12)
cbar.set_label('Trend in Tailwidth (°C per decade)', fontsize=14)
    
plt.title(f"a) HighResMIP multi-model mean tail width trends", fontsize=17, weight='bold')
plt.tight_layout()
save_figure(plt, f'HighResMIP_MMM_tailwidths_{sig}', 'Maps')
plt.show()

# Plot std
std = iris.load_cube(f'/data/users/samuel.sancar/Heat_extremes/Trends/HighResMIP_std_{t1}-{t2-1}.nc')
std_masked, _ = lat_land_mask(std, template, lat_minimum, lat_maximum)

fig = plt.figure(figsize=(10, 6))
ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=0))

map_features('Detail', ax, lat_minimum, lat_maximum)

map = iplt.pcolormesh(std_masked, cmap='Reds', vmin=0, vmax=0.15)

#cbar = plt.colorbar(map, ax = ax, orientation='horizontal', pad=0.05, aspect=50)
#cbar.ax.tick_params(labelsize=12)
#cbar.set_label('Standard deviation (°C per decade)', fontsize=14)
    
plt.title(f"b) HRMIP standard deviation in tail width trends", fontsize=17, weight='bold')
plt.tight_layout()
save_figure(plt, f'HighResMIP_std', 'Maps')
plt.show()

# Plot time series at a gridpoint (London) to check that trends make sense. 
# Extract from full data cube with time dimension, not the trend cube.

# Load full tail width data with time dimension
datadir = f'/data/users/samuel.sancar/Heat_extremes/Tail_widths/'
full_data_filename = f"HighResMIP_MMM_tailwidths_1958-2022.nc"
full_cube = iris.load_cube(datadir + full_data_filename)

# Choose a grid point (London is default)
lat_index = 55
lon_index = 37

longitudes = cube.coord('longitude')
latitudes = cube.coord('latitude')
 
lon_box = longitudes.nearest_neighbour_index(lat_index)
lat_box = latitudes.nearest_neighbour_index(lon_index)

longitudes1 = full_cube.coord('longitude')
latitudes1 = full_cube.coord('latitude')
 
lon_box1 = longitudes1.nearest_neighbour_index(lat_index)
lat_box1 = latitudes1.nearest_neighbour_index(lon_index)

# Extract time series and time points
time_series = full_cube.data[:, lat_box1, lon_box1]

trend = cube.data[lon_box, lon_box]
#trend_line = trend/10 * np.arange(len(time_series)) + time_series[0]

plt.figure(figsize=(10, 5))
plt.plot(np.arange(1958,2023), time_series, marker='o', linestyle='-', label='Tail Width')
plt.title(f"Tail Width timeseries for \nHighResMIP mmm\n at Grid Point {lat_index}, {lon_index} (London)", fontsize=14)
plt.xlabel("Time")
plt.ylabel("Tail Width (°C)")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()
print("Trend value:", trend)