"""
15/05/26
--------

Plot CANARI model mean tail widths trends and std.

To run: python CANARI_mean_trend_plotter.py
"""

import iris
import sys
import matplotlib.pyplot as plt
import iris.plot as iplt
import cartopy.crs as ccrs
import matplotlib.cm as mpl_cm
import datetime
import numpy as np
from Lat_Land_mask_function import lat_land_mask
from Map_features import map_features
from Save_figure_function import save_figure

# Define colour map
cmap_rb = mpl_cm.get_cmap('coolwarm', 15)

t1 = 1958
t2 = 2023

alpha =  0.05

# Latitude limits for plotting
lat_minimum = -60
lat_maximum = 90

# Load trends data
template = iris.load_cube('template_cube_1deg.nc')
filename = f"CANARI_em_tailwidth_trend_{t1}-{t2-1}.nc"
datadir = f'/data/users/samuel.sancar/Heat_extremes/Trends/CANARI_trends/'
cube = iris.load_cube(datadir + filename)
p87p5_trend = iris.load_cube(datadir + filename.replace('tailwidth','p87p5'))

cube = iris.util.squeeze(cube)
p87p5_trend = iris.util.squeeze(p87p5_trend)

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
lon = template.coord('longitude').points
lat = template.coord('latitude').points
lon2d, lat2d = np.meshgrid(lon, lat)

# Extract grid point locations where p87p5 trend is negative for overlaying white stippling
neg_mask = np.ma.filled(p87p5_masked.data < 0, False)
neg_indices = np.where(neg_mask)

neg_lats = lat2d[neg_indices]
neg_lons = lon2d[neg_indices]

# Plot trends
fig = plt.figure(figsize=(12, 6))
ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=0))

map_features('Detail', ax, lat_minimum, lat_maximum)

map = iplt.pcolormesh(masked_cube, cmap=cmap_rb, vmin=-0.2, vmax=0.2)

cbar = plt.colorbar(map, ax=ax, orientation='horizontal', pad=0.05, aspect=50)
cbar.ax.tick_params(labelsize=16)
cbar.set_label('Trend in Tailwidth (°C per decade)', fontsize=12)

plt.title(f"CANARI ensemble mean tail width trend ({t1}-t2-1)", fontsize=14, weight='bold')
plt.tight_layout()
save_figure(plt, f'CANARI_model_mean_tailwidth_trends', 'Maps')
plt.show()

# Plot std
std_filename = f"CANARI_std_{t1}-{t2-1}.nc"
std_datadir = f'/data/users/samuel.sancar/Heat_extremes/Trends/CANARI_trends/'
std_cube = iris.load_cube(std_datadir + std_filename)

masked_std, _ = lat_land_mask(std_cube, template, lat_minimum, lat_maximum)

fig = plt.figure(figsize=(12, 6))
ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=0))

map_features('Detail', ax, lat_minimum, lat_maximum)

map = iplt.pcolormesh(masked_std, cmap=cmap_rb, vmin=0, vmax=0.1)

# Add colourbar
cbar = plt.colorbar(map, ax=ax, orientation='horizontal', pad=0.05, aspect=50)
cbar.ax.tick_params(labelsize=16)
cbar.set_label('Std', fontsize=12)

plt.title(f"CANARI std ({t1}-{t2-1})", fontsize=14, weight='bold')
plt.tight_layout()
save_figure(plt, f'CANARI_std', 'Maps')
plt.show()

# Load full tail width data with time dimension
full_cube = iris.load_cube(f'/data/users/samuel.sancar/Heat_extremes/Tail_widths/canari/CANARI_em_tailwidths_{t1}-{t2-1}.nc')

# Choose a grid point (default is over London)
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

# Extract trend value from masked_cube
trend = cube.data[lon_box, lon_box]

# Generate trend line using the same x-axis scale
trend_line = trend/10 * np.arange(len(time_series)) + time_series[0]

# Plot the time series and trend line
plt.figure(figsize=(10, 5))
plt.plot(np.arange(1958,2023), time_series, marker='o', linestyle='-', label='Tail Width')

plt.title(f"Tail Width Over Time for CANARI member CANARI_em at Grid Point ({lat_index}, {lon_index})", fontsize=14)
plt.xlabel("Time")
plt.ylabel("Tail Width (°C)")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()
print("Trend value:", trend)