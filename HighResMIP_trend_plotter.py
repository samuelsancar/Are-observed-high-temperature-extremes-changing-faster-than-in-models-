"""

02/07/25

Plot land_trends in tail width calculated by Model_trends.py using Iris.

To run: python Model_trend_plotter.py [model] [experiment] [member]

"""

import iris
import sys
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import iris.plot as iplt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.cm as mpl_cm
import datetime
import matplotlib.dates as mdates
import numpy as np
import xarray as xr
from Lat_Land_mask_function import lat_land_mask
from Save_figure_function import save_figure
from Map_features import map_features

# Define colour map
cmap_rb = mpl_cm.get_cmap('coolwarm', 15)

# Command-line arguments
Mode = sys.argv[1]  # 'Coupled' or 'SST'
model = sys.argv[2]
experiment = sys.argv[3]
member = sys.argv[4]
sig = sys.argv[5]  # 'fdr' or 'p'

# Load trends data
datadir = f'/data/users/samuel.sancar/Heat_extremes/Trends/{Mode}_model_trends/'
filename = f"{model}_{experiment}_{member}_tailwidth_trend_1958-2022.nc"
cube = iris.load_cube(datadir + filename)
p87p5_trend = iris.load_cube(datadir + filename.replace('tailwidth','p87p5'))

# Apply land mask
masked_cube, _ = lat_land_mask(cube, -60, 90)
p87p5_masked, _ = lat_land_mask(p87p5_trend, -60, 90)

# Extract grid point locations for stippling
lon = cube.coord('longitude').points
lat = cube.coord('latitude').points
lon2d, lat2d = np.meshgrid(lon, lat)

# Find indices where p87p5_trend is negative
neg_indices = np.where(p87p5_masked.data < 0)
neg_lats = lat2d[neg_indices]
neg_lons = lon2d[neg_indices]

# Load significance mask (already masked for land and latitude)
if sig == 'fdr':
    sig_datadir = f'/data/users/samuel.sancar/Heat_extremes/Trends/Trend_significance/FDR_corrections_All/'
    statsig = iris.load_cube(f"{sig_datadir}{model}_{experiment}_{member}_fdr_significant_all.nc")
    statsig_bool =  statsig.data.astype(bool)
    stipple_lats, stipple_lons = lat2d[statsig_bool], lon2d[statsig_bool]
else:
    sig_datadir = f'/data/users/samuel.sancar/Heat_extremes/Trends/Trend_significance/P-vals/'
    statsig_ds = xr.open_dataset(f"{sig_datadir}{model}_{experiment}_{member}_tailwidths_p_1958.nc")
    statsig = list(statsig_ds.data_vars.values())[0]
    sig_mask_np = statsig.values < 0.05
    stipple_lats, stipple_lons = lat2d[sig_mask_np], lon2d[sig_mask_np]

# Create figure and axis
fig = plt.figure(figsize=(12, 6))
ax = plt.axes(projection=ccrs.Robinson(central_longitude=0))

iplt.pcolormesh(masked_cube, cmap=cmap_rb, vmin=-0.3, vmax=0.3)

# Overlay white dots where p87p5_trend is negative
ax.scatter(neg_lons, neg_lats, color='white', s=0.9, alpha=0.6, transform=ccrs.PlateCarree(), label='p87p5 trend < 0')

# Plot stippling
ax.scatter(stipple_lons, stipple_lats, color='black', s=0.5, alpha=0.6, transform=ccrs.PlateCarree(central_longitude=0), label='Stippling')

map_features('Rob', ax, -60, 90)

#cbar = plt.colorbar(ax.collections[0], orientation='horizontal', pad=0.05, aspect=50)
#cbar.ax.tick_params(labelsize=12)
#cbar.set_label('Trend in Tailwidth (°C per decade)', fontsize=14)

if experiment=='hist-1950':
    exp = 'Coupled'
else:
    exp = 'SST forced'
    
plt.title(f"{model} {exp} {member}", fontsize=17, weight='bold')

plt.tight_layout()
save_figure(plt, f'{model}_{experiment}_{member}_tailwidths_{sig}', 'Maps')
plt.show()

# Load full tail width data with time dimension
datadir = f'/data/users/samuel.sancar/Heat_extremes/Tail_widths/'
full_data_filename = f"{model}_{experiment}_{member}_tailwidths_1958-2022.nc"
full_cube = iris.load_cube(datadir + full_data_filename)

# Choose a grid point
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

plt.title(f"Tail Width timeseries for \n{model} {experiment} {member}\n at Grid Point {lat_index}, {lon_index} (London)", fontsize=14)
plt.xlabel("Time")
plt.ylabel("Tail Width (°C)")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()
print("Trend value:", trend)








