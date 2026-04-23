"""
20/07/25
--------

Plot HighResMIP trend in tail widths calculated by HighResMIP_trend_calc.py.

Used to produce plots for figure 1.

Note: HighResMIP significant trends are calculated in HighResMIP_MKtest.py and saved as netCDF files. These must be calculated first for code to work. 
These are loaded here to add stippling to the trend map.

To run: 'python HighResMIP_trend_plotter.py [Mode] [model [experiment] [member] [sig] [fdr|p]' # where fdr or p determines whether to use the 
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

# Define colour map
cmap_rb = mpl_cm.get_cmap('coolwarm', 15)

t1 = 1958
t2= 2023

alpha = 0.05

lat_minimum = -60
lat_maximum = 90

# Command-line arguments
Mode = sys.argv[1]  # 'Coupled' or 'SST'
model = sys.argv[2]
experiment = sys.argv[3]
member = sys.argv[4]
sig = sys.argv[5]  # 'fdr' or 'p'
datadir = f'/data/users/samuel.sancar/Heat_extremes/Trends/{Mode}_model_trends/'

if experiment=='hist-1950':
    exp = 'Coupled'
else:
    exp = 'SST forced'

template = iris.load_cube('template_cube_1deg.nc')
filename = f"{model}_{experiment}_{member}_tailwidth_trend_1958-2022.nc"
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
    sig_datadir = f'/data/users/samuel.sancar/Heat_extremes/Trends/Trend_significance/FDR_corrections_All/'
    statsig = iris.load_cube(f"{sig_datadir}{model}_{experiment}_{member}_fdr_significant_all.nc")
    data = statsig.data  # masked array
    stipple_mask = (~data.mask) & data.astype(bool)
    stipple_lats, stipple_lons = lat2d[stipple_mask], lon2d[stipple_mask]
else:
    sig_datadir = f'/data/users/samuel.sancar/Heat_extremes/Trends/Trend_significance/P-vals/'
    statsig = iris.load_cube(f"{sig_datadir}{model}_{experiment}_{member}_tailwidth_p_{t1}-{t2-1}.nc")
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

map = iplt.pcolormesh(masked_cube, cmap=cmap_rb, vmin=-0.3, vmax=0.3)

# Overlay white dots where p87p5_trend is negative
ax.scatter(neg_lons, neg_lats, color='white', s=0.9, alpha=0.6, transform=ccrs.PlateCarree(), label='p87p5 trend < 0')

# Plot significant trends
ax.scatter(stipple_lons, stipple_lats, color='black', s=0.5, alpha=0.6, transform=ccrs.PlateCarree(central_longitude=0), label='Stippling')

cbar = plt.colorbar(map, ax = ax, orientation='horizontal', pad=0.05, aspect=50)
cbar.ax.tick_params(labelsize=12)
cbar.set_label('Trend in Tailwidth (°C per decade)', fontsize=14)
    
plt.title(f"{model} {exp} {member}", fontsize=17, weight='bold')
plt.tight_layout()
save_figure(plt, f'{model}_{experiment}_{member}_tailwidths_{sig}', 'Maps')
plt.show()

# Plot time series at a gridpoint (London) to check that trends make sense. 
# Extract from full data cube with time dimension, not the trend cube.

# Load full tail width data with time dimension
datadir = f'/data/users/samuel.sancar/Heat_extremes/Tail_widths/'
full_data_filename = f"{model}_{experiment}_{member}_tailwidths_1958-2022.nc"
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
plt.title(f"Tail Width timeseries for \n{model} {experiment} {member}\n at Grid Point {lat_index}, {lon_index} (London)", fontsize=14)
plt.xlabel("Time")
plt.ylabel("Tail Width (°C)")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()
print("Trend value:", trend)








