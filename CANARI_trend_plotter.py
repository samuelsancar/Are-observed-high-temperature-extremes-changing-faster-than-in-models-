"""

04/11/25

Plot land_trends in tail width calculated by CANARI_trends.py.
Note: CANARI significant trends are calculated in CANARI_MKtest.py and saved as netCDF files. These must be calculated first for code to work. 
These are loaded here to add stippling to the trend map.

To run: 'python CANARI_trend_plotter.py [member] [fdr|p]' where fdr or p determines whether to use the 
FDR-corrected significance mask or the raw p-value mask for stippling.

"""

import iris
import sys
import matplotlib.pyplot as plt
import iris.plot as iplt
import cartopy.crs as ccrs
import matplotlib.cm as mpl_cm
import datetime
import numpy as np
import xarray as xr
from Lat_Land_mask_function import lat_land_mask
from Map_features import map_features
from Save_figure_function import save_figure

# Define colour map
cmap_rb = mpl_cm.get_cmap('coolwarm', 15)

t1 = 1958
t2 = 2023

alpha =  0.05

member = sys.argv[1]
sig = sys.argv[2]

# Latitude limits for plotting
lat_minimum = -60
lat_maximum = 90

# Load trends data
template = iris.load_cube('template_cube_1deg.nc')
filename = f"{member}_tailwidth_trend_{t1}-{t2-1}.nc"
datadir = f'/data/users/samuel.sancar/Heat_extremes/Trends/CANARI_trends/'
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
lon = template.coord('longitude').points
lat = template.coord('latitude').points
lon2d, lat2d = np.meshgrid(lon, lat)

# Load significance mask (already masked for land and latitude)
if sig == 'fdr':
    sig_datadir = f'/data/users/samuel.sancar/Heat_extremes/Trends/Trend_significance/CANARI/FDR_corrections/'
    statsig = iris.load_cube(f"{sig_datadir}{member}_fdr_corrected_p_1958-2022.nc")
    data = statsig.data  # masked array
    stipple_mask = (~data.mask) & data.astype(bool)
    stipple_lats, stipple_lons = lat2d[stipple_mask], lon2d[stipple_mask]
else:
    sig_datadir = f'/data/users/samuel.sancar/Heat_extremes/Trends/Trend_significance/CANARI/P-vals/'
    statsig = iris.load_cube(f"{sig_datadir}{member}_tailwidth_p_{t1}-{t2-1}.nc")
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

# Plot stippling
ax.scatter(stipple_lons, stipple_lats, color='black', s=0.4, alpha=0.5, transform=ccrs.PlateCarree(central_longitude=0), label='Stippling')

# Add colourbar
cbar = plt.colorbar(map, ax=ax, orientation='horizontal', pad=0.05, aspect=50)
cbar.ax.tick_params(labelsize=16)
cbar.set_label('Trend in Tailwidth (°C per decade)', fontsize=12)

plt.title(f"Trend in Tailwidth for\n CANARI member {member} ({t1}-t2-1)", fontsize=14, weight='bold')
plt.tight_layout()
save_figure(plt, f'CANARI_{member}_tailwidths_p', 'Maps')
plt.show()

# Load full tail width data with time dimension
datadir = f'/data/users/samuel.sancar/Heat_extremes/Tail_widths/canari/'
full_data_filename = f"{member}_tailwidths_{t1}-{t2-1}.nc "
full_cube = iris.load_cube(datadir + full_data_filename)

# Choose a grid point (default is over London)
lat_index = 55
lon_index =37

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
time_coord = full_cube.coord('time')
time_points = time_coord.units.num2date(time_coord.points)

# Convert cftime.DatetimeProlepticGregorian to datetime.datetime
converted_time_points = [datetime.datetime(tp.year, tp.month, tp.day) for tp in time_points]

# Extract trend value from masked_cube
trend = cube.data[lon_box, lon_box]

# Generate trend line using the same x-axis scale
trend_line = trend/10 * np.arange(len(time_series)) + time_series[0]

# Plot the time series and trend line
plt.figure(figsize=(10, 5))
plt.plot(converted_time_points, time_series, marker='o', linestyle='-', label='Tail Width')

plt.title(f"Tail Width Over Time for CANARI member {member} at Grid Point ({lat_index}, {lon_index})", fontsize=14)
plt.xlabel("Time")
plt.ylabel("Tail Width (°C)")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()
print("Trend value:", trend)