"""

02/07/25

Plot land_trends in tail width calculated by Model_trends.py using Iris.

To run: python Model_trend_plotter.py [model] [experiment] [member]

"""

import iris
import sys
import matplotlib.pyplot as plt
import iris.plot as iplt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.cm as mpl_cm
import datetime
import matplotlib.dates as mdates
import numpy as np
import xarray as xr
sys.path.append('/home/users/samuel.sancar/Heat_Extremes/Code/Functions/')
from Lat_Land_mask_function import lat_land_mask

# Define colour map
cmap_rb = mpl_cm.get_cmap('RdBu_r', 15)

# Command-line arguments
Mode = sys.argv[1]  # 'Coupled' or 'SST'
model = sys.argv[2]
experiment = sys.argv[3]
member = sys.argv[4]
sig = sys.argv[5]  # 'fdr' or 'p'

# Load trends data
datadir = f'/data/users/samuel.sancar/Heat_extremes/Trends/{Mode}_model_trends/'
filename = f"{model}_{experiment}_{member}_tailwidth_trend.nc"
cube = iris.load_cube(datadir + filename)

# Apply land mask
masked_cube, _ = lat_land_mask(cube, -60, 90)

# Extract grid point locations for stippling
lon = cube.coord('longitude').points
lat = cube.coord('latitude').points
lon2d, lat2d = np.meshgrid(lon, lat)

# Load significance mask (already masked for land and latitude)
#if sig == 'fdr':
#    statsig = xr.DataArray(
#        data=np.load(f"{datadir}{model}_{experiment}_{member}_tailwidth_fdr_1958-2023.npy"),
#        dims=["latitude", "longitude"],
#        coords={"latitude": lat, "longitude": lon},
#        name="significance_mask"
#    )
#    stipple_lats, stipple_lons = lat2d[statsig], lon2d[statsig]
#else:
sig_datadir = f'/data/users/samuel.sancar/Heat_extremes/Trends/Trend_significance/'
statsig_ds = xr.open_dataset(f"{sig_datadir}{model}_{experiment}_{member}_tailwidths_tailwidth_p_1958.nc")
statsig = list(statsig_ds.data_vars.values())[0]
sig_mask_np = statsig.values < 0.05
stipple_lats, stipple_lons = lat2d[sig_mask_np], lon2d[sig_mask_np]

# Plotting
fig = plt.figure(figsize=(12, 6))
ax = plt.axes(projection=ccrs.PlateCarree())

# Plot using Iris
mesh = iplt.pcolormesh(masked_cube, cmap=cmap_rb, vmin=-0.4, vmax=0.4)
ax.coastlines(resolution='50m', linewidth=1)
ax.add_feature(cfeature.BORDERS, linestyle=':')
ax.add_feature(cfeature.LAND, facecolor='lightgray')
ax.add_feature(cfeature.OCEAN, facecolor='white')

# Plot stippling
ax.scatter(stipple_lons, stipple_lats, color='black', s=0.5, alpha=0.5, transform=ccrs.PlateCarree(), label='Stippling')

gl = ax.gridlines(draw_labels=True, linewidth=0.5, color='gray', alpha=0.5, linestyle='--')
gl.top_labels = False
gl.right_labels = False

cbar = plt.colorbar(mesh, orientation='horizontal', pad=0.05, aspect=50)
cbar.set_label('Trend in Tail Width (°C per decade)', fontsize=12)

plt.title(f"Trend in Tail Width\n{model} {experiment} {member}", fontsize=14, weight='bold')
plt.tight_layout()
plt.show()

# Load full tail width data with time dimension
datadir = f'/data/users/samuel.sancar/Heat_extremes/Tail_widths/'
full_data_filename = f"{model}_{experiment}_{member}_tailwidths.nc"
full_cube = iris.load_cube(datadir + full_data_filename)

# Choose a grid point
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

if not np.isnan(trend):
    plt.plot(converted_time_points, trend_line, color='red', linestyle='--', label='Trend Line')
else:
    print("Trend is NaN at the selected grid point")

plt.title(f"Tail Width Over Time\n{model} {experiment} {member} at Grid Point ({lat_index}, {lon_index})", fontsize=14)
plt.xlabel("Time")
plt.ylabel("Tail Width (°C)")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()
print("Trend value:", trend)








