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
sys.path.append('/home/users/samuel.sancar/Heat_Extremes/Code/Functions/')
from Lat_Land_mask_function import lat_land_mask

# Define colour map
cmap_rb = mpl_cm.get_cmap('RdBu_r', 15)

# Command-line arguments
Mode = sys.argv[1]  # 'Coupled' or 'SST'
model = sys.argv[2]
experiment = sys.argv[3]
member = sys.argv[4]

# Load trends data
datadir = f'/data/users/samuel.sancar/Heat_extremes/Trends/{Mode}_model_trends/'
filename = f"{model}_{experiment}_{member}_tailwidth_trend.nc"
cube = iris.load_cube(datadir + filename)

# Apply land mask
masked_cube, _ = lat_land_mask(cube, -60)

# Plotting
fig = plt.figure(figsize=(12, 6))
ax = plt.axes(projection=ccrs.PlateCarree())

# Plot using Iris
mesh = iplt.pcolormesh(masked_cube, cmap=cmap_rb, vmin=-0.4, vmax=0.4)
ax.coastlines(resolution='50m', linewidth=1)
ax.add_feature(cfeature.BORDERS, linestyle=':')
ax.add_feature(cfeature.LAND, facecolor='lightgray')
ax.add_feature(cfeature.OCEAN, facecolor='white')

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

# Choose a grid point (e.g., latitude index 50, longitude index 100)
lat_index = 50
lon_index = 100

# Extract time series for the selected grid point
time_series = full_cube.data[:, lat_index, lon_index]
time_coord = full_cube.coord('time')
time_points = time_coord.units.num2date(time_coord.points)


# Convert cftime.DatetimeProlepticGregorian to datetime.datetime
converted_time_points = [datetime.datetime(tp.year, tp.month, tp.day) for tp in time_points]


# Plot the time series
plt.figure(figsize=(10, 5))
plt.plot(time_points, time_series, marker='o', linestyle='-')
plt.title(f"Tail Width Over Time\n{model} {experiment} {member} at Grid Point ({lat_index}, {lon_index})", fontsize=14)
plt.xlabel("Time")
plt.ylabel("Tail Width (°C)")
plt.grid(True)
plt.tight_layout()
plt.show()




