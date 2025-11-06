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
from Save_figure_function import save_figure

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
filename = f"{model}_{experiment}_{member}_tailwidth_trend.nc"
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
    sig_datadir = f'/data/users/samuel.sancar/Heat_extremes/Trends/Trend_significance/P-vals/'
    statsig = xr.DataArray(
        data=np.load(f"{sig_datadir}{model}_{experiment}_{member}_tailwidth_fdr_1958.npy"),
        dims=["latitude", "longitude"],
        coords={"latitude": lat, "longitude": lon},
        name="significance_mask"
    )
    stipple_lats, stipple_lons = lat2d[statsig], lon2d[statsig]
else:
    sig_datadir = f'/data/users/samuel.sancar/Heat_extremes/Trends/Trend_significance/P-vals/'
    statsig_ds = xr.open_dataset(f"{sig_datadir}{model}_{experiment}_{member}_tailwidths_p_1958.nc")
    statsig = list(statsig_ds.data_vars.values())[0]
    sig_mask_np = statsig.values < 0.05
    stipple_lats, stipple_lons = lat2d[sig_mask_np], lon2d[sig_mask_np]

# Create figure and axis
fig = plt.figure(figsize=(12, 6))
ax = plt.axes(projection=ccrs.Robinson(central_longitude=0))
ax.set_extent([-180, 180, -60, 90], crs=ccrs.PlateCarree())

# Plot using Iris
iplt.pcolormesh(masked_cube, cmap=cmap_rb, vmin=-0.3, vmax=0.3)

# Overlay white dots where p87p5_trend is negative
ax.scatter(neg_lons, neg_lats, color='white', s=0.9, alpha=0.6, transform=ccrs.PlateCarree(), label='p87p5 trend < 0')

# Add map features
ax.coastlines(resolution='50m', linewidth=0.5)
ax.add_feature(cfeature.NaturalEarthFeature(category='cultural', name='admin_0_boundary_lines_land', facecolor='none', scale='50m', edgecolor='0', lw=.05))
ax.add_feature(cfeature.LAKES.with_scale('110m'), edgecolor='0', facecolor='none', lw=.05)
ax.add_feature(cfeature.STATES.with_scale('110m'), edgecolor='0', facecolor='none', lw=.2)
ax.add_feature(cfeature.RIVERS.with_scale('110m'), edgecolor='tab:blue', facecolor='none', lw=.4)
ax.add_feature(cfeature.LAND, facecolor='lightgray')
ax.add_feature(cfeature.OCEAN, facecolor='white')

# Plot stippling
ax.scatter(stipple_lons, stipple_lats, color='black', s=0.4, alpha=0.5, transform=ccrs.PlateCarree(central_longitude=0), label='Stippling')

# Add gridlines
gl = ax.gridlines(draw_labels=True, linewidth=0.5, color='gray', alpha=0.5, linestyle='--')
gl.top_labels = False
gl.right_labels = False

# Add colourbar
cbar = plt.colorbar(ax.collections[0], orientation='horizontal', pad=0.05, aspect=50)
cbar.set_label('Trend in Tailwidth (°C per decade)', fontsize=12)

# Add title
plt.title(f"Trend in Tailwidth for\n{model} {experiment} {member} (1958-2022)", fontsize=14, weight='bold')

# Layout
plt.tight_layout()
save_figure(plt, f'{model}_{experiment}_{member}_tailwidths_p', 'Maps')
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

#if not np.isnan(trend):
    #plt.plot(converted_time_points, trend_line, color='red', linestyle='--', label='Trend Line')
#else:
    #print("Trend is NaN at the selected grid point")

plt.title(f"Tail Width Over Time\n{model} {experiment} {member} at Grid Point ({lat_index}, {lon_index})", fontsize=14)
plt.xlabel("Time")
plt.ylabel("Tail Width (°C)")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()
print("Trend value:", trend)








