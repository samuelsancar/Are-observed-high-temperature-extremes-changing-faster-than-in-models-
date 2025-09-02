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
sys.path.append('/home/users/samuel.sancar/Heat_Extremes/Code/Functions/')
from Lat_Land_mask_function import lat_land_mask

# Define colour map
cmap_rb = mpl_cm.get_cmap('RdBu_r', 15)

# Command-line arguments
model = sys.argv[1]
experiment = sys.argv[2]
member = sys.argv[3]

# Load trends data
datadir = '/data/users/samuel.sancar/Heat_extremes/trends/Coupled_model_trends/'
filename = f"{model}_{experiment}_{member}_tailwidth_trend.nc"
cube = iris.load_cube(datadir + filename)

# Apply land mask
masked_cube = lat_land_mask(cube, -60)

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



