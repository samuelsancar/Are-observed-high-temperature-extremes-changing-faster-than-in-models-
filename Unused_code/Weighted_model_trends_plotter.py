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

cmap_rb = mpl_cm.get_cmap('coolwarm', 15)

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

# Load p-values for weighting
sig_datadir = f'/data/users/samuel.sancar/Heat_extremes/Trends/Trend_significance/P-vals/'
statsig_ds = xr.open_dataset(f"{sig_datadir}{model}_{experiment}_{member}_tailwidths_p_1958.nc")
statsig = list(statsig_ds.data_vars.values())[0]

# Compute Weighted Significance-Magnitude Index
trend_data = masked_cube.data 
p_values = statsig.values     

# Avoid log(0) by clipping p-values
p_values = np.clip(p_values, 1e-10, 1)
significance_score = -np.log10(p_values)

weighted_index = np.abs(trend_data) * significance_score

weighted_cube = masked_cube.copy()
weighted_cube.data = weighted_index

# Plot Weighted Significance-Magnitude Index
fig2 = plt.figure(figsize=(12, 6))
ax2 = plt.axes(projection=ccrs.Robinson(central_longitude=0))
ax2.set_extent([-180, 180, -60, 90], crs=ccrs.PlateCarree())

iplt.pcolormesh(weighted_cube, cmap='coolwarm', vmin=0, vmax=np.nanpercentile(weighted_index, 99))

# Add stippling for significant points (p < 0.05)
sig_mask_np = p_values < 0.05
stipple_lats, stipple_lons = lat2d[sig_mask_np], lon2d[sig_mask_np]
ax2.scatter(stipple_lons, stipple_lats, color='black', s=0.4, alpha=0.5, transform=ccrs.PlateCarree(), label='Significant')

ax2.coastlines(resolution='50m', linewidth=0.5)
ax2.add_feature(cfeature.LAND, facecolor='lightgray')
ax2.add_feature(cfeature.OCEAN, facecolor='white')

cbar2 = plt.colorbar(ax2.collections[0], orientation='horizontal', pad=0.05, aspect=50)
cbar2.set_label('Weighted Significance-Magnitude Index', fontsize=12)

plt.title(f"Weighted Significance-Magnitude Index\n{model} {experiment} {member}", fontsize=14, weight='bold')
plt.tight_layout()
save_figure(plt, f'{model}_{experiment}_{member}_weighted_index', 'Maps')
plt.show()