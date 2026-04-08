"""
17/07/25
---------
Plot trend in tail width calculated by ERA5_tailwidth.py
To run: python plot_tailwidth_trend.py [1958] [2023] [fdr|p]
"""

import iris
import sys
import matplotlib.pyplot as plt
import iris.plot as iplt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.cm as mpl_cm
import numpy as np
import xarray as xr
from matplotlib.colors import ListedColormap
from Lat_Land_mask_function import lat_land_mask
from Save_figure_function import save_figure
from ERA5_MKtest import significance_plot
from Map_features import map_features

reds = mpl_cm.get_cmap('Reds')

t1 = 1958
t2 = 2023
sig = sys.argv[1] # either 'fdr' or 'pvals'
alpha = 0.05 # For MK_test
lat_minimum = 60
lat_maximum = 90

reanalysis = 'era5'
variable = 'daily_maximum_2m_temperature'
datadir = '/data/users/samuel.sancar/Heat_extremes/Trends/Reanalysis_trends/'
trend = iris.load_cube(f"{datadir}{reanalysis}_{variable}_tailwidth_trends_{t1}-{t2}.nc")
p87p5_trend = iris.load_cube(f"{datadir}{reanalysis}_{variable}_p87p5_trends_{t1}-{t2}.nc")
template = iris.load_cube('template_cube_1deg.nc')

# Apply latitude and land mask using function
model_mask, _ = lat_land_mask(trend, lat_minimum, lat_maximum)
p87p5_masked, _ = lat_land_mask(p87p5_trend, lat_minimum, lat_maximum)

# Extract grid point locations for stippling
lon = template.coord('longitude').points
lat = template.coord('latitude').points
lon2d, lat2d = np.meshgrid(lon, lat)

# Load significance mask (already masked for land and latitude)
if sig == 'fdr':
    statsig = iris.load_cube(f"{datadir}{reanalysis}_{variable}_tailwidth_fdr_all_{t1}-{t2}.nc")
    statsig_bool =  statsig.data.astype(bool)
    stipple_lats, stipple_lons = lat2d[statsig_bool], lon2d[statsig_bool]
    
else:
    statsig_ds = xr.open_dataset(f"{datadir}{reanalysis}_{variable}_tailwidth_p_{t1}-{t2}.nc")
    statsig = list(statsig_ds.data_vars.values())[0]
    sig_mask_np = statsig < alpha
    stipple_lats, stipple_lons = lat2d[sig_mask_np], lon2d[sig_mask_np]

neg_indices = np.where(p87p5_masked.data < 0)
neg_lats = lat2d[neg_indices]
neg_lons = lon2d[neg_indices]


strong = model_mask.data.copy()
strong[strong <= 0.6] = np.nan   # everything ≤ 0.6 becomes white (transparent)
strong_cube = model_mask.copy(data=strong)


from matplotlib.colors import LinearSegmentedColormap

def truncate_cmap(cmap, vmin=0.2, vmax=1.0):
    new_colors = cmap(np.linspace(vmin, vmax, 256))
    return LinearSegmentedColormap.from_list("trunc", new_colors)

red_trunc = truncate_cmap(reds,0.3,1.0)

# Plot trends
fig = plt.figure(figsize=(12, 6))
ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=0))

map = iplt.pcolormesh(strong_cube, cmap=red_trunc, vmin=0.6, vmax=1)

ax.set_facecolor("white")

map_features('Plate', ax, 60, 90)

cbar = plt.colorbar(map, ax=ax, orientation='horizontal', pad=0.05, aspect=50)
cbar.ax.tick_params(labelsize=14)
cbar.set_label('Trend in Tailwidth (°C per decade)', fontsize=15)

plt.title(f"ERA5 strongest tail width trends ({t1}-{t2-1})", fontsize=18, weight='bold')

gl = ax.gridlines(draw_labels=True, linewidth=0.5, color='gray', alpha=0.5, linestyle='--')

gl.top_labels = False
gl.right_labels = False

plt.tight_layout()
save_figure(plt, f'ERA5_taiwidth_trends_({t1}-{t2-1}({sig}))', 'Maps')
plt.show()