"""
Code to plot where ERA5 and MERRA significant trends agree.
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
from Lat_Land_mask_function import lat_land_mask
from Save_figure_function import save_figure
from ERA5_MKtest import significance_plot
from Map_features import map_features

cmap_rb = mpl_cm.get_cmap('coolwarm', 12)

t1 = 1980
t2 = 2023
sig = sys.argv[1] # either 'fdr' or 'pvals'
alpha = 0.05 # For MK_test
lat_minimum = -60
lat_maximum = 90

datadir = '/data/users/samuel.sancar/Heat_extremes/Trends/Reanalysis_trends/'
MERRA_datadir = '/data/users/samuel.sancar/Heat_extremes/Trends/Reanalysis_trends/'
era5_name = 'era5_daily_maximum_2m_temperature'
merra_name = 'MERRA_daily_maximum_2m_temperature'

template = iris.load_cube('template_cube_1deg.nc')

# Extract grid point locations for stippling
lon = template.coord('longitude').points
lat = template.coord('latitude').points
lon2d, lat2d = np.meshgrid(lon, lat)

# Load significance mask (already masked for land and latitude)
def significance(datadir, variable, sig):
    
    if sig == 'fdr':
        statsig = iris.load_cube(f"{datadir}{variable}_tailwidth_fdr_all_{t1}-{t2}.nc")
        statsig_bool =  statsig.data.astype(bool)
        stipple_lats, stipple_lons = lat2d[statsig_bool], lon2d[statsig_bool]
        
    else:
        statsig_ds = xr.open_dataset(f"{datadir}{variable}_tailwidth_p_{t1}-{t2}.nc")
        statsig = list(statsig_ds.data_vars.values())[0]
        statsig_bool = statsig < alpha
        
    stipple_lats, stipple_lons = lat2d[statsig_bool], lon2d[statsig_bool]
        
    return statsig_bool, stipple_lats, stipple_lons
    
era5_mask, era5_lat, era5_lon = significance(datadir, era5_name, sig)
merra_mask, merra_lat, merra_lon = significance(MERRA_datadir, merra_name, sig)

era5_only_mask = era5_mask & ~merra_mask
both_mask = era5_mask & merra_mask
merra_only_mask = merra_mask & ~era5_mask

era5_only_lat, era5_only_lon = lat2d[era5_only_mask], lon2d[era5_only_mask]
merra_only_lat, merra_only_lon = lat2d[merra_only_mask], lon2d[merra_only_mask]
both_lat, both_lon = lat2d[both_mask], lon2d[both_mask]

# Plot trends
fig = plt.figure(figsize=(12,12))
ax =plt.axes(projection = ccrs.PlateCarree(central_longitude=0))

# Background
ax.add_feature(cfeature.LAND, facecolor='white', zorder=0)
ax.add_feature(cfeature.OCEAN, facecolor='white', zorder=0)
ax.coastlines()

    # Plot significant trends (data are in geographic coords)
ax.scatter(
    both_lon, both_lat, color='red', marker='o', s=6, alpha=1,
    transform=ccrs.PlateCarree(), label='Both'
)

ax.scatter(
    era5_only_lon, era5_only_lat, color='black', marker='o', s=6, alpha=1,
    transform=ccrs.PlateCarree(), label='ERA5'
)
ax.scatter(
    merra_only_lon, merra_only_lat, color='blue', marker='o', s=6, alpha=1,
    transform=ccrs.PlateCarree(), label='MERRA'
)

# Apply region extent: [west, east, south, north]

ax.set_extent([-7.5, 10.5, 42.5, 62.5], crs=ccrs.PlateCarree())

# Figure-level title
plt.title(f"ERA5 and MERRA NW Europe significant trends ({t1}-{t2-1})", fontsize=18, weight='bold')

# Single, figure-level legend (avoid per-axes duplicates)
# Build proxy artists to ensure order is clear
import matplotlib.patches as mpatches
handles = [
    mpatches.Patch(facecolor='red',   label='Both'),
    mpatches.Patch(facecolor='black', label='ERA5'),
    mpatches.Patch(facecolor='blue',  label='merra-55'),
]
fig.legend(handles=handles, loc='lower center', ncol=3)

#save_figure(plt, f'ERA5_MERRA_sig_trends_panels_({t1}-{t2-1}({sig}))', 'Maps')
plt.show()