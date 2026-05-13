"""
Code to plot where ERA5 and MERRA significant trends agree.
"""

import iris
import sys
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.cm as mpl_cm
import numpy as np
from Lat_Land_mask_function import lat_land_mask
from Save_figure_function import save_figure

cmap_rb = mpl_cm.get_cmap('coolwarm', 12)

# Beginning & end year+1
t1 = 1980
t2 = 2023

sig = sys.argv[1] # either 'fdr' or 'pvals'
alpha = 0.05 # For MK_test

lat_minimum = -60
lat_maximum = 80

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
        statsig = iris.load_cube(f"{datadir}{variable}_fdr_corrected_p_{t1}-{t2-1}.nc")
        data = statsig.data  # masked array
        sig_bool = (~data.mask) & data.astype(bool)
        stipple_lats, stipple_lons = lat2d[sig_bool], lon2d[sig_bool]
    else:
        statsig = iris.load_cube(f"{datadir}{variable}_tailwidth_p_{t1}-{t2-1}.nc")
        statsig_masked, _ = lat_land_mask(statsig, template, lat_minimum, lat_maximum)
        sig_bool = np.ma.filled(statsig_masked.data < alpha, False)
        stipple_lats, stipple_lons = lat2d[sig_bool], lon2d[sig_bool]
        
    return sig_bool, stipple_lats, stipple_lons
    
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

# Plot significant trends 
ax.scatter(
    both_lon, both_lat, color="#ff1e22", marker='o', s=10, alpha=1,transform=ccrs.PlateCarree(), label='Both')

ax.scatter(
    era5_only_lon, era5_only_lat,  marker='s', s=4, alpha=0.6, facecolors='black',edgecolors='white', linewidths=0.5,transform=ccrs.PlateCarree(), label='ERA5')

ax.scatter(
    merra_only_lon, merra_only_lat, marker='s', s=4, alpha=0.5,facecolors='blue',edgecolors='white', linewidths=0.5,transform=ccrs.PlateCarree(), label='MERRA2')

ax.set_extent([-135, 160, -60, 90], crs=ccrs.PlateCarree())

plt.title(f"b) ERA5 and MERRA-2 significant trends ({t1}-{t2-1})", fontsize=18, weight='bold')

import matplotlib.patches as mpatches
handles = [
    mpatches.Patch(facecolor='red',   label='Both'),
    mpatches.Patch(facecolor='black', label='ERA5'),
    mpatches.Patch(facecolor='blue',  label='MERRA-2'),
]
fig.legend(handles=handles, loc='lower center', ncol=3)

save_figure(plt, f'ERA5_MERRA-2_sig_trends({t1}-{t2-1}({sig}))', 'Maps')
plt.show()