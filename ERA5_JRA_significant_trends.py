"""
Code to plot where ERA5 and JRA-55 significant trends agree.
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
from Save_figure_function import save_figure
from Map_features import map_features

cmap_rb = mpl_cm.get_cmap('coolwarm', 12)

t1 = 1958
t2 = 2023
sig = sys.argv[1] # either 'fdr' or 'pvals'
alpha = 0.05 # For MK_test
lat_minimum = -60
lat_maximum = 90

datadir = '/data/users/samuel.sancar/Heat_extremes/Trends/Reanalysis_trends/'
era5_name = 'era5_daily_maximum_2m_temperature'
era5_p87p5 = iris.load_cube(f"{datadir}{era5_name}_p87p5_trends_{t1}-{t2}.nc")
jra_name = 'jra55__minmax_surf.015_tmax.reg_tl319.'
jra_p87p5 = iris.load_cube(f"{datadir}{jra_name}_p87p5_trends_{t1}-{t2}.nc")
canari_datadir = '/data/users/samuel.sancar/Heat_extremes/Trends/Trend_significance/CANARI/'

members= np.arange(1,41)

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
jra_mask, jra_lat, jra_lon = significance(datadir, jra_name, sig)

era5_only_mask = era5_mask & ~jra_mask
both_mask = era5_mask & jra_mask
jra_only_mask = jra_mask & ~era5_mask

era5_only_lat, era5_only_lon = lat2d[era5_only_mask], lon2d[era5_only_mask]
jra_only_lat, jra_only_lon = lat2d[jra_only_mask], lon2d[jra_only_mask]
both_lat, both_lon = lat2d[both_mask], lon2d[both_mask]

print(f"ERA5 and JRA-55 agree on only {len(both_lat)} grid points.")

# Define regions as [west, east, south, north]
#regions,  = (NW Europe, [-7.5, 12.5, 42.5, 62.5]),     # NW Europe (square)
    #('Maritime Continent', [112.5, 132.5, -10, 10]),       # Maritime Continent (square)
    #('Argentina', [-74, -54, -47.5, -27.5]),      # Argentina (square)
    #('Persian Gulf', [43, 63, 17, 37]),              # Persian Gulf (square)


fig = plt.figure(figsize=(12,12))
ax =plt.axes(projection = ccrs.PlateCarree(central_longitude=0))

# Background
ax.add_feature(cfeature.LAND, facecolor='white', zorder=0)
ax.add_feature(cfeature.OCEAN, facecolor='white', zorder=0)
ax.coastlines()

    # Plot significant trends (data are in geographic coords)
ax.scatter(
    both_lon, both_lat, color="#ff1e22", marker='o', s=10, alpha=1,
    transform=ccrs.PlateCarree(), label='Both'
)

ax.scatter(
    era5_only_lon, era5_only_lat, marker='s', s=4, alpha=0.6, 
    facecolors='black',
    edgecolors='white',       
    linewidths=0.5,
    transform=ccrs.PlateCarree(), label='ERA5'
)

ax.scatter(
    jra_only_lon, jra_only_lat, marker='s', s=4, alpha=0.5,
    facecolors='blue',
    edgecolors='white',       
    linewidths=0.5,
    transform=ccrs.PlateCarree(), label='JRA-55'
)

# Apply region extent: [west, east, south, north]

#ax.set_extent([-7.5, 10.5, 42.5, 62.5], crs=ccrs.PlateCarree())
ax.set_extent([-135, 160, -60, 90], crs=ccrs.PlateCarree())

# Figure-level title
plt.title(f"ERA5 and JRA-55 significant trends ({t1}-{t2-1})", fontsize=18, weight='bold')

# Single, figure-level legend (avoid per-axes duplicates)
# Build proxy artists to ensure order is clear
import matplotlib.patches as mpatches
handles = [
    mpatches.Patch(facecolor='Red',   label='Both'),
    mpatches.Patch(facecolor='black', label='ERA5'),
    mpatches.Patch(facecolor='blue',  label='JRA-55'),
]
fig.legend(handles=handles, loc='lower center', ncol=3)

save_figure(plt, f'ERA5_JRA-55_sig_trends_panels_({t1}-{t2-1}({sig}))', 'Maps')
plt.show()


