"""
Docstring for Untitled-1
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

def canari_sig_trends(i, sig):
    
    if sig == 'fdr':
        statsig = iris.load_cube(f"{canari_datadir}/FDR_corrections_All/{i}_tailwidths_p_1958_fdr_significant_all.nc")
        statsig_bool =  statsig.data.astype(bool)
        stipple_lats, stipple_lons = lat2d[statsig_bool], lon2d[statsig_bool]
        
    else:
        statsig_ds = xr.open_dataset(f"{canari_datadir}/P-vals/{i}_tailwidths_p_1958.nc")
        statsig = list(statsig_ds.data_vars.values())[0]
        statsig_bool = statsig < alpha
        
    stipple_lats, stipple_lons = lat2d[statsig_bool], lon2d[statsig_bool]
        
    return statsig_bool, stipple_lats, stipple_lons


canari_stat_sig = np.logical_or.reduce(
    [np.asarray(canari_sig_trends(i, sig)[0]).astype(bool) for i in members]
)

era5_and_canari_mask = era5_mask & canari_stat_sig
era5_canari_lat, era5_canari_lon = lat2d[era5_and_canari_mask], lon2d[era5_and_canari_mask]
canari_mask = canari_stat_sig & ~era5_mask
canari_lat, canari_lon = lat2d[canari_mask], lon2d[canari_mask]
era5_not_canari_mask = era5_mask & ~canari_stat_sig
era5_not_lat, era5_not_lon = lat2d[era5_not_canari_mask], lon2d[era5_not_canari_mask]

# Plot trends
fig = plt.figure(figsize=(12, 6))
ax = plt.axes(projection=ccrs.Robinson(central_longitude=0))

ax.scatter(era5_canari_lon, era5_canari_lat, color='orange', s=1, alpha=0.8, transform=ccrs.PlateCarree(central_longitude=0), label='Both')
ax.scatter(canari_lon, canari_lat, color='green', s=1, alpha=0.8, transform=ccrs.PlateCarree(central_longitude=0), label='CANARI only')
ax.scatter(era5_not_lon, era5_not_lat, color = 'blue', s=1, alpha=0.8, transform=ccrs.PlateCarree(central_longitude=0), label='ERA5 only')

ax.set_extent([-180, 180, -60, 90], crs=ccrs.PlateCarree(central_longitude=0))
ax.coastlines(resolution='50m', linewidth=0.5)
ax.add_feature(cfeature.NaturalEarthFeature(category='cultural', name='admin_0_boundary_lines_land', facecolor='none', scale='50m', edgecolor='0', lw=.05))
ax.add_feature(cfeature.LAKES.with_scale('110m'), edgecolor='0', facecolor='none', lw=.05)
ax.add_feature(cfeature.STATES.with_scale('110m'), edgecolor='0', facecolor='none', lw=.2)
#ax.add_feature(cfeature.RIVERS.with_scale('110m'), edgecolor='tab:blue', facecolor='none', lw=.4)
ax.add_feature(cfeature.LAND, facecolor='white')
ax.add_feature(cfeature.OCEAN, facecolor='white')

plt.title(f"ERA5 and JRA-55 significant trends ({t1}-{t2-1})", fontsize=18, weight='bold')
plt.legend()
plt.tight_layout()
save_figure(plt, f'ERA5_JRA-55_sig_trends_({t1}-{t2-1}({sig}))', 'Maps')
plt.show()