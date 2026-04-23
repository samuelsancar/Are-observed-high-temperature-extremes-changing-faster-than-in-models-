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
#era5_and_canari_mask.data = canari_stat_sig.data.astype(bool)
era5_canari_lat, era5_canari_lon = lat2d[era5_and_canari_mask], lon2d[era5_and_canari_mask]
canari_mask = canari_stat_sig & ~era5_mask
canari_lat, canari_lon = lat2d[canari_mask], lon2d[canari_mask]
era5_not_canari_mask = era5_mask & ~canari_stat_sig
era5_not_lat, era5_not_lon = lat2d[era5_not_canari_mask], lon2d[era5_not_canari_mask]

era5_only_mask = era5_mask & ~jra_mask
both_mask = era5_mask & jra_mask
jra_only_mask = jra_mask & ~era5_mask

era5_only_lat, era5_only_lon = lat2d[era5_only_mask], lon2d[era5_only_mask]
jra_only_lat, jra_only_lon = lat2d[jra_only_mask], lon2d[jra_only_mask]
both_lat, both_lon = lat2d[both_mask], lon2d[both_mask]

print(f"ERA5 and JRA-55 agree on only {len(both_lat)} grid points.")

import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from matplotlib.colors import ListedColormap, BoundaryNorm

# --- Ensure your boolean masks exist as in your script ---
# era5_only_mask, jra_only_mask, both_mask should be boolean arrays with shape (lat, lon)
# They should align with template grid (lat2d, lon2d).

# 1) Build categorical grid
cat = np.zeros_like(lat2d, dtype=np.uint8)
cat[era5_only_mask] = 1
cat[jra_only_mask]  = 2
cat[both_mask]      = 3

# Optional: mask out regions you don’t want to render (e.g., ocean or outside lat range)
# If your masks were already land-only within lat [-60, 90], you can skip this.
# Example:
# ocean_or_outside = ~land_mask | (lat2d < -60) | (lat2d > 90)
# cat = np.ma.array(cat, mask=ocean_or_outside)

# 2) Build lon/lat *edges* for pcolormesh (assuming regular 1° spacing)
lon_centers = template.coord('longitude').points  # e.g., [-179.5 ... 179.5] or [0..359] or [-180..180]
lat_centers = template.coord('latitude').points   # e.g., [-89.5 ... 89.5] or [-90..90]

def make_edges(centers):
    # Works for monotonic, nearly uniform spacing
    d = np.diff(centers)
    d = np.r_[d[0], d]  # assume first spacing for the lower edge
    edges = centers - d/2
    # Add one last edge
    edges = np.r_[edges, centers[-1] + d[-1]/2]
    return edges

lon_edges = make_edges(lon_centers)
lat_edges = make_edges(lat_centers)

# If your longitude wraps (e.g., 0..360), pcolormesh + PlateCarree can still work,
# but many prefer centering at 0°; use cartopy’s transform + set_extent as below.

# 3) Colormap & normalization
cmap = ListedColormap([
    (1,1,1,0),   # 0: transparent/white (no signal) -> alpha 0 helps see background
    (0,0,0,1),   # 1: black for ERA5-only
    (0,0,1,1),   # 2: blue  for JRA-only
    (1,0,0,1),   # 3: red   for Both
])
bounds = [ -0.5, 0.5, 1.5, 2.5, 3.5 ]  # bins centered on integers
norm = BoundaryNorm(bounds, cmap.N)

# 4) Plot
fig = plt.figure(figsize=(14, 12))
ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=0))
ax.set_extent([-140, 160, -55, 85], crs=ccrs.PlateCarree())

# Background
ax.add_feature(cfeature.LAND, facecolor='white', zorder=0)
ax.add_feature(cfeature.OCEAN, facecolor='white', zorder=0)
ax.coastlines(linewidth=0.7, zorder=2)

# pcolormesh (rasterized helps PDF/AI size; edgecolors off for speed)
mesh = ax.pcolormesh(
    lon_edges, lat_edges, cat,
    cmap=cmap, norm=norm,
    transform=ccrs.PlateCarree(),
    shading='flat', rasterized=True,
    hatch = 'O'
)

# Legend proxy handles (manual since categorical)
import matplotlib.patches as mpatches
handles = [
    mpatches.Patch(facecolor=(1,0,0,1), label='Both (ERA5 & JRA-55)'),
    mpatches.Patch(facecolor=(0,0,0,1), label='ERA5 only'),
    mpatches.Patch(facecolor=(0,0,1,1), label='JRA-55 only'),
]
ax.legend(handles=handles, loc='lower left')

plt.title(f"ERA5 and JRA-55 significant trends ({t1}-{t2-1}) [{sig}]", fontsize=16, weight='bold')
plt.tight_layout()
save_figure(plt, f'ERA5_JRA-55_sig_trends_pcolormesh_({t1}-{t2-1}({sig}))', 'Maps')
plt.show()