""" 
17/11/2025
----------
CANARI regional analysis

"""

import iris
import sys
import glob
import matplotlib.pyplot as plt
import iris.plot as iplt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.cm as mpl_cm
import numpy as np
import xarray as xr
from matplotlib.patches import Rectangle
import iris.analysis as ia
import iris.analysis.cartography as iac
from iris.coord_categorisation import add_year
import warnings
import datetime
import matplotlib.dates as mdates

# Custom functions
sys.path.append('/home/users/samuel.sancar/Heat_Extremes/Code/Functions/')
from Lat_Land_mask_function import lat_land_mask
from Save_figure_function import save_figure

# CANARI member
member = sys.argv[1]

# Load reanalysis cubes
era5 = iris.load_cube('/data/users/samuel.sancar/Heat_extremes/Trends/Reanalysis_trends/era5_daily_maximum_2m_temperature_tailwidth_trends_1958-2023.nc')
jra55 = iris.load_cube('/data/users/samuel.sancar/Heat_extremes/Trends/Reanalysis_trends/jra55__minmax_surf.015_tmax.reg_tl319._tailwidth_trends_1958-2023.nc')

# Load CANARI files
CANARI_files = sorted(glob.glob('/data/users/samuel.sancar/Heat_extremes/Trends/CANARI_trends/*_tailwidth_trend.nc'))

# Regions and titles
regtitles = ['Northeast United States', 'Southwest Canada', 'Southeast Brazil', 'Paraguay', 'Northwest Africa', 'Eastern Central China', 'Northeast Russia']
reglats = [[32, 46], [51, 61], [-22, -12], [-29, -19],[15, 26],[29,33], [63, 71]]
reglons = [[-81, -70], [-125, -107], [-47, -39], [-61, -49], [-17,-6],[111,118], [131, 180]]

cmap_rb = mpl_cm.get_cmap('coolwarm', 15)

datadir = '/data/users/samuel.sancar/Heat_extremes/Trends/CANARI_trends/'
filename = f"{member}_tailwidth_trend.nc"
cube = iris.load_cube(datadir + filename)
p87p5_trend = iris.load_cube(datadir + filename.replace('tailwidth', 'p87p5'))

masked_cube, _ = lat_land_mask(cube, -60, 90)
p87p5_masked, _ = lat_land_mask(p87p5_trend, -60, 90)

# Extract grid points for stippling
lon = cube.coord('longitude').points
lat = cube.coord('latitude').points
lon2d, lat2d = np.meshgrid(lon, lat)

# Negative p87p5 trend points
neg_indices = np.where(p87p5_masked.data < 0)
neg_lats = lat2d[neg_indices]
neg_lons = lon2d[neg_indices]

# Significance mask
sig_datadir = '/data/users/samuel.sancar/Heat_extremes/Trends/Trend_significance/CANARI/P-vals/'
statsig_ds = xr.open_dataset(f"{sig_datadir}{member}_tailwidths_p_1958.nc")
statsig = list(statsig_ds.data_vars.values())[0]
sig_mask_np = statsig.values < 0.05
stipple_lats, stipple_lons = lat2d[sig_mask_np], lon2d[sig_mask_np]

fig = plt.figure(figsize=(12, 6))
ax = plt.axes(projection=ccrs.Robinson(central_longitude=0))
ax.set_extent([-180, 180, -60, 90], crs=ccrs.PlateCarree())

# Plot tailwidth trend
iplt.pcolormesh(masked_cube, cmap=cmap_rb, vmin=-0.3, vmax=0.3)

# Overlay negative trend points
ax.scatter(neg_lons, neg_lats, color='white', s=0.9, alpha=0.6, transform=ccrs.PlateCarree(), label='p87p5 trend < 0')

# Draw region boxes
for title, lat_range, lon_range in zip(regtitles, reglats, reglons):
    rect = Rectangle((lon_range[0], lat_range[0]),
                     lon_range[1] - lon_range[0],
                     lat_range[1] - lat_range[0],
                     linewidth=1, edgecolor='black', facecolor='none',
                     transform=ccrs.PlateCarree())
    ax.add_patch(rect)

ax.coastlines(resolution='50m', linewidth=0.5)
ax.add_feature(cfeature.LAND, facecolor='lightgray')
ax.add_feature(cfeature.OCEAN, facecolor='white')
ax.add_feature(cfeature.LAKES.with_scale('110m'), edgecolor='0', facecolor='none', lw=.05)
ax.add_feature(cfeature.RIVERS.with_scale('110m'), edgecolor='tab:blue', facecolor='none', lw=.4)

ax.scatter(stipple_lons, stipple_lats, color='black', s=0.4, alpha=0.5, transform=ccrs.PlateCarree(), label='Significant')

gl = ax.gridlines(draw_labels=True, linewidth=0.5, color='gray', alpha=0.5, linestyle='--')
gl.top_labels = False
gl.right_labels = False

cbar = plt.colorbar(ax.collections[0], orientation='horizontal', pad=0.05, aspect=50)
cbar.set_label('Trend in Tailwidth (°C per decade)', fontsize=12)

plt.title(f"Trend in Tailwidth for CANARI member {member} (1958–2022)", fontsize=14, weight='bold')

plt.tight_layout()
save_figure(plt, f'CANARI_{member}_tailwidths_p', 'Maps')
plt.show()

def get_area_weights(cube):
    """
    Compute area weights for a cube, ensuring latitude and longitude bounds exist.
    """
    for coord in ['latitude', 'longitude']:
        if not cube.coord(coord).has_bounds():
            cube.coord(coord).guess_bounds()
    return iac.area_weights(cube)


def areaMean(cube):
    """
    Compute weighted area mean for a cube, ignoring NaNs and masked values.
    """
    
    cube.data = np.ma.masked_invalid(cube.data)
    
    collapse_dims = [dim for dim in ['latitude', 'longitude'] if dim in [coord.name() for coord in cube.coords()]]

    if collapse_dims:
        # Compute weights
        weights = get_area_weights(cube)
      
        if np.ma.is_masked(cube.data):
            weights = np.ma.masked_where(cube.data.mask, weights)

        cube = cube.collapsed(collapse_dims, ia.MEAN, weights=weights)

    return cube

# Land-sea mask
lsm = xr.open_dataset('lsm_1279l4_0.1x0.1.grb_v4_unpack.nc')['lsm'].squeeze()

# Prepare figure
n_regions = len(regtitles)

plt.style.use('ggplot')

fig, axes = plt.subplots(n_regions, 2, figsize=(14, 4 * n_regions), sharex='col')
years = np.arange(1958, 2023)

positions = [1, 2, 3, 4]
labels = ['CANARI', 'JRA-55', 'ERA5', f'Member {member}']
member_cubes = [(iris.load_cube(fname), fname) for fname in CANARI_files]

all_trends = []

for i, region_name in enumerate(regtitles):
    ax_line = axes[i, 0]
    ax_box = axes[i, 1]

    # Reanalysis regional averages
    jra55_region = jra55.intersection(latitude=reglats[i], longitude=reglons[i])
    era5_region = era5.intersection(latitude=reglats[i], longitude=reglons[i])
    jra55_avg = areaMean(jra55_region)
    era5_avg = areaMean(era5_region)

    trend_array = []
    pseudo_obs = None

    # Plot CANARI members
    for trend_cube, fname in member_cubes:
        member_name = fname.split('/')[-1].split('_')[0]
        lat_points = trend_cube.coord('latitude').points
        lon_points = trend_cube.coord('longitude').points
        lsm_interp = lsm.interp(latitude=lat_points, longitude=lon_points, method="nearest")
        trend_cube.data = np.ma.masked_where(lsm_interp.values <= 0.5, trend_cube.data)

        region = trend_cube.intersection(latitude=reglats[i], longitude=reglons[i], ignore_bounds=True)
        avg_trend = areaMean(region)

        if avg_trend.data.size != 0:
            fit_line = avg_trend.data / 10 * (years - years[0])
            trend_array.append(avg_trend.data)
            all_trends.append(avg_trend.data)

            if member_name == member:
                ax_line.plot(years, fit_line, linewidth=2.5, color='red', label=f'Member: {member_name}')
                pseudo_obs = avg_trend
            else:
                ax_line.plot(years, fit_line, linewidth=0.8, alpha=0.3, color='grey')

    # Plot ERA5/JRA55
    if jra55_avg.data.size != 0:
        jra55_line = jra55_avg.data / 10 * (years - years[0])
        ax_line.plot(years, jra55_line, linewidth=2.0, color='green', label='JRA-55')
        all_trends.append(jra55_avg.data)

    if era5_avg.data.size != 0:
        era5_line = era5_avg.data / 10 * (years - years[0])
        ax_line.plot(years, era5_line, linewidth=2.0, color='orange', label='ERA5')
        all_trends.append(era5_avg.data)

    # Line plot styling
    ax_line.set_title(region_name, fontsize=14, fontweight='bold')
    ax_line.tick_params(axis='both', labelsize=10)
    ax_line.grid(True)

    # Boxplots
    if len(trend_array) > 0:
        bp = ax_box.boxplot(trend_array, positions=[1], widths=0.5)     

    if jra55_avg.data.size > 0:
        ax_box.scatter(2, jra55_avg.data, 
                       color='green', marker='x')
        
    if era5_avg.data.size > 0:
        ax_box.scatter(3, era5_avg.data,
                       color='orange', marker = 'x')
    if pseudo_obs is not None:
        ax_box.scatter(4, pseudo_obs.data,
                       color='red', marker = 'x')

    ax_box.set_xticks(positions)
    ax_box.set_xticklabels(labels, fontsize=10, rotation=0)
    ax_box.tick_params(axis='both', labelsize=10)
    ax_box.grid(True)

# Shared x-axis for line plots
axes[-1, 0].set_xlabel('Year', fontsize=12)
axes[2, 0].set_ylabel('Average Yearly tailwidth (°C)')
axes[2, 1].set_ylabel('Average Trend (°C/decade)')

# Consistent y-axis for boxplots
global_ylim = (min(all_trends) - 0.05, max(all_trends) + 0.05)
for i in range(n_regions):
    axes[i, 1].set_ylim(global_ylim)

# Global legend
handles, lbls = axes[0, 0].get_legend_handles_labels()
fig.legend(handles, lbls, loc='outside upper left', bbox_to_anchor=(0.735, 0.995), ncol=4, fontsize=10)

# Figure title
fig.suptitle(f'Regional Tailwidth Trends (CANARI member {member})', fontsize=14, fontweight='bold')

#plt.tight_layout(rect=[0, 0.03, 1, 0.95])
plt.savefig(f"Beautiful_Tailwidth_Trends_{member}.png", dpi=300)
plt.show()
































