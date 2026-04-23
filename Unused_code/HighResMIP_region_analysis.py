""" 
HighResMIP regional analysis.
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
from Region_average_Function import regional_tailwidths, MIP_regional_tailwidths
from Trend_calculation_functions import regional_trends
from MK_functions import regional_mann_kendall_test

# Custom functions
sys.path.append('/home/users/samuel.sancar/Heat_Extremes/Code/Functions/')
from Lat_Land_mask_function import lat_land_mask
from Save_figure_function import save_figure

# HighResMIP member
member = 'EC-Earth3P-HR_highresSST-present_r3i1p1f1'

# Load reanalysis cubes
era5 = iris.load_cube('/data/users/samuel.sancar/Heat_extremes/Trends/Reanalysis_trends/era5_daily_maximum_2m_temperature_tailwidth_trends_1958-2023.nc')
jra55 = iris.load_cube('/data/users/samuel.sancar/Heat_extremes/Trends/Reanalysis_trends/jra55__minmax_surf.015_tmax.reg_tl319._tailwidth_trends_1958-2023.nc')

coupled = sorted(glob.glob('/data/users/samuel.sancar/Heat_extremes/Trends/Coupled_model_trends/*tailwidth_trend.nc'))
sst = sorted(glob.glob('/data/users/samuel.sancar/Heat_extremes/Trends/SST_model_trends/*tailwidth_trend.nc'))

coupled_files = [(iris.load_cube(fname), fname) for fname in coupled]
sst_files = [(iris.load_cube(fname), fname) for fname in sst]

all_files = [("Coupled", cube, fname) for cube, fname in coupled_files] + \
            [("SST", cube, fname) for cube, fname in sst_files]
            
# Regions and titles
regtitles = ['Eastern Europe','Southern Niger', 'Alaska', 'Great Plains',' Southern Brazil & Paraguay', 'Papua New Guinea']
reglats = [[49, 54.5], [12.5, 17.5], [59, 66.5], [36.5, 44], [-25,-20], [-9.5, -4.5]]
reglons = [[22, 36.5], [-5.5, 14.5], [-166, -141], [-98.5, -89.5], [-65, -49.5], [137, 144.5]]

cmap_rb = mpl_cm.get_cmap('coolwarm', 15)

datadir = '/data/users/samuel.sancar/Heat_extremes/Trends/SST_model_trends/'
filename = f"EC-Earth3P-HR_highresSST-present_r3i1p1f1_tailwidth_trend.nc"
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
sig_datadir = '/data/users/samuel.sancar/Heat_extremes/Trends/Trend_significance/P-vals/'
statsig_ds = xr.open_dataset(f"{sig_datadir}EC-Earth3P-HR_highresSST-present_r3i1p1f1_tailwidths_p_1958.nc")
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

plt.title(f"Trend in Tailwidth for {member} (1958–2022)", fontsize=14, weight='bold')

plt.tight_layout()
save_figure(plt, f'{member}', 'Maps')
plt.show()

# Land-sea mask
lsm = xr.open_dataset('lsm_1279l4_0.1x0.1.grb_v4_unpack.nc')['lsm'].squeeze()

# Prepare figure
n_regions = len(regtitles)

plt.style.use('ggplot')

# Wide horizontal layout: 2 rows (top: lines, bottom: boxplots), columns = regions
# IMPORTANT: do NOT share x across rows (boxplots shouldn't use "years" scale)
fig, axes = plt.subplots(
    2, n_regions,
    figsize=(4 * n_regions + 4, 8),   # tweak width/height to taste
    sharex=False, sharey=False
)

years = np.arange(1958, 2023)
positions = [1, 2, 3]
labels = ['HighResMIP', 'ERA5', 'JRA-55']
member_cubes = all_files

all_trends = []
region_pvals = []

for i, region_name in enumerate(regtitles):
    ax_line = axes[0, i]   # top row: line plots
    ax_box  = axes[1, i]   # bottom row: boxplots

    # Reanalysis regional averages
    era5_series  = regional_tailwidths('era5',  reglats[i], reglons[i]) 
    jra55_series = regional_tailwidths('jra55', reglats[i], reglons[i])
    era5_trends  = regional_trends(era5_series)
    jra55_trends = regional_trends(jra55_series)

    trend_array = []
    pseudo_obs = None

    # Plot HighResMIP members fits
    for _, trend_cube, fname in member_cubes:
        member_name = fname.split('/')[-1].split('_tailwidth_trend.nc')[0]
        lat_points = trend_cube.coord('latitude').points
        lon_points = trend_cube.coord('longitude').points
        lsm_interp = lsm.interp(latitude=lat_points, longitude=lon_points, method="nearest")
        trend_cube.data = np.ma.masked_where(lsm_interp.values <= 0.5, trend_cube.data)

        region = MIP_regional_tailwidths(member_name, reglats[i], reglons[i])
        avg_trend = regional_trends(region)

        if avg_trend.size != 0:
            fit_line = avg_trend / 10 * (years - years[0])
            trend_array.append(avg_trend)
            all_trends.append(avg_trend)

            if member_name == member:
                _, p_val = regional_mann_kendall_test(region)
                region_pvals.append((region_name, np.round(p_val, 3)))

                # Fit line + pseudo-obs aligned series
                ax_line.plot(years, fit_line, linewidth=2.5, color='red', label=f'{member_name}')
                mean_fit  = np.nanmean(fit_line)
                mean_data = np.nanmean(region.data)
                delta = mean_data - mean_fit
                ax_line.plot(years, region.data - delta, linewidth=1.0, color='red')

                pseudo_obs = avg_trend
            else:
                ax_line.plot(years, fit_line, linewidth=0.8, alpha=0.3, color='grey')

    # ERA5/JRA55 lines
    if jra55_trends.size != 0:
        jra55_line = jra55_trends / 10 * (years - years[0])
        ax_line.plot(years, jra55_line, linewidth=2.0, color='green', label='JRA-55')
        all_trends.append(jra55_trends)

    if era5_trends.size != 0:
        era5_line = era5_trends / 10 * (years - years[0])
        ax_line.plot(years, era5_line, linewidth=2.0, color='orange', label='ERA5')
        all_trends.append(era5_trends)

    # Style line axes
    
    # Inside the loop, after plotting for each region:
    ax_line.set_title(region_name, fontsize=13, fontweight='bold', pad=15)

    ax_line.tick_params(axis='both', labelsize=10)
    ax_line.grid(True)

    # --- Boxplots (bottom row) ---
    if len(trend_array) > 0:
        ax_box.boxplot(trend_array, positions=[1], widths=0.6)

    # Slightly bigger scatter markers for reference points
    scatter_size = 60  # increased from typical defaults
    if jra55_trends.size > 0:
        ax_box.scatter(3, jra55_trends, color='green', marker='x', s=scatter_size)
    if era5_trends.size > 0:
        ax_box.scatter(2, era5_trends, color='orange', marker='x', s=scatter_size)
    if pseudo_obs is not None:
        # Place pseudo-obs as a separate marker at position 4 (to match label order)
        ax_box.scatter(1, pseudo_obs, color='red', marker='x', s=scatter_size)

    ax_box.set_xticks(positions)
    ax_box.set_xticklabels(labels, fontsize=10, rotation=0)
    ax_box.tick_params(axis='both', labelsize=10)
    ax_box.grid(True)

# --- Column titles centered above each region ---
fig.canvas.draw() 

# Axis labels (left-most column only)
axes[0, 0].set_ylabel('Average Yearly tailwidth (°C)', fontsize=13)
axes[1, 0].set_ylabel('Average Trend (°C/decade)', fontsize=13)

# Global legend from the first line axis
handles, lbls = axes[0, 0].get_legend_handles_labels()
fig.legend(handles, ['Pseudo Obs', 'JRA-55', 'ERA5'],
           loc='upper left', bbox_to_anchor=(0.7, 0.995), ncol=3, fontsize=12)

# Figure title
fig.suptitle(f'Pseudo observations regional trends:\n{member}', fontsize=15, fontweight='bold')

# Layout and save
plt.tight_layout(rect=[0, 0.03, 1, 0.93])  # leave space at top for column titles + legend
plt.savefig(f"HighResMIP_regional_trends_{member}.png", dpi=300)
# plt.show()



print(f"Regional Mann-Kendall p-values: {region_pvals}")