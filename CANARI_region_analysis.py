""" 
17/11/2025
----------

CANARI regional analysis. Current set up of code is specific to CANARI member 29.
Bounds for latitudes and longitudes to analyse over are determined by CANARI_regional_selection code.
Boxes are determined based on regions considered as robust heatwave hot spots, passing all conditions, 
as specified in the methods.

To run: python CANARI_region_analysis.py [<member number>]
"""

import iris
import sys
import glob
import matplotlib.pyplot as plt
import iris.plot as iplt
import cartopy.crs as ccrs
import matplotlib.cm as mpl_cm
import numpy as np
import xarray as xr
import string
import matplotlib.patheffects as pe
from matplotlib.patches import Rectangle
from Region_average_Function import regional_tailwidths
from MK_functions import regional_mann_kendall_test
from Trend_calculation_functions import regional_trends
from Save_figure_function import save_figure
from Lat_Land_mask_function import lat_land_mask
from Map_features import map_features

def index_to_letters(i: int) -> str:
    alphabet = string.ascii_lowercase
    s = ""
    n = i + 1   
    while True:
        n, r = divmod(n, 26)
        s = alphabet[r] + s
        if n == 0:
            break
        n -= 1
    return s

# Beginning and end year+1
t1 = 1958
t2 = 2023

lat_minimum = -60
lat_maximum = 90

member = int(sys.argv[1])

template = iris.load_cube('template_cube_1deg.nc')
era5 = iris.load_cube(f'/data/users/samuel.sancar/Heat_extremes/Trends/Reanalysis_trends/era5_daily_maximum_2m_temperature_tailwidth_trends_{t1}-{t2-1}.nc')
jra55 = iris.load_cube(f'/data/users/samuel.sancar/Heat_extremes/Trends/Reanalysis_trends/jra55__minmax_surf.015_tmax.reg_tl319._tailwidth_trends_{t1}-{t2-1}.nc')

# Regions and titles for member 32
#regtitles = ['Northeast United States', 'Southwest Canada', 'Southeast Brazil', 'Paraguay', 'Northwest Africa', 'Eastern Central China', 'Northeast Russia'] # , 'Northern Summatra'
#reglats = [[32, 46], [51.5, 58.5], [-21.5, -14.5], [-29, -22], [19, 25.5], [29.5,33.5], [63, 71]] # , [0.5, 4.5]
#reglons = [[-81, -70], [-120.5, -107.5], [-46.5, -40.5], [-61, -49], [-16.5,-9.5], [111.5,117.5], [131, 149]] #, [96.5,102.5]

# Regions and titles from member 8 [fdr]
#regtitles = ['Western Mauritania', 'Southern China', 'Northern Madagascar', 'Eastern USA']
#reglats = [[16, 20], [25, 29], [-20, -16], [39, 43]]
#reglons = [[-16, -7], [111, 116], [44, 49], [-86, -75]]

# Regions and title from member 29
regtitles = ['b. Patagonia', 'c. Northern Argentina', 'd. DRC' , 'e. SE Africa', 'f. Britain and Ireland', 'g. Eastern Europe', 'h. Nepal', 'i. NW Siberia']
reglats = [[-56, -42.5], [-34.5, -22], [-5.5, 1], [-20.5, -11], [50.5, 57.5], [48, 54], [29.5, 34.5], [57.5, 64.5]]
reglons = [[-76.5, -65.5], [-68.5, -60], [13.5, 22.5], [24.5, 36.5], [-10.5, -1.5], [24.5, 34.5], [81.5, 89.5], [128, 145.5]]

cmap_rb = mpl_cm.get_cmap('coolwarm', 10)

datadir = '/data/users/samuel.sancar/Heat_extremes/Trends/CANARI_trends/'
filename = f"{member}_tailwidth_trend_{t1}-{t2-1}.nc"

cube = iris.load_cube(datadir + filename)
p87p5_trend = iris.load_cube(datadir + filename.replace('tailwidth', 'p87p5'))

masked_cube, _ = lat_land_mask(cube, template, lat_minimum, lat_maximum)
p87p5_masked, _ = lat_land_mask(p87p5_trend, template, lat_minimum, lat_maximum)

# Extract grid points for stippling
lon = cube.coord('longitude').points
lat = cube.coord('latitude').points
lon2d, lat2d = np.meshgrid(lon, lat)

# Extract grid point locations where p87p5 trend is negative for overlaying white stippling
neg_mask = np.ma.filled(p87p5_masked.data < 0, False)
neg_indices = np.where(neg_mask)

neg_lats = lat2d[neg_indices]
neg_lons = lon2d[neg_indices]

# Significance mask
statsig = iris.load_cube(f"/data/users/samuel.sancar/Heat_extremes/Trends/Trend_significance/CANARI/P-vals/{member}_p_{t1}-{t2-1}.nc")
statsig_masked, _ = lat_land_mask(statsig, template, lat_minimum, lat_maximum)
sig_mask_np = np.ma.filled(statsig_masked.data < 0.05, False)
stipple_lats, stipple_lons = lat2d[sig_mask_np], lon2d[sig_mask_np]

# Plot tailwidth trend
fig = plt.figure(figsize=(12, 6))
ax = plt.axes(projection=ccrs.Robinson(central_longitude=0))

map_features('Detail', ax, lat_minimum, lat_maximum)

map = iplt.pcolormesh(masked_cube, cmap=cmap_rb, vmin=-0.3, vmax=0.3)

# For plotting label
offset_pts = 3

# Draw region boxes
for i, (title, lat_range, lon_range) in enumerate(zip(regtitles, reglats, reglons)):
    rect = Rectangle((lon_range[0], lat_range[0]),
                     lon_range[1] - lon_range[0],
                     lat_range[1] - lat_range[0],
                     linewidth=1, edgecolor='black', facecolor='none',
                     transform=ccrs.PlateCarree())
    ax.add_patch(rect)
    
    lon_center = 0.5 * (lon_range[0] + lon_range[1])
    lat_top    = max(lat_range)
    label = index_to_letters(i)

    ax.annotate(
        label,
        xy=(lon_center, lat_top),
        xycoords=ccrs.PlateCarree()._as_mpl_transform(ax),
        xytext=(0, offset_pts), textcoords='offset points',
        ha='center', va='bottom',
        fontsize=11, fontweight='bold', color='black',
        path_effects=[pe.withStroke(linewidth=2, foreground='white')]
    )

# Plot significant trends
ax.scatter(stipple_lons, stipple_lats, color='black', s=0.6, alpha=0.6, transform=ccrs.PlateCarree(), label='Stippling')

# Pot white stippling where p87p5_trend is negative
ax.scatter(neg_lons, neg_lats, color="white", s=0.6, alpha=0.6, marker = 's', transform=ccrs.PlateCarree(), label='p87p5 trend < 0')

cbar = plt.colorbar(map, orientation='horizontal', pad=0.05, aspect=50)
cbar.set_label('Trend in Tail width (°C per decade)', fontsize=13)

plt.title(f"a) CANARI-LE member {member} tail width trends ({t1}-{t2-1})", fontsize=15, weight='bold')
plt.tight_layout()
save_figure(plt, f'CANARI_{member}_tailwidths_p', 'Maps')
plt.show()

# Land-sea mask
lsm = xr.open_dataset('lsm_1279l4_0.1x0.1.grb_v4_unpack.nc')['lsm'].squeeze()

template_lat_points = template.coord('latitude').points
template_lon_points = template.coord('longitude').points
lsm_interp = lsm.interp(latitude=template_lat_points, longitude=template_lon_points, method="nearest")

n_regions = len(regtitles)

plt.style.use('ggplot')

ncols = 4
nrows = int(np.ceil(n_regions / ncols))

fig, axes = plt.subplots(
    nrows, ncols,
    figsize=(5 * ncols + 4, 11),
    sharex=False, sharey=False
)

plt.subplots_adjust(wspace=0.45, hspace=0.35)

years = np.arange(t1, t2)

positions = [1,2]
labels = ['CANARI-LE', 'Reanalyses']

all_trends = []
region_pvals = []
percentiles  = []
box_axes_used = []
line_axes_used = []

# Colours
line_color_member = 'red'
line_color_era    = 'black'
line_color_jra    = 'black'
line_color_others = 'grey'

first_line_ax = None

# Iterate over all regions 
for i, region_name in enumerate(regtitles):
    
    print(f'Processing region {region_name}')
    
    if i <= 3:
        ax_box  = axes[0, i]  
    else:
        ax_box = axes[1, i-4]

    # Reanalysis regional averages (series and trend)
    era5_series  = regional_tailwidths('era5',  reglats[i], reglons[i], lat_minimum, lat_maximum, t1, t2)
    jra55_series = regional_tailwidths('jra55', reglats[i], reglons[i], lat_minimum, lat_maximum, t1, t2)
    era5_trends  = regional_trends(era5_series)   
    jra55_trends = regional_trends(jra55_series)  

    trend_array = []   # CANARI member trends for the box
    pseudo_obs  = None

    #Plot CANARI members (lines) and collect trends for box 
    for member_name in np.arange(1, 41):

        print(f'Processing {member_name} for region {region_name}')
        
        # Regional series and trend
        region_series = regional_tailwidths(member_name, reglats[i], reglons[i], lat_minimum, lat_maximum, t1, t2)
        avg_trend     = regional_trends(region_series)  
        
        if np.size(avg_trend) != 0:
            if member_name != member:
                # Line for fitted trend across years
                fit_line = avg_trend / 10.0 * (years - years[0])  
                trend_array.append(avg_trend)
                all_trends.append(avg_trend)

            if member_name == member:
                # Mann-Kendall p-value for the series
                _, p_val = regional_mann_kendall_test(region_series)
                region_pvals.append((region_name, np.round(p_val, 4)))

                # Trend line 
                # Overlay series (centered to match the fitted line’s mean)
                mean_fit  = np.nanmean(fit_line)
                mean_data = np.nanmean(region_series.data)
                delta = mean_data - mean_fit
                pseudo_obs = avg_trend

    # ERA5 & JRA-55 fitted lines
    if np.size(jra55_trends) != 0:
        jra55_line = jra55_trends / 10.0 * (years - years[0])
        all_trends.append(jra55_trends)

    if np.size(era5_trends) != 0:
        era5_line = era5_trends / 10.0 * (years - years[0])
        all_trends.append(era5_trends)

    # Boxplot
    ax_box.set_title(region_name, fontsize=17, fontweight='bold', pad=8)
    if len(trend_array) > 0:
        ax_box.boxplot(
            trend_array,
            positions=[1],
            widths=0.85,  # wider box
            boxprops={'linewidth': 2},
            whiskerprops={'linewidth': 1.6},
            capprops={'linewidth': 1.6},
            medianprops={'linewidth': 2, 'color': 'black'}
        )

    # Plot ERA5, JRA-55 and pseudo-observation
    if np.size(era5_trends) != 0:
        ax_box.scatter(2, era5_trends, color=line_color_era, marker='D', s=250, linewidths=1.8, label = 'ERA5')
    if np.size(jra55_trends) != 0:
        ax_box.scatter(2, jra55_trends, color=line_color_jra, marker='*', s=250, linewidths=1.8, label = 'JRA-55')
    if pseudo_obs is not None:
        ax_box.scatter(1, pseudo_obs, color=line_color_member, marker='*', s=260, linewidths=1.8, label = 'Member 29')

    ax_box.set_xticks(positions)
    ax_box.set_xticklabels(labels, fontsize=20, rotation=0)
    ax_box.tick_params(axis='both', labelsize=20)
    ax_box.grid(True)
    
    # Calculate JRA-55 and ERA5 percentile within model spread.
    era5_area_percentile = (np.sum(trend_array < era5_trends))/ len(trend_array)*100
    jra_are_percentile =  (np.sum(trend_array < jra55_trends))/ len(trend_array)*100
    pseudo_obs_percentile = (np.sum(trend_array < pseudo_obs))/ len(trend_array)*100
    
    percentiles.append((region_name, 
                        np.round(era5_area_percentile,2), 
                        np.round(jra_are_percentile,2),
                        np.round(pseudo_obs_percentile,2)))

axes[0, 0].set_ylabel('Average tail width \n trend (°C/decade)', fontsize=22)
axes[1, 0].set_ylabel('Average tail width\n trend (°C/decade)', fontsize=22)

handles, lbls = axes[0, 0].get_legend_handles_labels()
fig.legend(handles, ['ERA5', 'JRA-55', 'Member 29'],
           loc='upper left', bbox_to_anchor=(0.70, 0.995), ncol=3, fontsize=20)

fig.suptitle(f'Regional Tailwidth Trends ({t1}-{t2-1})', fontsize=18, fontweight='bold', y=1.06)

length_str = str(n_regions)
save_figure(plt, f"CANARI_regional_trends_member_{member}.png", 'Graphs')
plt.show()

print(f"Regional Mann-Kendall p-values: {region_pvals}")
print(f"Regional Percentiles (ERA5, JRA-55, Member {member}): {percentiles}")


















