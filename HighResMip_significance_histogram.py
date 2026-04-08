"""

26/08/25

Plot histogram for significant tail widths in Models and ERA5 trends. Land masks have already been applied to models.

To run: python plot_tailwidth_trend.py [1958] [2023]

"""

import iris
import iris.analysis
import cf_units
import sys
import iris.quickplot as iplt
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import numpy as np
import glob
import os
from Save_figure_function import save_figure
from Lat_Land_mask_function import lat_land_mask

variable = 'daily_maximum_2m_temperature'
reanalysis = 'era5'
significant = sys.argv[1]  # Choose from 'positive' or 'all', 

def units_and_coords(cube):
    if not cube.coord('latitude').has_bounds():
        cube.coord('latitude').guess_bounds()
    if not cube.coord('longitude').has_bounds():
        cube.coord('longitude').guess_bounds()
        
    # Fix units if missing or unknown
    if cube.coord('latitude').units.is_unknown():
        cube.coord('latitude').units = cf_units.Unit('degrees')
    if cube.coord('longitude').units.is_unknown():
        cube.coord('longitude').units = cf_units.Unit('degrees')
        
    return cube

# Load template cube
template_cube = iris.load_cube('template_cube_1deg.nc')

if significant == 'positive':
    graph_name = 'positive trends'
    significance_type = 'Positive_Pvals'
    Model_name = 'pos_p'
    jra_name = 'jra55__minmax_surf.015_tmax.reg_tl319._tailwidth_pos_p_1958-2023.nc'
    era5_name = 'tailwidth_pos_p' 
else:
    graph_name = 'trends'
    significance_type = 'P-vals'
    Model_name = 'p'
    jra_name = 'jra55__minmax_surf.015_tmax.reg_tl319._tailwidth_p_1958-2023.nc'
    era5_name = 'tailwidth_p'
    
# Load and regrid ERA5 trend data
datadir = '/data/users/samuel.sancar/Heat_extremes/Trends/Reanalysis_trends/'
era5_cube = iris.load_cube(f"{datadir}{reanalysis}_{variable}_{era5_name}_1958-2023.nc")

datadir_jra = '/data/users/samuel.sancar/Heat_extremes/Trends/Reanalysis_trends/'
jra_cube = iris.load_cube(f"{datadir_jra}{jra_name}")

era5_masked_cube = era5_cube
era5_data = era5_masked_cube.data.flatten()
era5_data = era5_data[~np.isnan(era5_data)]

# Load and regrid model cubes
coupled_filenames = sorted(glob.glob(f'/data/users/samuel.sancar/Heat_extremes/Trends/Trend_significance/{significance_type}/*_hist-1950_*_{Model_name}_1958.nc'))
sst_filenames = sorted(glob.glob(f'/data/users/samuel.sancar/Heat_extremes/Trends/Trend_significance/{significance_type}/*_highresSST-present_*_{Model_name}_1958.nc'))

# Total valid land area (Based off ERA5)
benchmark_unmasked = iris.load_cube(f"/data/users/samuel.sancar/Heat_extremes/ERA5_tailwidths/era5_daily_maximum_2m_temperature_1958_tail_width.nc")
benchmark, total_gridpoints = lat_land_mask(benchmark_unmasked, -60, 90)
benchmark = units_and_coords(benchmark)
area_weights = iris.analysis.cartography.area_weights(benchmark)
masked_area = np.where(np.isnan(benchmark.data), 0, area_weights)
total_area = masked_area.sum()

print(f'Total are = {total_area}')

Total_models = len(coupled_filenames) + len(sst_filenames)
significant_gridpoints = []
significant_area =[]

for filename in coupled_filenames:
    
    cube = iris.load_cube(filename)
    model_name = os.path.basename(filename)

    cube = units_and_coords(cube)

    area_weights = iris.analysis.cartography.area_weights(cube)
    sig_mask = cube.data < 0.05
    sig_area = np.ma.array(area_weights, mask=~sig_mask).sum()
    sum_sig_gridpoints = np.count_nonzero(sig_mask)
    
    coupled_proportion_area = (sig_area / total_area) * 100
    
    print(f"Significant grid points in {model_name}: {sum_sig_gridpoints}")
    print(f"Significant area in {model_name}: {coupled_proportion_area}%")
    significant_gridpoints.append(sum_sig_gridpoints)
    significant_area.append(sig_area)

for filename in sst_filenames:
    
    cube = iris.load_cube(filename)
    model_name = os.path.basename(filename)

    cube = units_and_coords(cube)

    area_weights = iris.analysis.cartography.area_weights(cube)
    sig_mask = cube.data < 0.05
    sig_area = np.ma.array(area_weights, mask=~sig_mask).sum()
    sum_sig_gridpoints = np.count_nonzero(sig_mask)
    
    sst_proportion_area = (sig_area / total_area) * 100
    
    print(f"Significant grid points in {model_name}: {sum_sig_gridpoints}")
    print(f"Significant area in {model_name}: {sst_proportion_area}%")
    significant_gridpoints.append(sum_sig_gridpoints)
    significant_area.append(sig_area)

# Looking at proportion of positive trends
trends = iris.load_cube('/data/users/samuel.sancar/Heat_extremes/Trends/Reanalysis_trends/era5_daily_maximum_2m_temperature_tailwidth_trends_1958-2023.nc')
trends, total_gridpoints = lat_land_mask(trends, -60, 90)

trends = units_and_coords(trends)

positive_area_weighting = iris.analysis.cartography.area_weights(trends)
positive_sig_mask = trends.data > 0 
positive_area = np.ma.array(positive_area_weighting, mask=~positive_sig_mask).sum()
positive_era5_sig_gridpoints = np.count_nonzero(positive_sig_mask)

# Process ERA5
era5_cube = units_and_coords(era5_cube)

era5_area_weights = iris.analysis.cartography.area_weights(era5_cube)
era5_sig_mask = era5_cube.data < 0.05
era5_sig_area = np.ma.array(era5_area_weights, mask=~era5_sig_mask).sum()
sum_era5_sig_gridpoints = np.count_nonzero(era5_sig_mask)

#significant_gridpoints.append(sum_era5_sig_gridpoints)

jra_cube = units_and_coords(jra_cube)

jra_area_weights = iris.analysis.cartography.area_weights(jra_cube)
jra_sig_mask = jra_cube.data < 0.05
jra_sig_area = np.ma.array(jra_area_weights, mask=~jra_sig_mask).sum()
sum_jra_sig_gridpoints = np.count_nonzero(jra_sig_mask)
#significant_gridpoints.append(sum_jra_sig_gridpoints)

both_mask = era5_sig_mask & jra_sig_mask

both_sig_cube = jra_cube.copy()
both_sig_cube.data = both_mask.astype(int)
sum_both_sig_gridpoints = np.count_nonzero(both_mask)

both_area_weights = iris.analysis.cartography.area_weights(both_sig_cube)

# Area
both_sig_area = both_sig_cube.collapsed(
    ['latitude', 'longitude'],
    iris.analysis.SUM,
    weights=both_area_weights
).data

sum_both_sig_gridpoints = np.count_nonzero(both_sig_cube)

# Calculate JRA-55 and ERA5 percentile within model spread.
significant_gridpoints_array = np.array(significant_gridpoints)
significant_area_array = np.array(significant_area)
era5_percentile = (np.sum(significant_gridpoints_array < sum_era5_sig_gridpoints))/ len(significant_gridpoints_array)*100
jra_percentile =  (np.sum(significant_gridpoints_array < sum_jra_sig_gridpoints))/ len(significant_gridpoints_array)*100
era5_area_percentile = (np.sum(significant_area_array < era5_sig_area))/ len(significant_area_array)*100
jra_are_percentile =  (np.sum(significant_area_array < jra_sig_area))/ len(significant_area_array)*100
both_area_percentile = (np.sum(significant_area_array < both_sig_area))/ len(significant_area_array)*100

proportion_gridpoints = (sum_era5_sig_gridpoints / total_gridpoints) * 100
proportion_area = (era5_sig_area / total_area) * 100
jra_proportion_area = (jra_sig_area / total_area) * 100
proportion_positive = (positive_area / total_area) * 100
model_average_sig_area = np.mean(significant_area_array)
proportion_model_average_area = (model_average_sig_area / total_area)*100


print(f"ERA5 significant grid points: {sum_era5_sig_gridpoints}, Percentile: {era5_percentile:.2f}%")
print(f"JRA-55 significant grid points: {sum_jra_sig_gridpoints}, Percentile: {jra_percentile:.2f}%")
print(f"Proportion of land gridpoints that is significant for ERA5: {np.round(proportion_gridpoints, 2)}%")
print(f"Proportion of land area that is significant for ERA5: {np.round(proportion_area, 2)}%")
print(f"Proportion of land area that is significant for JRA-55: {np.round(jra_proportion_area, 2)}%")
print(f"Proportion of land area that has positive trends for ERA5: {np.round(proportion_positive)}%")
print(f"ERA5 significant area: {era5_sig_area}, Percentile: {era5_area_percentile:.2f}%")
print(f"JRA-55 significant area: {jra_sig_area}, Percentile: {jra_are_percentile:.2f}%")
print(f"Both significant area: {both_sig_area}, Percentile: {both_area_percentile:.2f}%")
print(f"Model average significant area: {model_average_sig_area}, Proportion of land area : {np.round(proportion_model_average_area)}%")
print("Total number of models:", Total_models)

plt.style.use('ggplot')

plt.figure(figsize=(12, 7))
counts, bins, patches = plt.hist(
    significant_gridpoints,
    bins=10,
    range=(significant_gridpoints_array.min(), significant_gridpoints_array.max()),
    histtype='bar',
    color='#4682B4', 
    alpha=0.85,
    edgecolor='white',
    linewidth=2
)

plt.axvline(x=sum_era5_sig_gridpoints, color='#D7263D', linestyle='--', linewidth=2.5, label='ERA5 Value')
plt.text(
    sum_era5_sig_gridpoints + 50, max(counts) * 0.92, 'ERA5',
    color='#D7263D', fontsize=15, fontweight='bold',
    bbox=dict(facecolor='white', edgecolor='#D7263D', boxstyle='round,pad=0.3', alpha=0.8)
)

plt.axvline(x=sum_jra_sig_gridpoints, color='#1B998B', linestyle='--', linewidth=2.5, label='JRA-55 Value')
plt.text(
    sum_jra_sig_gridpoints + 50, max(counts) * 0.82, 'JRA-55',
    color='#1B998B', fontsize=15, fontweight='bold',
    bbox=dict(facecolor='white', edgecolor='#1B998B', boxstyle='round,pad=0.3', alpha=0.8)
)

plt.axvline(x=sum_both_sig_gridpoints, color='#1B998B', linestyle='--', linewidth=2.5, label='ERA5 & JRA-55 agreement')
plt.text(
    sum_both_sig_gridpoints + 50, max(counts) * 0.82, 'ERA5 & JRA-55 agreement',
    color='#1B998B', fontsize=15, fontweight='bold',
    bbox=dict(facecolor='white', edgecolor='#1B998B', boxstyle='round,pad=0.3', alpha=0.8)
)

plt.xlabel('Number of Significant Grid Points', fontsize=18, fontweight='bold')
plt.ylabel('Number of Models', fontsize=18, fontweight='bold')
plt.title(f'Histogram of Grid Points with significant {graph_name} (1958–2022)', fontsize=20, fontweight='bold')

plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.grid(True, linestyle='--', alpha=0.5, linewidth=1.2)

#plt.legend(fontsize=13, loc='upper right', frameon=False)
plt.tight_layout()
plt.show()

# Histogram for significant area
plt.figure(figsize=(12, 7))
counts, bins, patches = plt.hist(
    significant_area,
    bins=10,
    range=(significant_area_array.min(), significant_area_array.max()),
    histtype='bar',
    color='#4682B4', 
    alpha=0.85,
    edgecolor='white',
    linewidth=2
)

plt.axvline(x=era5_sig_area, color='#D7263D', linestyle='--', linewidth=2.5, label='ERA5 Value')
plt.text(
    era5_sig_area + 50, max(counts) * 0.92, 'ERA5',
    color='#D7263D', fontsize=15, fontweight='bold',
    bbox=dict(facecolor='white', edgecolor='#D7263D', boxstyle='round,pad=0.3', alpha=0.8)
)

# Label histogram at JRA-55 value
plt.axvline(x=jra_sig_area, color='#1B998B', linestyle='--', linewidth=2.5, label='JRA-55 Value')
plt.text(
    jra_sig_area + 50, max(counts) * 0.82, 'JRA-55',
    color='#1B998B', fontsize=15, fontweight='bold',
    bbox=dict(facecolor='white', edgecolor='#1B998B', boxstyle='round,pad=0.3', alpha=0.8)
)

plt.axvline(x=both_sig_area, color='#1B998B', linestyle='--', linewidth=2.5, label='ERA5 & JRA-55 Agreement')
plt.text(
    both_sig_area + 50, max(counts) * 0.82, 'ERA5 & JRA-55 Agreement',
    color='#1B998B', fontsize=15, fontweight='bold',
    bbox=dict(facecolor='white', edgecolor='#1B998B', boxstyle='round,pad=0.3', alpha=0.8)
)

plt.xlabel('Area of significant trends (m\u00B2)', fontsize=18, fontweight='bold')
plt.ylabel('Number of Models', fontsize=18, fontweight='bold')
plt.title(f'HighResMIP distribution area-weighted\nsignificant {graph_name} ', fontsize=20, fontweight='bold')

plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.grid(True, linestyle='--', alpha=0.5, linewidth=1.2)

#plt.legend(fontsize=13, loc='upper right', frameon=False)
plt.tight_layout()
save_figure(plt, f'HighResMIP_Histogram_significant_gridpoints', 'Graphs')
plt.show()



