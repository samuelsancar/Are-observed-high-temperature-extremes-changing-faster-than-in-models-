"""
23/09/2025
-----------

Histogram of FDR correction for tailwidth trends for all models, ERA5 and JRA-55 reanalysis data.

To run: python Model_FDR_histogram.py
"""

import iris
import iris.analysis
import cf_units
import sys
import matplotlib.pyplot as plt
import numpy as np
import glob
import os
from Save_figure_function import save_figure

t1 = 1958
t2 = 2023

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

variable = 'daily_maximum_2m_temperature'
reanalysis = 'era5'
significant = sys.argv[1]  # Choose from 'positive' or 'all', 

# Load template cube
template_cube = iris.load_cube('template_cube_1deg.nc')

if significant == 'positive':
    graph_name = 'positive trends'
    significance_type = 'FDR_corrections_PositiveOnly'
    selection = 'tailwidth_fdr_positive'
    jra_name = f'jra55__minmax_surf.015_tmax.reg_tl319._tailwidth_fdr_positive_{t1}-{t2}.nc'
    era5_name = 'tailwidth_fdr_positive' 

else:
    graph_name = 'trends'
    significance_type = 'FDR_corrections_All'
    selection = 'fdr_significant_all'
    jra_name = f'jra55__minmax_surf.015_tmax.reg_tl319._tailwidth_fdr_all_{t1}-{t2}.nc'
    era5_name = 'tailwidth_fdr_all'
    
# Load and regrid ERA5 trend data
datadir = '/data/users/samuel.sancar/Heat_extremes/Trends/Reanalysis_trends/'
era5_cube = iris.load_cube(f"{datadir}{reanalysis}_{variable}_{era5_name}_{t1}-{t2}.nc")

datadir_jra = '/data/users/samuel.sancar/Heat_extremes/Trends/Reanalysis_trends/'
jra_cube = iris.load_cube(f"{datadir_jra}{jra_name}")

# Load and regrid model cubes
CANARI_filenmaes = sorted(glob.glob(f'/data/users/samuel.sancar/Heat_extremes/Trends/Trend_significance/CANARI/{significance_type}/*_{selection}.nc'))

significant_gridpoints = []
significant_area =[]

for filename in CANARI_filenmaes:
    
    cube = iris.load_cube(filename)
    model_name = os.path.basename(filename)

    cube = units_and_coords(cube)

    # Area weights
    sig_mask = cube.data.astype(bool)
    sig_cube = cube.copy()
    sig_cube.data = sig_mask.astype(int)  
    sum_sig_gridpoints = np.count_nonzero(sig_mask)
    
    area_weights = iris.analysis.cartography.area_weights(sig_cube)
    
    sig_area = sig_cube.collapsed(
        ['latitude', 'longitude'],
        iris.analysis.SUM,
        weights=area_weights
    ).data

    print(f"Significant grid points in {model_name}: {sum_sig_gridpoints}")
    print(f"Significant area in {model_name}: {sig_area}")

    significant_gridpoints.append(sum_sig_gridpoints)
    significant_area.append(sig_area)

# --- ERA5 processing ---
era5_cube = units_and_coords(era5_cube)

# Boolean significance mask
era5_sig_mask = era5_cube.data.astype(bool)

era5_sig_cube = era5_cube.copy()
era5_sig_cube.data = era5_sig_mask.astype(int)
sum_era5_sig_gridpoints = np.count_nonzero(era5_sig_mask)

era5_area_weights = iris.analysis.cartography.area_weights(era5_sig_cube)
# Area
era5_sig_area = era5_sig_cube.collapsed(
    ['latitude', 'longitude'],
    iris.analysis.SUM,
    weights=era5_area_weights
).data

# --- JRA-55 processing ---
jra_cube = units_and_coords(jra_cube)

# Boolean significance mask
jra_sig_mask = jra_cube.data.astype(bool)

jra_sig_cube = jra_cube.copy()
jra_sig_cube.data = jra_sig_mask.astype(int)
sum_jra_sig_gridpoints = np.count_nonzero(jra_sig_mask)

jra_area_weights = iris.analysis.cartography.area_weights(jra_sig_cube)

# Area
jra_sig_area = jra_sig_cube.collapsed(
    ['latitude', 'longitude'],
    iris.analysis.SUM,
    weights=jra_area_weights
).data

jra_sig_gridpoints = jra_cube.data #< 0.05
sum_jra_sig_gridpoints = np.count_nonzero(jra_sig_gridpoints)

# Calculate JRA-55 and ERA5 percentile within model spread.
significant_gridpoints_array = np.array(significant_gridpoints)
significant_area_array = np.array(significant_area)
era5_percentile = (np.sum(significant_gridpoints_array < sum_era5_sig_gridpoints))/ len(significant_gridpoints_array)*100
jra_percentile =  (np.sum(significant_gridpoints_array < sum_jra_sig_gridpoints))/ len(significant_gridpoints_array)*100
era5_area_percentile = (np.sum(significant_area_array < era5_sig_area))/ len(significant_area_array)*100
jra_are_percentile =  (np.sum(significant_area_array < jra_sig_area))/ len(significant_area_array)*100

print(f"ERA5 significant grid points: {sum_era5_sig_gridpoints}, Percentile: {era5_percentile:.2f}%")
print(f"JRA-55 significant grid points: {sum_jra_sig_gridpoints}, Percentile: {jra_percentile:.2f}%")

print(f"ERA5 significant area: {era5_sig_area}, Percentile: {era5_area_percentile:.2f}%")
print(f"JRA-55 significant area: {jra_sig_area}, Percentile: {jra_are_percentile:.2f}%")

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

plt.xlabel('Area of significant trends (m\u00B2)', fontsize=18, fontweight='bold')
plt.ylabel('Number of Models', fontsize=18, fontweight='bold')
plt.title(f'CANARI distribution area-weighted\nsignificant {graph_name} ', fontsize=20, fontweight='bold')

plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.grid(True, linestyle='--', alpha=0.5, linewidth=1.2)

#plt.legend(fontsize=13, loc='upper right', frameon=False)
plt.tight_layout()
save_figure(plt, f'CANARI_fdr_{graph_name}_histogram', 'Graphs')
plt.show()