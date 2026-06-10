"""
26/08/25
--------

Plot histogram for significant tail widths in Models and ERA5 trends. Land masks have already been applied to models.

To run: python plot_tailwidth_trend.py 
"""

import iris
import iris.analysis
import cf_units
import sys
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
from Save_figure_function import save_figure
from Lat_Land_mask_function import lat_land_mask

# Beginning & end year+1
t1 = 1958
t2 = 2023

lat_minimum = -60
lat_maximum = 80

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
era5 = 'era5'
merra = 'MERRA'
jra = f'jra55__minmax_surf.015_tmax.reg_tl319.'
variable = 'daily_maximum_2m_temperature'

# Load template cube
template_cube = iris.load_cube('template_cube_1deg.nc')

if significant == 'positive':
    graph_name = 'positive trends'
    selection = f'_tailwidth_pos_p_{t1}-{t2-1}.nc'
    model_selection = f'_pos_p_{t1}-{t2-1}.nc'

else:
    graph_name = 'trends'
    selection = f'_tailwidth_p_{t1}-{t2-1}.nc'
    model_selection = f'_p_{t1}-{t2-1}.nc'
            
# Load and regrid ERA5 trend data
datadir = '/data/users/samuel.sancar/Heat_extremes/Trends/Reanalysis_trends/'

era5_cube = iris.load_cube(f"{datadir}{era5}_{variable}{selection}")
jra_cube = iris.load_cube(f"{datadir}{jra}{selection}")

# Total land area
template = iris.load_cube('template_cube_1deg.nc')
benchmark_unmasked = iris.load_cube(f"/data/users/samuel.sancar/Heat_extremes/ERA5_tailwidths/era5_daily_maximum_2m_temperature_1958_tail_width.nc")
benchmark, total_gridpoints = lat_land_mask(benchmark_unmasked, template, -60, 90)
benchmark = units_and_coords(benchmark)

weights2d = iris.analysis.cartography.area_weights(benchmark)
mask = np.ma.getmaskarray(benchmark.data)
weights = np.asarray(weights2d)[~mask].astype(float)
total_area = weights.sum()

print(f'Total are = {total_area}')

significant_area =[]

# Load and regrid model cubes
model_file = f'/data/users/samuel.sancar/Heat_extremes/Trends/Trend_significance/P-vals/'

# Read in excel sheet containing model names
df = pd.read_excel("/data/users/samuel.sancar/Model_table.xlsx", sheet_name='Model Stack', engine="openpyxl")

# ensure consistent column names (strip spaces)
df.columns = df.columns.str.strip()

significant_area =[]

for row in np.arange(0,63):
    
    row_idx = row  
    row = df.iloc[row_idx]

    # Map columns to variables 
    Mode       = str(row["Coupled/SST"]).strip()
    model      = str(row["Model"]).strip()
    experiment = str(row["Base experiment"]).strip()
    member     = str(row["Member"]).strip()

    print(model, experiment, member)

    cube = iris.load_cube(f'{model_file}{model}_{experiment}_{member}{model_selection}')
    cube, _ = lat_land_mask(cube, template, lat_minimum, lat_maximum)
    cube = units_and_coords(cube)

    # Apply significance test level
    sig = np.ma.filled(cube.data < 0.05, False)  # True only where significant and valid

    weights2d = iris.analysis.cartography.area_weights(cube)
    sig_area = np.sum(np.asarray(weights2d)[sig])

    proportion_area = (sig_area / total_area) * 100
    
    print(f"Significant area in {model}_{experiment}_{member}: {proportion_area}%")
    significant_area.append(sig_area)

# Process ERA5
era5_cube, _ = lat_land_mask(era5_cube, template, lat_minimum, lat_maximum)
era5_cube = units_and_coords(era5_cube)
    
# Apply significance test level
sig = np.ma.filled(era5_cube.data < 0.05, False)
weights2d = iris.analysis.cartography.area_weights(era5_cube)
era5_sig_area = np.sum(np.asarray(weights2d)[sig])

# Process JRA-55
jra_cube, _ = lat_land_mask(jra_cube, template, lat_minimum, lat_maximum)
jra_cube = units_and_coords(jra_cube)
    
# Apply significance test level
sig = np.ma.filled(jra_cube.data < 0.05, False)
weights2d = iris.analysis.cartography.area_weights(jra_cube)
jra_sig_area = np.sum(np.asarray(weights2d)[sig])

# Find area that both ERA5 and JRA-55 agree
era5_sig = np.ma.filled(era5_cube.data < 0.05, False)
jra_sig  = np.ma.filled(jra_cube.data  < 0.05, False)

both_sig = era5_sig & jra_sig

weights2d = iris.analysis.cartography.area_weights(era5_cube)  # same grid assumed
both_sig_area = np.sum(np.asarray(weights2d)[both_sig])


# Calculate JRA-55 and ERA5 percentile within model spread.
significant_area_array = np.array(significant_area)
era5_area_percentile = (np.sum(significant_area_array < era5_sig_area))/ len(significant_area_array)*100
jra_are_percentile =  (np.sum(significant_area_array < jra_sig_area))/ len(significant_area_array)*100
both_area_percentile = (np.sum(significant_area_array < both_sig_area))/ len(significant_area_array)*100

proportion_area = (era5_sig_area / total_area) * 100
jra_proportion_area = (jra_sig_area / total_area) * 100
model_average_sig_area = np.mean(significant_area_array)
proportion_model_average_area = (model_average_sig_area / total_area)*100

print(f"Proportion of land area that is significant for ERA5: {np.round(proportion_area, 2)}%")
print(f"Proportion of land area that is significant for JRA-55: {np.round(jra_proportion_area, 2)}%")
print(f"ERA5 significant area: {era5_sig_area}, Percentile: {era5_area_percentile:.2f}%")
print(f"JRA-55 significant area: {jra_sig_area}, Percentile: {jra_are_percentile:.2f}%")
print(f"Both significant area: {both_sig_area}, Percentile: {both_area_percentile:.2f}%")
print(f"Model average significant area: {model_average_sig_area}, Proportion of land area : {np.round(proportion_model_average_area)}%")

plt.style.use('ggplot')

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
save_figure(plt, f'HighResMIP_significant_gridpoints_histogram', 'Graphs')
plt.show()



