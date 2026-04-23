"""
23/09/2025
-----------

Histogram of FDR correction for tailwidth trends for all models, ERA5 and JRA-55 Reanlysis data.

To run: python Model_FDR_histogram.py
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

# Beginning & end year+1
t1 = 1980
t2 = 2023

significant = sys.argv[1]  # Choose from 'positive' or 'all', 

era5 = 'era5'
merra = 'MERRA'
jra = f'jra55__minmax_surf.015_tmax.reg_tl319.'
variable = 'daily_maximum_2m_temperature'

# Load template cube
template_cube = iris.load_cube('template_cube_1deg.nc')

if significant == 'positive':
    graph_name = 'positive trends'
    selection = f'_positive_fdr_corrected_p_{t1}-{t2-1}.nc'

else:
    graph_name = 'trends'
    selection = f'_fdr_corrected_p_{t1}-{t2-1}.nc'

            
# Load and regrid ERA5 trend data
datadir = '/data/users/samuel.sancar/Heat_extremes/Trends/Reanalysis_trends/'

era5_cube = iris.load_cube(f"{datadir}{era5}_{variable}{selection}")
jra_cube = iris.load_cube(f"{datadir}{jra}{selection}")

if t1 == 1980:
    merra_cube = iris.load_cube(f"{datadir}{merra}_{variable}{selection}")

# Load and regrid model cubes
model_file = f'/data/users/samuel.sancar/Heat_extremes/Trends/Trend_significance/FDR_corrections/'

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

    cube = iris.load_cube(f'{model_file}{model}_{experiment}_{member}{selection}')
    cube = units_and_coords(cube)

    # Area weights
    sig_mask = cube.data.astype(bool)
    sig_cube = cube.copy()
    sig_cube.data = sig_mask.astype(int)  
    
    area_weights = iris.analysis.cartography.area_weights(sig_cube)
    
    sig_area = sig_cube.collapsed(['latitude', 'longitude'],iris.analysis.SUM,weights=area_weights)

    significant_area.append(float(sig_area.data))
    
# Total land area
template = iris.load_cube('template_cube_1deg.nc')
benchmark_unmasked = iris.load_cube(f"/data/users/samuel.sancar/Heat_extremes/ERA5_tailwidths/era5_daily_maximum_2m_temperature_1958_tail_width.nc")
benchmark, total_gridpoints = lat_land_mask(benchmark_unmasked, template, -60, 80)
benchmark = units_and_coords(benchmark)

weights2d = iris.analysis.cartography.area_weights(benchmark)
mask = np.ma.getmaskarray(benchmark.data)
weights = np.asarray(weights2d)[~mask].astype(float)
total_area = weights.sum()

print(f'Total are = {total_area}')

# ERA5 
era5_cube = units_and_coords(era5_cube)

# Boolean significance mask
era5_sig_mask = era5_cube.data.astype(bool)
era5_sig_cube = era5_cube.copy()
era5_sig_cube.data = era5_sig_mask.astype(int)
sum_era5_sig_gridpoints = np.count_nonzero(era5_sig_mask)

era5_area_weights = iris.analysis.cartography.area_weights(era5_sig_cube)

# Area
era5_sig_area = era5_sig_cube.collapsed(['latitude', 'longitude'],iris.analysis.SUM,weights=era5_area_weights).data

# JRA-55 
jra_cube = units_and_coords(jra_cube)

# Boolean significance mask
jra_sig_mask = jra_cube.data.astype(bool)
jra_sig_cube = jra_cube.copy()
jra_sig_cube.data = jra_sig_mask.astype(int)
sum_jra_sig_gridpoints = np.count_nonzero(jra_sig_mask)

jra_area_weights = iris.analysis.cartography.area_weights(jra_sig_cube)

# Area
jra_sig_area = jra_sig_cube.collapsed(['latitude', 'longitude'],iris.analysis.SUM,weights=jra_area_weights).data

jra_sig_gridpoints = jra_cube.data 
sum_jra_sig_gridpoints = np.count_nonzero(jra_sig_gridpoints)

if t1 == 1980:
    
    print('Also calculating Merra significant area')
    # MERRA
    merra_cube = units_and_coords(merra_cube)

        # Boolean significance mask
    merra_sig_mask = merra_cube.data.astype(bool)
    merra_sig_cube = merra_cube.copy()
    merra_sig_cube.data = merra_sig_mask.astype(int)
    sum_merra_sig_gridpoints = np.count_nonzero(merra_sig_mask)

    merra_area_weights = iris.analysis.cartography.area_weights(merra_sig_cube)

    # Area
    merra_sig_area = merra_sig_cube.collapsed(['latitude', 'longitude'],iris.analysis.SUM,weights=merra_area_weights).data
    
    both_merra = era5_sig_mask & merra_sig_mask
    
    both_merra_sig_cube = merra_cube.copy()
    both_merra_sig_cube.data = both_merra.astype(int)
    sum_both_merra_sig_gridpoints = np.count_nonzero(both_merra)
    
    both_merra_area_weights = iris.analysis.cartography.area_weights(both_merra_sig_cube)
    
    both_merra_sig_area = both_merra_sig_cube.collapsed(['latitude', 'longitude'], iris.analysis.SUM, weights=both_merra_area_weights).data
    
# JRA-55 and ERA5 agreed sig points

both_mask = era5_sig_mask & jra_sig_mask

both_sig_cube = jra_cube.copy()
both_sig_cube.data = both_mask.astype(int)
sum_both_sig_gridpoints = np.count_nonzero(both_mask)

both_area_weights = iris.analysis.cartography.area_weights(both_sig_cube)

# Area
both_sig_area = both_sig_cube.collapsed(['latitude', 'longitude'],iris.analysis.SUM,weights=both_area_weights).data

sum_both_sig_gridpoints = np.count_nonzero(both_sig_cube)

# Calculate JRA-55 and ERA5 percentile within model spread.
significant_area_array = np.array(significant_area)

era5_area_percentile = (np.sum(significant_area_array < era5_sig_area))/ len(significant_area_array)*100
jra_area_percentile =  (np.sum(significant_area_array < jra_sig_area))/ len(significant_area_array)*100
both_area_percentile = (np.sum(significant_area_array < both_sig_area))/ len(significant_area_array)*100

if t1 == 1980:
    
    merra_area_percentile = (np.sum(significant_area_array < merra_sig_area))/ len(significant_area_array)*100
    both_merra_area_percentile = (np.sum(significant_area_array < both_merra_sig_area))/ len(significant_area_array)*100
    merra_proportion = merra_sig_area/total_area * 100

proportion = era5_sig_area/total_area * 100
jra_proportion = jra_sig_area/total_area * 100
difference = era5_sig_area - jra_sig_area

print(f"ERA5 significant area: {era5_sig_area}, Percentile: {era5_area_percentile:.2f}%")
print(f"Proportion of valid land that is significant for ERA5: {np.round(proportion, 2)}%")

print(f"JRA-55 significant area: {jra_sig_area}, Percentile: {jra_area_percentile:.2f}%")
print(f"Proportion of valid land that is significant for JRA-55: {np.round(jra_proportion, 2)}%")

print(f"Difference in significant area between JRA-55 and ERA5 {difference}")
print(f"Both significant area: {both_sig_area}, Percentile: {both_area_percentile:.2f}%")

if t1 ==1980:
    print(f"Proportion of valid land that is significant for MERRA-2: {np.round(merra_proportion, 2)}%")
    print(f"MERRA significant area: {merra_sig_area}, Percentile: {merra_area_percentile:.2f}%")
    print(f"MERRA-2 and ERA5 overlap area: {both_merra_sig_area}, Percentile: {both_merra_area_percentile:.2f}%")

# Histogram for significant area
plt.figure(figsize=(12, 5))
counts, bins, patches = plt.hist(
    significant_area,
    bins=10,
    range=(significant_area_array.min(), significant_area_array.max()),
    histtype='bar',
    color='#4682B4', 
    alpha=0.85,
    edgecolor='white',
    linewidth=2,
    label = 'HighResMIP'
)

plt.axvline(x=era5_sig_area, color='black', linestyle='-', linewidth=2.5, label='ERA5')

# Label histogram at JRA-55 value
plt.axvline(x=jra_sig_area, color='black', linestyle='--', linewidth=2.5, label='JRA-55')

if t1 == 1980:
    
    plt.axvline(x=merra_sig_area, color='#1B998B', linestyle='--', linewidth=2.5, label='MERRA')
    plt.axvline(x=both_merra_sig_area, color='#1B998B', linestyle='-.', linewidth=2.5, label='ERA5 MERRA agreement')
    
# Label histogram at JRA-55 value
plt.axvline(x=both_sig_area, color="black", linestyle='-.', linewidth=2.5)
plt.text(
        both_sig_area + 50, max(counts) * 0.82, 'ERA5 & JRA-55 agreement',
        color='black', fontsize=15, fontweight='bold',
        bbox=dict(facecolor='white', edgecolor='black', boxstyle='round' ,pad=0.3, alpha=0.8)
    )

plt.xlabel('Area of significant trends (\u00D710\u00B9\u00B2 m\u00B2)', fontsize=18, fontweight='bold')
plt.ylabel('Frequency', fontsize=18, fontweight='bold')
plt.title(f'HRMIP area distribution of significant {graph_name} ', fontsize=20, fontweight='bold')

plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.grid(True, linestyle='--', alpha=0.5, linewidth=1.2)

plt.legend(fontsize=13, loc='upper right')
plt.tight_layout()
save_figure(plt, f'HighResMIP_fdr_{graph_name}_histogram', 'Graphs')
plt.show()



