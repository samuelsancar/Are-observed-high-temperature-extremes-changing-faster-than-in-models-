"""
04/11/25
---------
Plot histogram in tail widths trends of ERA5 data.
To run: python plot_tailwidth_trend.py [1958] [2023]

"""

import iris
import sys
import matplotlib.pyplot as plt
import numpy as np
import glob
import os
from Lat_Land_mask_function import lat_land_mask
from Save_figure_function import save_figure
import cf_units


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

# Read command-line arguments
reanalysis = 'era5'
t1 = int(sys.argv[1])
t2 = int(sys.argv[2])
variable = 'daily_maximum_2m_temperature'

jra = 'jra55'
variable_jra = '_minmax_surf.015_tmax.reg_tl319.'

# Load ERA5 trend data
datadir = '/data/users/samuel.sancar/Heat_extremes/Trends/Reanalysis_trends/'
trend = iris.load_cube(f"{datadir}{reanalysis}_{variable}_tailwidth_trends_{t1}-{t2}.nc")
era5_trend, _ = lat_land_mask(trend, -60, 70)
era5_data = era5_trend.data.astype(float).flatten()
valid_era5 = ~np.isnan(era5_data)
era5_data = era5_data[valid_era5]
max_value = era5_data.max()
    
era5_masked = units_and_coords(era5_trend)
era5_area_weights_2d = iris.analysis.cartography.area_weights(era5_masked)
era5_area_weights = era5_area_weights_2d.flatten()[valid_era5]
   
# Load JRA-55 trend data
datadir_jra = '/data/users/samuel.sancar/Heat_extremes/Trends/Reanalysis_trends/'
trend_jra = iris.load_cube(f"{datadir_jra}{jra}_{variable_jra}_tailwidth_trends_{t1}-{t2}.nc")
jra_trend, _ = lat_land_mask(trend_jra, -60, 70)
jra_data = jra_trend.data.astype(float).flatten()
valid_jra = ~np.isnan(jra_data)
jra_data = jra_data[valid_jra]
max_value = jra_data.max()
    
jra_masked = units_and_coords(jra_trend)
jra_area_weights_2d = iris.analysis.cartography.area_weights(jra_masked)
jra_area_weights = jra_area_weights_2d.flatten()[valid_jra]

# Load model data
CANARI_files = sorted(glob.glob('/data/users/samuel.sancar/Heat_extremes/Trends/CANARI_trends/*_tailwidth_trend.nc'))
CANARI_data = [(iris.load_cube(f), f) for f in CANARI_files]

# Define bins
bins = np.arange(-0.5, 1.5, 0.1)

plt.figure(figsize=(14, 12))

for cube, f in CANARI_data:
    
    masked_model, _ = lat_land_mask(cube, -60, 70)
    model_data = masked_model.data.astype(float).flatten()
    valid = ~np.isnan(model_data)
    model_data = model_data[valid]
    max_value = model_data.max()
    
    masked = units_and_coords(masked_model)
    area_weights_2d = iris.analysis.cartography.area_weights(masked)
    area_weights = area_weights_2d.flatten()[valid]
   
    fname = os.path.basename(f)
    print(f'{fname} has a max tail width trend of {max_value} degrees C/decade.')
    plt.hist(model_data, bins=bins, weights=area_weights, histtype='stepfilled', alpha=0.2, color='skyblue', edgecolor='black')

plt.hist(era5_data, bins=bins, weights = era5_area_weights, histtype='step', color='red', linewidth=2.5, linestyle='--', label='ERA5')
plt.hist(jra_data, bins=bins, weights = jra_area_weights, histtype='step', color='blue', linewidth=2.5, linestyle='-.', label='JRA-55')
plt.xlabel('Trend Magnitude (°C per Decade)', fontsize=10, fontweight='bold')
plt.ylabel('Area (m^2)', fontsize=10, fontweight='bold')
plt.title(f'Histogram of Tail-Width Trends (CANARI LE, ERA5 & JRA-55, {t1}–{t2})', fontsize=15, fontweight='bold')
plt.legend(loc='upper right', fontsize=10)
plt.xlim(-0.5, 1.5)
plt.grid(True, linestyle='--', alpha=0.6)
plt.tight_layout()
plt.show()

