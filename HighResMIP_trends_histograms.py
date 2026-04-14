"""
17/07/25
---------
Plot histogram in tail widths trends of HighResMIP, JRA-55, and ERA5 data.
To run: python plot_tailwidth_trend.py [1958] [2023]

"""

import iris
import sys
import matplotlib.pyplot as plt
import numpy as np
import glob
import cf_units
from Lat_Land_mask_function import lat_land_mask
from Save_figure_function import save_figure

def area_in_bin(data, i):
    
    sig_mask = units_and_coords(data)
    
    area_weights = iris.analysis.cartography.area_weights(data)
    sig_mask = (data.data >=i) & (data.data < i + 0.05)

    sig_area = np.ma.array(area_weights, mask=~sig_mask).sum()
    return sig_area

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

# Command-line arguments
t1 = 1958
t2 = 2023

reanalysis = 'era5'
variable = 'daily_maximum_2m_temperature'
jra = 'jra55'
variable_jra = '_minmax_surf.015_tmax.reg_tl319.'

# Load ERA5 trend data
datadir = '/data/users/samuel.sancar/Heat_extremes/Trends/Reanalysis_trends/'
trend = iris.load_cube(f"{datadir}{reanalysis}_{variable}_tailwidth_trends_{t1}-{t2}.nc")
era5_trend, _ = lat_land_mask(trend, -60, 80)
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
jra_trend, _ = lat_land_mask(trend_jra, -60, 80)
jra_data = jra_trend.data.astype(float).flatten()
valid_jra = ~np.isnan(jra_data)
jra_data = jra_data[valid_jra]
max_value = jra_data.max()
    
jra_masked = units_and_coords(jra_trend)
jra_area_weights_2d = iris.analysis.cartography.area_weights(jra_masked)
jra_area_weights = jra_area_weights_2d.flatten()[valid_jra]

coupled_files = sorted(glob.glob(f'/data/users/samuel.sancar/Heat_extremes/Trends/Coupled_model_trends/*_tailwidth_trend_{t1}-2022.nc'))
sst_files = sorted(glob.glob(f'/data/users/samuel.sancar/Heat_extremes/Trends/SST_model_trends/*_tailwidth_trend_{t1}-2022.nc'))
coupled_data = [iris.load_cube(f) for f in coupled_files]
sst_data = [iris.load_cube(f) for f in sst_files]

bins = np.arange(-0.9, 1.2, 0.05)

ntot = 0

# SST and Coupled models 
plt.figure(figsize=(10, 8))
for cube in sst_data:
    
    ntot = ntot+1
    
    masked_model, _ = lat_land_mask(cube, -60, 80)
    model_data = masked_model.data.astype(float).flatten()
    valid = ~np.isnan(model_data)
    model_data = model_data[valid]
    max_value = model_data.max()
    
    masked = units_and_coords(masked_model)
    area_weights_2d = iris.analysis.cartography.area_weights(masked)
    area_weights = area_weights_2d.flatten()[valid]

    if ntot == 1:
        plt.hist(
            model_data,
            bins=bins,
            weights=area_weights,
            density=True,
            alpha = 0.5,
            range = (-0.9, 1.2),
            histtype='step',
            label = 'HRMIP members',
            edgecolor='red'
        )
    else:
        plt.hist(
            model_data,
            bins=bins,
            weights=area_weights,
            density=True,
            alpha = 0.5,
            range = (-0.9, 1.2),
            histtype='step',
            edgecolor='red'
        )

for cube in coupled_data:
    masked_model, _ = lat_land_mask(cube, -60, 80)
    model_data = masked_model.data.astype(float).flatten()
    valid = ~np.isnan(model_data)
    model_data = model_data[valid]
    max_value = model_data.max()
    
    masked = units_and_coords(masked_model)
    area_weights_2d = iris.analysis.cartography.area_weights(masked)
    area_weights = area_weights_2d.flatten()[valid]

    plt.hist(
        model_data,
        bins=bins,
        weights=area_weights,
        density=True,
        alpha = 0.5,
        range = (-0.9, 1.2),
        histtype='step',
        edgecolor='red'
    )
    
ax = plt.gca()
for spine in ax.spines.values():
    spine.set_linewidth(2)

plt.hist(era5_data, weights=era5_area_weights, bins=bins, range = (-0.9, 1.2), density = True, histtype='step', color='black', linewidth=3, label='ERA5') #, linestyle='-'
plt.hist(jra_data, weights=jra_area_weights, bins=bins, range = (-0.9, 1.2), density = True, histtype='step', color='black', linewidth=3, linestyle='-.', label='JRA-55')
plt.xlabel('Tail width trend (°C per Decade)', fontsize=17, fontweight='bold')
plt.xticks(fontsize = 16)
plt.ylabel(f'Density', fontsize=17, fontweight='bold')
plt.yticks(fontsize = 16)
plt.yscale('log')
plt.xlim(-0.9, 1.2)
plt.title('Tail width trend distribution', fontsize=19, fontweight='bold')
plt.legend(loc='upper right', fontsize=17)
plt.grid(True, linestyle='--', alpha=0.6)
plt.tight_layout()
save_figure(plt, 'Tail_width_trend_densities', 'Graphs')
plt.show()

era5_percentile = []
jra_percentile = []

for i in enumerate(bins):
    
    model_areas = []
    
    for cube in coupled_data:
        
        masked_model, _ = lat_land_mask(cube, -60, 80)
    
        sig_area = area_in_bin(masked_model, i[1])
        model_areas.append(sig_area)
        
    for cube in sst_data:
        
        masked_model, _ = lat_land_mask(cube, -60, 80)
        
        sig_area = area_in_bin(masked_model, i[1])
        model_areas.append(sig_area)
        
    model_area_array = np.array(model_areas)
    
    era5_sig_area = area_in_bin(era5_trend, i[1])
    jra_sig_area = area_in_bin(jra_trend, i[1])
    
    era5 = (np.sum(model_area_array <= era5_sig_area))/ len(model_area_array)*100
    era5_percentile.append(era5)
    
    jra = (np.sum(model_area_array <= jra_sig_area)) / len(model_area_array) * 100
    jra_percentile.append(jra)
    
    print(f"For bin {i[1]:.2f} - {i[1]+0.05:.2f}: ERA5 Percentile = {era5:.2f}%, JRA-55 Percentile  = {jra:.2f}%")
    
plt.figure(figsize=(10, 10))
plt.plot(bins, era5_percentile, label='ERA5 Percentile', color='red', marker='o')
plt.plot(bins, jra_percentile, label='JRA-55 Percentile', color = 'blue', marker='o')
plt.xlabel('Trend magnitude (°C per Decade)', fontsize=15, fontweight='bold')
plt.ylabel("Percentile within mode distribution (%)", fontsize = 15, fontweight='bold')
plt.title("Percentile of ERA5 and JRA-55 trends within model distribution", fontsize=15, fontweight='bold')
plt.legend()
plt.tight_layout()
save_figure(plt, 'ERA5_JRA-55_trend_percentiles', 'Graphs')
plt.show()



