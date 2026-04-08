"""
17/07/25
---------
Plot histogram in tail widths trends of ERA5 data.
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
t1 = int(sys.argv[1])
t2 = int(sys.argv[2])

reanalysis = 'era5'
variable = 'daily_maximum_2m_temperature'
jra = 'jra55'
variable_jra = '_minmax_surf.015_tmax.reg_tl319.'

# Load ERA5 trend data
datadir = '/data/users/samuel.sancar/Heat_extremes/Trends/Reanalysis_trends/'
trend = iris.load_cube(f"{datadir}{reanalysis}_{variable}_tailwidth_trends_{t1}-{t2}.nc")
era5_trend, _ = lat_land_mask(trend, -60,80)
era5_data = era5_trend.data.astype(float).flatten()
valid_era5 = ~np.isnan(era5_data)
era5_data = era5_data[valid_era5]
era5_pos = era5_data[era5_data >= 0]
era5_neg  = era5_data[era5_data <= 0]
max_value = era5_data.max()
    
era5_masked = units_and_coords(era5_trend)
era5_area_weights_2d = iris.analysis.cartography.area_weights(era5_masked)
era5_area_weights = era5_area_weights_2d.flatten()[valid_era5]
era5_pos_w  = era5_area_weights[era5_data >= 0]
era5_neg_w = era5_area_weights[era5_data <= 0]

# Load JRA-55 trend data
datadir_jra = '/data/users/samuel.sancar/Heat_extremes/Trends/Reanalysis_trends/'
trend_jra = iris.load_cube(f"{datadir_jra}{jra}_{variable_jra}_tailwidth_trends_{t1}-{t2}.nc")
jra_trend, _ = lat_land_mask(trend_jra, -60,80)
jra_data = jra_trend.data.astype(float).flatten()
valid_jra = ~np.isnan(jra_data)
jra_data = jra_data[valid_jra]
jra_pos = jra_data[jra_data >= 0]
jra_neg  = jra_data[jra_data < 0]
max_value = jra_data.max()
    
jra_masked = units_and_coords(jra_trend)
jra_area_weights_2d = iris.analysis.cartography.area_weights(jra_masked)
jra_area_weights = jra_area_weights_2d.flatten()[valid_jra]
jra_pos_w  = jra_area_weights[jra_data >= 0]
jra_neg_w = jra_area_weights[jra_data < 0]

coupled_files = sorted(glob.glob(f'/data/users/samuel.sancar/Heat_extremes/Trends/Coupled_model_trends/*_tailwidth_trend_{t1}-2022.nc'))
sst_files = sorted(glob.glob(f'/data/users/samuel.sancar/Heat_extremes/Trends/SST_model_trends/*_tailwidth_trend_{t1}-2022.nc'))
coupled_data = [iris.load_cube(f) for f in coupled_files]
sst_data = [iris.load_cube(f) for f in sst_files]

bins = np.arange(0, 1, 0.025)

#zero = 0.0
#bins_neg = np.concatenate((bins[bins < zero], [zero]))     
#bins = np.concatenate(([zero], bins[bins > zero]))   
 
# SST and Coupled models 
plt.figure(figsize=(10, 10))

ntot = 0

for cube in sst_data:
    
    ntot = ntot+1
    
    masked_model, _ = lat_land_mask(cube, -60,80)
    model_data = masked_model.data.astype(float).flatten()
    valid = ~np.isnan(model_data)
    model_data = model_data[valid]
    pos = model_data[model_data > 0]
    neg  = model_data[model_data < 0]
    max_value = model_data.max()
    
    masked = units_and_coords(masked_model)
    area_weights_2d = iris.analysis.cartography.area_weights(masked)
    area_weights = area_weights_2d.flatten()[valid]
    pos_w  = area_weights[model_data > 0]
    neg_w = area_weights[model_data < 0]
    plt.hist(model_data, bins=bins, weights=(area_weights/10**13), density= False, cumulative=-1, histtype='step', alpha =0.5, edgecolor='red')
    #plt.hist(neg, bins=bins_neg, weights=neg_w, density= True, cumulative = True, histtype='stepfilled', alpha=0.2, color='skyblue', edgecolor='black')
    
    if ntot == 1: # To ensure correct labelling
        
        plt.hist(model_data, bins=bins, weights=(area_weights/10**13), density= False, cumulative=-1, histtype='step', alpha =0.5, edgecolor='red', label='HighResMIP members')
        
    else:
        plt.hist(model_data, bins=bins, weights=(area_weights/10**13), density= False, cumulative=-1, histtype='step', alpha =0.5, edgecolor='red')

for cube in coupled_data:
    masked_model, _ = lat_land_mask(cube, -60,80)
    model_data = masked_model.data.astype(float).flatten()
    valid = ~np.isnan(model_data)
    model_data = model_data[valid]
    pos = model_data[model_data > 0]
    neg  = model_data[model_data < 0]
    max_value = model_data.max()
    
    masked = units_and_coords(masked_model)
    area_weights_2d = iris.analysis.cartography.area_weights(masked)
    area_weights = area_weights_2d.flatten()[valid]
    pos_w  = area_weights[model_data > 0]
    neg_w = area_weights[model_data < 0]
    
    plt.hist(model_data, bins=bins, weights=(area_weights/10**13), density= False, cumulative=-1, histtype='step', alpha =0.5, edgecolor='red')
    #plt.hist(neg, bins=bins_neg, weights=neg_w, density= True, cumulative = True, histtype='stepfilled', alpha=0.2, color='skyblue', edgecolor='black')

plt.hist(era5_data, weights=(era5_area_weights/10**13), density= False, cumulative =-1, bins=bins, histtype='step', color='black', linewidth=3.5, label='ERA5') #, linestyle='-'
#plt.hist(era5_neg, weights=era5_neg_w, density= True, cumulative = True, bins=bins_neg, range = (-0.5, 0), histtype='step', color='red', linewidth=2.5) #, linestyle='-'
#plt.hist(jra_neg, weights=jra_neg_w, density= True, cumulative = True, bins=bins_neg, range = (-0.5, 0), histtype='step', color='blue', linewidth=2.5, linestyle='-.')
plt.hist(jra_data, weights=(jra_area_weights/10**13), density= False, cumulative =-1, bins=bins, histtype='step', color='black', linewidth=3.5, linestyle='-.', label='JRA-55')
plt.xlabel('Tail width trend (°C per Decade)', fontsize=16, fontweight='bold')
plt.xticks(fontsize = 15)
plt.ylabel(f'Cumulative area (\u00D710\u00B9\u00B3 m\u00B2)', fontsize=16, fontweight='bold')
plt.yticks(fontsize = 15)

plt.xlim(0, 1)

plt.title('Tail width trend cumulative area distributions', fontsize=17, fontweight='bold')
plt.legend(loc='upper right', fontsize=15)
plt.grid(True, linestyle='--', alpha=0.6)
plt.tight_layout()
save_figure(plt, 'Tail_width_trend_cumulative_distributions', 'Graphs')
plt.show()