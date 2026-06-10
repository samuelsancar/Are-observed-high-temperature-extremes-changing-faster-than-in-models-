"""
21/05/26
---------

Plot probability density histogram in tail widths trends of CANARI and reanalysis data.
To run: python plot_tailwidth_trend.py 

To run: python CANARI_trends_histogram.py
"""

import iris
import sys
import matplotlib.pyplot as plt
import numpy as np
import glob
import os
from Lat_Land_mask_function import lat_land_mask
import cf_units


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

def extract_data_and_weights(cube):
    """
    Extract data and area weights from a masked 2D cube,
    respecting the mask.
    """
    cube = units_and_coords(cube)
    weights2d = iris.analysis.cartography.area_weights(cube)

    mask = np.ma.getmaskarray(cube.data)

    data = np.asarray(cube.data)[~mask].astype(float)
    weights = np.asarray(weights2d)[~mask].astype(float)

    return data, weights

# Beginning & end year+1
t1 = 1958
t2 = 2023

lat_minimum = -60
lat_maximum = 60

reanalysis = 'era5'
variable = 'daily_maximum_2m_temperature'
jra = 'jra55'
variable_jra = '_minmax_surf.015_tmax.reg_tl319.'
Berkeley = 'BerkeleyEarth'

if t1 == 1980:
    merra = 'MERRA'

template = iris.load_cube('template_cube_1deg.nc')

# Load ERA5 trend data
datadir = '/data/users/samuel.sancar/Heat_extremes/Trends/Reanalysis_trends/'
trend = iris.load_cube(f"{datadir}{reanalysis}_{variable}_tailwidth_trends_{t1}-{t2-1}.nc")
era5_trend, _ = lat_land_mask(trend, template, lat_minimum, lat_maximum)

era5_data, era5_weights = extract_data_and_weights(era5_trend)

# Load JRA-55 trend data
datadir_jra = '/data/users/samuel.sancar/Heat_extremes/Trends/Reanalysis_trends/'
trend_jra = iris.load_cube(f"{datadir_jra}{jra}_{variable_jra}_tailwidth_trends_{t1}-{t2-1}.nc")
jra_trend, _ = lat_land_mask(trend_jra, template, lat_minimum, lat_maximum)

jra_data, jra_weights = extract_data_and_weights(jra_trend)

# Load BerkeleyEarth trend data
trend_berkeley = iris.load_cube(f"{datadir}{Berkeley}_{variable}_tailwidth_trends_{t1}-{t2-1}.nc")
berkeley_trend, _ = lat_land_mask(trend_berkeley, template, lat_minimum, lat_maximum)

berkeley_data, berkeley_weights = extract_data_and_weights(berkeley_trend)

# Load MERRA2 trend data
if t1 == 1980:
    datadir_merra = '/data/users/samuel.sancar/Heat_extremes/Trends/Reanalysis_trends/'
    merra_trend = iris.load_cube(f"{datadir_merra}{merra}_{variable}_tailwidth_trends_{t1}-{t2-1}.nc")
    merra_trend, _ = lat_land_mask(merra_trend, template, lat_minimum, lat_maximum)

    merra_data, merra_weights = extract_data_and_weights(merra_trend)

bins = np.arange(-0.9, 1.2, 0.05)

ntot = 0

if t1==1980:
    x_min = -1.0
else:
    x_min = -0.9

plt.figure(figsize=(14, 12))

for i in np.arange(1,41):
    
    ntot = ntot+1
    
    cube = iris.load_cube(f'/data/users/samuel.sancar/Heat_extremes/Trends/CANARI_trends/{i}_tailwidth_trend_{t1}-{t2-1}.nc')
    
    masked_model, _ = lat_land_mask(cube, template, lat_minimum, lat_maximum)

    model_data, model_weights = extract_data_and_weights(masked_model)

    if ntot == 1:
        plt.hist(model_data,bins=bins,weights=model_weights,density=True,alpha = 0.5,range = (x_min, 1.2),histtype='step',label = 'HRMIP members',edgecolor='red')
    else:
        plt.hist(model_data,bins=bins,weights=model_weights,density=True,alpha = 0.5,range = (x_min, 1.2),histtype='step',edgecolor='red')
    
plt.hist(era5_data, weights=era5_weights, bins=bins, range = (x_min, 1.2), density = True, histtype='step', color='black', linewidth=3, label='ERA5') #, linestyle='-'
plt.hist(jra_data, weights=jra_weights, bins=bins, range = (x_min, 1.2), density = True, histtype='step', color='black', linewidth=3, linestyle='-.', label='JRA-55')
plt.hist(berkeley_data, weights=berkeley_weights, bins=bins, range = (x_min, 1.2), density = True, histtype='step', color='blue', linewidth=3, label='Berkeley') #, linestyle='-'

if t1 == 1980:
    plt.hist(merra_data, weights = merra_weights, bins=bins, range = (x_min, 1.2), density = True, histtype='step', color= "blue", linewidth=3, linestyle = '-.', label='MERRA-2')

plt.xlabel('Tail width trend (°C per Decade)', fontsize=17, fontweight='bold')
plt.xticks(fontsize = 16)
plt.xlim(x_min, 1.2)

plt.ylabel(f'Density', fontsize=17, fontweight='bold')
plt.yticks(fontsize = 16)
plt.yscale('log')

plt.title('a) Tail width trend probability distributions (CANARI)', fontsize=19, fontweight='bold')
plt.legend(loc='upper right', fontsize=17)
plt.grid(True, linestyle='--', alpha=0.6)
plt.tight_layout()
plt.show()
