"""
17/07/25
---------

Plot histogram in tail widths trends of ERA5 data.
Produces figure 2b.

To run: python HighResMIP_cumulative_trends.py
"""

import iris
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import cf_units
from Lat_Land_mask_function import lat_land_mask
from Save_figure_function import save_figure

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
merra = 'MERRA'
Berkeley = 'BerkeleyEarth'

template = iris.load_cube('template_cube_1deg.nc')
datadir = '/data/users/samuel.sancar/Heat_extremes/Trends/Reanalysis_trends/'

# Load ERA5 trend data
trend = iris.load_cube(f"{datadir}{reanalysis}_{variable}_tailwidth_trends_{t1}-{t2-1}.nc")
era5_trend, _ = lat_land_mask(trend, template, lat_minimum, lat_maximum)

era5_data, era5_weights = extract_data_and_weights(era5_trend)

# Load JRA-55 trend data
trend_jra = iris.load_cube(f"{datadir}{jra}_{variable_jra}_tailwidth_trends_{t1}-{t2-1}.nc")
jra_trend, _ = lat_land_mask(trend_jra, template, lat_minimum, lat_maximum)

jra_data, jra_weights = extract_data_and_weights(jra_trend)

# Load BerkeleyEarth trend data
trend_berkeley = iris.load_cube(f"{datadir}{Berkeley}_{variable}_tailwidth_trends_{t1}-{t2-1}.nc")
berkeley_trend, _ = lat_land_mask(trend_berkeley, template, lat_minimum, lat_maximum)

berkeley_data, berkeley_weights = extract_data_and_weights(berkeley_trend)

# Load MERRA2 trend data
if t1 == 1980:
    merra_trend = iris.load_cube(f"{datadir}{merra}_{variable}_tailwidth_trends_{t1}-{t2-1}.nc")
    merra_trend, _ = lat_land_mask(merra_trend, template, lat_minimum, lat_maximum)

    merra_data, merra_weights = extract_data_and_weights(merra_trend)

# Read in excel sheet containing model names
df = pd.read_excel("/data/users/samuel.sancar/Model_table.xlsx", sheet_name='Model Stack', engine="openpyxl")

# ensure consistent column names (strip spaces)
df.columns = df.columns.str.strip()

bins = np.arange(0, 1.2, 0.025)
 
# SST and Coupled models 
plt.figure(figsize=(10, 8))

coupled = 0
sst = 0

coupled_tot = 0

for row in np.arange(0,63):
    
    row_idx = row  
    row = df.iloc[row_idx]

    # Map columns to variables 
    Mode       = str(row["Coupled/SST"]).strip()
    model      = str(row["Model"]).strip()
    experiment = str(row["Base experiment"]).strip()
    member     = str(row["Member"]).strip()

    print(model, experiment, member)
    
    datadir = f'/data/users/samuel.sancar/Heat_extremes/Trends/{Mode}_model_trends/'
    cube  = iris.load_cube(f'{datadir}{model}_{experiment}_{member}_tailwidth_trend_{t1}-{t2-1}.nc')
        
    masked_model, _ = lat_land_mask(cube, template, lat_minimum, lat_maximum)

    model_data, model_weights = extract_data_and_weights(masked_model)
    
    #coupled_tot = 0

    #if Mode == 'Coupled':
    #    color = 'blue'
    #    label = 'Coupled'
    #    coupled = coupled+1
    #else:
    #    color = 'red'
    #    label = 'SST'
    #    sst = sst+1
    
    if coupled == 1: # To ensure only one label
        plt.hist(model_data, bins=bins, weights= model_weights, density= False, cumulative=-1, histtype='step', alpha = 0.5, edgecolor='red', label='HighResMIP')#/10**13
        coupled_tot = coupled_tot + 1 
#   elif sst == 1:
#        plt.hist(model_data, bins=bins, weights= model_weights, density= False, cumulative=-1, histtype='step', alpha = 1, edgecolor=color, label=label)#/10**13
    else:
        plt.hist(model_data, bins=bins, weights= model_weights, density= False, cumulative=-1, histtype='step', alpha = 0.5, edgecolor='red')#/10**13

ax = plt.gca()

for spine in ax.spines.values():
    spine.set_linewidth(2)
    
plt.hist(era5_data, weights=era5_weights, density= False, cumulative =-1, bins=bins, histtype='step', color='black', linewidth=3, label='ERA5') #, linestyle='-'
plt.hist(jra_data, weights=jra_weights, density= False, cumulative =-1, bins=bins, histtype='step', color='black', linewidth=3, linestyle='-.', label='JRA-55')#/10**13
plt.hist(berkeley_data, weights=berkeley_weights, density= False, cumulative =-1, bins=bins, histtype='step', color='blue', linewidth=3, label='Berkeley') #, linestyle='-'

if t1 == 1980:
    plt.hist(merra_data, weights = merra_weights, density = False, cumulative = -1, bins=bins, histtype='step', color= "black", linewidth=3, linestyle = '-.', label='MERRA-2')

plt.xlabel('Tail width trend (°C per Decade)', fontsize=18, fontweight='bold')
plt.xticks(fontsize = 16)
plt.xlim(0, 1.2)

plt.ylabel(f'Cumulative area (m\u00B2)', fontsize=18, fontweight='bold') #\u00D710\u00B9\u00B3
plt.yscale('log')
plt.yticks(fontsize = 16)

plt.title('b) Tail width trend cumulative area distributions', fontsize=20, fontweight='bold')
plt.legend(loc='upper right', fontsize=18)
plt.grid(True, linestyle='--', alpha=0.6)
plt.tight_layout()
save_figure(plt, f'Tail_width_trend_cumulative_distributions_CoupledvsSST({t1}-{t2-1})', 'Graphs')
plt.show()


