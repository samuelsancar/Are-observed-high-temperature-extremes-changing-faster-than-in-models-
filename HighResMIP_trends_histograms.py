"""
17/07/25
---------
Plot histogram in tail widths trends of HighResMIP, JRA-55, and ERA5 data.

Produces figure 2a.

"""

import iris
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
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
lat_maximum = 80

reanalysis = 'era5'
variable = 'daily_maximum_2m_temperature'
jra = 'jra55'
variable_jra = '_minmax_surf.015_tmax.reg_tl319.'

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

# Load MERRA2 trend data
if t1 == 1980:
    datadir_merra = '/data/users/samuel.sancar/Heat_extremes/Trends/Reanalysis_trends/'
    merra_trend = iris.load_cube(f"{datadir_merra}{merra}_{variable}_tailwidth_trends_{t1}-{t2-1}.nc")
    merra_trend, _ = lat_land_mask(merra_trend, template, lat_minimum, lat_maximum)

    merra_data, merra_weights = extract_data_and_weights(merra_trend)

#Read in excel sheet containing model names
df = pd.read_excel("/data/users/samuel.sancar/Model_table.xlsx", sheet_name='Model Stack', engine="openpyxl")

# ensure consistent column names (strip spaces)
df.columns = df.columns.str.strip()

bins = np.arange(-0.9, 1.2, 0.05)

ntot = 0

if t1==1980:
    x_min = -1.0
else:
    x_min = -0.9
    
# SST and Coupled models 
plt.figure(figsize=(10, 8))
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
    
    ntot = ntot+1
    
    masked_model, _ = lat_land_mask(cube, template, lat_minimum, lat_maximum)

    model_data, model_weights = extract_data_and_weights(masked_model)

    if ntot == 1:
        plt.hist(model_data,bins=bins,weights=model_weights,density=True,alpha = 0.5,range = (x_min, 1.2),histtype='step',label = 'HRMIP members',edgecolor='red')
    else:
        plt.hist(model_data,bins=bins,weights=model_weights,density=True,alpha = 0.5,range = (x_min, 1.2),histtype='step',edgecolor='red')
    
ax = plt.gca()
for spine in ax.spines.values():
    spine.set_linewidth(2)

plt.hist(era5_data, weights=era5_weights, bins=bins, range = (x_min, 1.2), density = True, histtype='step', color='black', linewidth=3, label='ERA5') #, linestyle='-'
plt.hist(jra_data, weights=jra_weights, bins=bins, range = (x_min, 1.2), density = True, histtype='step', color='black', linewidth=3, linestyle='-.', label='JRA-55')

if t1 == 1980:
    plt.hist(merra_data, weights = merra_weights, bins=bins, range = (x_min, 1.2), density = True, histtype='step', color= "blue", linewidth=3, linestyle = '-.', label='MERRA-2')

plt.xlabel('Tail width trend (°C per Decade)', fontsize=17, fontweight='bold')
plt.xticks(fontsize = 16)
plt.xlim(x_min, 1.2)

plt.ylabel(f'Density', fontsize=17, fontweight='bold')
plt.yticks(fontsize = 16)
plt.yscale('log')

plt.title('a) Tail width trend probability distributions', fontsize=19, fontweight='bold')
plt.legend(loc='upper right', fontsize=17)
plt.grid(True, linestyle='--', alpha=0.6)
plt.tight_layout()
save_figure(plt, f'Tail_width_trend_densities({t1}-{t2-1})', 'Graphs')
plt.show()

#NOT USED
print('Next section NOT USED!! Only kept if needed for future reference')
# In each bin, what percentile does ERA5 lie within the HighResMIP spread
era5_percentile = []
jra_percentile = []

for i in enumerate(bins):
    
    model_areas = []
    
    for cube in coupled_data:
        
        masked_model, _ = lat_land_mask(cube, template, lat_minimum, lat_maximum)
    
        sig_area = area_in_bin(masked_model, i[1])
        model_areas.append(sig_area)
        
    for cube in sst_data:
        
        masked_model, _ = lat_land_mask(cube, template, lat_minimum, lat_maximum)
        
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



