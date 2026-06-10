''''
10/03/2026
----------

Plot location of HighResMIP tail width trends above a certain threshold, and compare with ERA5.

To run: python HighResMIP_tailwidth_threshold.py
'''

import iris 
import numpy as np 
import pandas as pd
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import matplotlib.cm as mpl_cm
from Save_figure_function import save_figure
import matplotlib.colors as mcolors
import cartopy.feature as cfeature

from Lat_Land_mask_function import lat_land_mask

# Beginning & end years+1
t1 = 1958
t2 = 2023

lat_minimum = -60
lat_maximum = 90

template = iris.load_cube('template_cube_1deg.nc')

df = pd.read_excel("/data/users/samuel.sancar/Model_table.xlsx", sheet_name='Model Stack', engine="openpyxl")

df.columns = df.columns.str.strip()


thresholds = np.linspace(0.1, 0.35, 6)

for y in thresholds:
    
    threshold_count = []
    
    for row in np.arange(0,63):
        
        row_idx = row
        row = df.iloc[row_idx]
        
        # Map columns to variables 
        Mode = str(row['Coupled/SST']).strip()
        model = str(row['Model']).strip()
        experiment = str(row['Base experiment']).strip() 
        member = str(row['Member']).strip()
        
        trend = iris.load_cube(f'/data/users/samuel.sancar/Heat_extremes/Trends/{Mode}_model_trends/{model}_{experiment}_{member}_tailwidth_trend_{t1}-{t2-1}.nc')
        p87p5 = iris.load_cube(f'/data/users/samuel.sancar/Heat_extremes/Trends/{Mode}_model_trends/{model}_{experiment}_{member}_p87p5_trend_{t1}-{t2-1}.nc')
        
        trend_masked, _ = lat_land_mask(trend, template, lat_minimum, lat_maximum)
        p87p5_masked, _ = lat_land_mask(p87p5, template, lat_minimum, lat_maximum)
        
        #Threshold tail width trend:
        threshold = np.ma.filled(trend_masked.data > y, False)
        threshold_count.append(threshold)
        
    lat = template.coord('latitude').points
    lon = template.coord('longitude').points

    lon_grid, lat_grid = np.meshgrid(lon, lat)    
        
    threshold_stack = np.stack(threshold_count, axis = 0)
    threshold_sum = np.sum(threshold_stack, axis = 0)

    cmap = mcolors.ListedColormap([
        "#D9D9D9",  # 0 members
        "#5AB4AC",  # 1
        "#99D594",  # 2
        "#E6F5B0",  # 3
        "#F6B661",  # 4
        "#F88C34",  # 5
        "#F85B34",  # 6
        "#E31414",  # 7+
    ])

    bounds = [0,1,2,3,4,5,6,7,8]
    norm = mcolors.BoundaryNorm(bounds, cmap.N)
    midpoints = [(bounds[i] + bounds[i+1]) / 2 for i in range(len(bounds)-1)]
        
    plt.figure(figsize = (12,6))
    ax = plt.axes(projection = ccrs.PlateCarree(central_longitude=0))
    ax.set_extent([-180,180, lat_minimum, lat_maximum], crs = ccrs.PlateCarree(central_longitude=0))

    for spine in ax.spines.values():
        spine.set_linewidth(3)   
        spine.set_edgecolor('black')
        
    ax.coastlines()
    ax.add_feature(cfeature.LAND, facecolor='lightgray')
    ax.add_feature(cfeature.OCEAN, facecolor='white', zorder= 3)

    mesh = ax.pcolormesh(lon_grid, lat_grid, threshold_sum, cmap = cmap, norm = norm, shading = 'auto', transform = ccrs.PlateCarree(central_longitude=0))

    cb = plt.colorbar(mesh, orientation = 'horizontal', pad = 0.05, aspect = 50)
    cb.set_label(f'Number of HRMIP members with tail width trend (> {y})', fontsize=13)
    cb.ax.tick_params(labelsize=12)
    cb.set_ticks(midpoints)
    cb.set_ticklabels(['0', '1', '2', '3', '4', '5', '6', '7+'])

    plt.title(f'Extreme tail width trends in HighResMIP members', fontsize=18, fontweight='bold')
    save_figure(plt, f'HighResMIP tail width trend threshold ({y})', 'Maps')
    plt.show()

    


    