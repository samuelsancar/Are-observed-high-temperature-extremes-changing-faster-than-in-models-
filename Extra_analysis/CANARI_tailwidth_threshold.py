''''
15/03/2026
----------

Plot location of CANARI tail width trends above a certain threshold, and compare with ERA5.

To run: python CANARI_tailwidth_threshold.py
'''

import iris 
import numpy as np 
import pandas as pd
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import matplotlib.cm as mpl_cm
import matplotlib.colors as mcolors
from Save_figure_function import save_figure
import cartopy.feature as cfeature

from Lat_Land_mask_function import lat_land_mask

# Beginning & end years+1
t1 = 1958
t2 = 2023

lat_minimum = -60
lat_maximum = 90

template = iris.load_cube('template_cube_1deg.nc')

thresholds = np.linspace(0.3, 0.35, 2)

era5_trend = iris.load_cube(f'/data/users/samuel.sancar/Heat_extremes/Trends/Reanalysis_trends/era5_daily_maximum_2m_temperature_tailwidth_trends_{t1}-{t2-1}.nc')
era5_trend_masked, _ = lat_land_mask(era5_trend, template, lat_minimum, lat_maximum)

canari = []

for n, i in enumerate(np.arange(1,41)):
                
    trend = iris.load_cube(f'/data/users/samuel.sancar/Heat_extremes/Trends/CANARI_trends/{i}_tailwidth_trend_{t1}-{t2-1}.nc')         
    trend_masked, _ = lat_land_mask(trend, template, lat_minimum, lat_maximum)
    
    canari.append(trend_masked)
    
for y in thresholds:
    
    canari_list = []
    
    full_list_threshold = []

    canari_full_threshold = [np.ma.filled(cube.data > y, False) for cube in canari]
    era5_threshold = np.ma.filled(era5_trend_masked.data > y, False)
        
    for k in range(len(canari)):
        left_out = canari[k]
        remaining = canari[:k] + canari[k+1:]
                
        #Threshold tail width trend:
        threshold = [np.ma.filled(cube.data > y, False) for cube in remaining]
        left_out_threshold = np.ma.filled(left_out.data > y, False)

        threshold_stack = np.stack(threshold, axis = 0)
        threshold_sum = np.sum(threshold_stack, axis = 0)    
        condition = np.ma.filled(threshold_sum == 0, False)
            
        selected_lonely_trends = left_out_threshold & condition   
        selected_sum = selected_lonely_trends.sum()
        
        canari_list.append(selected_sum)
        
        lat = template.coord('latitude').points
        lon = template.coord('longitude').points

        lon_grid, lat_grid = np.meshgrid(lon, lat)    

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
        
        # Single member
        selec_lats = lat_grid[left_out_threshold]
        selec_lons = lon_grid[left_out_threshold]  
        
        plt.figure(figsize = (12,6))
        ax = plt.axes(projection = ccrs.PlateCarree(central_longitude=0))
        ax.set_extent([-180,180, lat_minimum, lat_maximum], crs = ccrs.PlateCarree(central_longitude=0))

        for spine in ax.spines.values():
            spine.set_linewidth(3)   
            spine.set_edgecolor('black')
            
        ax.coastlines()
        ax.add_feature(cfeature.LAND, facecolor='lightgray')
        ax.add_feature(cfeature.OCEAN, facecolor='white', zorder= 3)

        # Plot stippling
        ax.scatter(selec_lons, selec_lats, color='red', s=1, alpha=1, transform=ccrs.PlateCarree(central_longitude=0), label='Stippling')

        plt.title(f'Extreme tail width trends in CANARI member {k} (>{y})', fontsize=18, fontweight='bold')
        save_figure(plt, f'CANARI member {k} tail width trend threshold ({y})', 'Maps')
        plt.show()
        plt.close()
        
        #print(sum_era5)
        #print(sum_random)
          
    full_threshold_stack = np.stack(canari_full_threshold, axis = 0)  
    full_threshold_sum = np.sum(full_threshold_stack, axis = 0)        
    full_condition = np.ma.filled(full_threshold_sum == 0, False)        
    era5_lonely_trends = era5_threshold & full_condition         
    sum_era5 = era5_lonely_trends.sum()
    
    ax = plt.gca()
    for spine in ax.spines.values():
        spine.set_linewidth(2)

    plt.hist(canari_list, bins=10, histtype='stepfilled', color='blue', linewidth=3, linestyle='-.', label='CANARI_members')
    plt.axvline(x=sum_era5, color='black', linestyle='-', linewidth=2.5, label='ERA5')
    plt.xlabel('Number of lonely tail widths', fontsize=17, fontweight='bold')
    plt.xticks(fontsize = 16)

    plt.ylabel(f'Frequency', fontsize=17, fontweight='bold')
    plt.yticks(fontsize = 16)

    plt.title(f'ERA5 vs CANARI lonely trends for (trends > {y})', fontsize=19, fontweight='bold')
    plt.legend(loc='upper right', fontsize=17)
    plt.grid(True, linestyle='--', alpha=0.6)
    plt.tight_layout()
    save_figure(plt, f'Lonely_tail_width_trends_({y})', 'Graphs')
    plt.show()
    plt.close()
    
        
        #lat = template.coord('latitude').points
        #lon = template.coord('longitude').points

        #lon_grid, lat_grid = np.meshgrid(lon, lat)    

        #cmap = mcolors.ListedColormap([
        #    "#D9D9D9",  # 0 members
        #    "#5AB4AC",  # 1
        #    "#99D594",  # 2
        #    "#E6F5B0",  # 3
        #    "#F6B661",  # 4
        #    "#F88C34",  # 5
        #    "#F85B34",  # 6
        #    "#E31414",  # 7+
        #])

        #bounds = [0,1,2,3,4,5,6,7,8]
        #norm = mcolors.BoundaryNorm(bounds, cmap.N)
        #midpoints = [(bounds[i] + bounds[i+1]) / 2 for i in range(len(bounds)-1)]
            
        #plt.figure(figsize = (12,6))
        #ax = plt.axes(projection = ccrs.PlateCarree(central_longitude=0))
        #ax.set_extent([-180,180, lat_minimum, lat_maximum], crs = ccrs.PlateCarree(central_longitude=0))

        #for spine in ax.spines.values():
        #    spine.set_linewidth(3)   
        #    spine.set_edgecolor('black')
            
        #ax.coastlines()
        #ax.add_feature(cfeature.LAND, facecolor='lightgray')
        #ax.add_feature(cfeature.OCEAN, facecolor='white', zorder= 3)

        #mesh = ax.pcolormesh(lon_grid, lat_grid, threshold_sum, cmap = cmap, norm = norm, shading = 'auto', transform = ccrs.PlateCarree(central_longitude=0))

        #cb = plt.colorbar(mesh, orientation = 'horizontal', pad = 0.05, aspect = 50)
        #cb.set_label(f'Number of members with tail width trend > {y}', fontsize=13)
        #cb.ax.tick_params(labelsize=12)
        #cb.set_ticks(midpoints)
        #cb.set_ticklabels(['0', '1', '2', '3', '4', '5', '6', '7+'])

        #plt.title(f'Extreme tail width trends in CANARI members (>{y})', fontsize=18, fontweight='bold')
        #save_figure(plt, f'CANARI tail width trend threshold ({y})', 'Maps')
        #plt.show()
        #plt.close()
        
        # Single member
        #selec_lats = lat_grid[selected]
        #selec_lons = lon_grid[selected]  
        
        #plt.figure(figsize = (12,6))
        #ax = plt.axes(projection = ccrs.PlateCarree(central_longitude=0))
        #ax.set_extent([-180,180, lat_minimum, lat_maximum], crs = ccrs.PlateCarree(central_longitude=0))

        #for spine in ax.spines.values():
        #    spine.set_linewidth(3)   
        #    spine.set_edgecolor('black')
            
        #ax.coastlines()
        #ax.add_feature(cfeature.LAND, facecolor='lightgray')
        #ax.add_feature(cfeature.OCEAN, facecolor='white', zorder= 3)

        # Plot stippling
        #ax.scatter(selec_lons, selec_lats, color='red', s=1, alpha=1, transform=ccrs.PlateCarree(central_longitude=0), label='Stippling')

        #plt.title(f'Extreme tail width trends in CANARI member {random} (>{y})', fontsize=18, fontweight='bold')
        #save_figure(plt, f'CANARI member {random} tail width trend threshold ({y})', 'Maps')
        #plt.show()
        #plt.close()
        
        #print(sum_era5)
        #print(sum_random)