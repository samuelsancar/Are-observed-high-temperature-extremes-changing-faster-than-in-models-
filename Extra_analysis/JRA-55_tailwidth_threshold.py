''''
14/03/2026
----------

Plot location of JRA-55 tail width trends above a certain threshold.

To run: python JRA-55_tailwidth_threshold.py
'''

import iris 
import numpy as np 
import pandas as pd
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import iris.plot as iplt
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
    
trend = iris.load_cube(f'/data/users/samuel.sancar/Heat_extremes/Trends/Reanalysis_trends/jra55__minmax_surf.015_tmax.reg_tl319._tailwidth_trends_{t1}-{t2-1}.nc')
    
trend_masked, _ = lat_land_mask(trend, template, lat_minimum, lat_maximum)
    
threshold = np.linspace(0.1, 0.35, 6)

for i in threshold:
    #Threshold tail width trend:
    threshold = np.ma.filled(trend_masked.data > i, False)
    #indices = np.where(trend_masked)
        
    lat = template.coord('latitude').points
    lon = template.coord('longitude').points

    lon_grid, lat_grid = np.meshgrid(lon, lat)     
    
    thresh_lats = lat_grid[threshold]
    thresh_lons = lon_grid[threshold]
        
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
    ax.scatter(thresh_lons, thresh_lats, color='red', s=1, alpha=1, transform=ccrs.PlateCarree(central_longitude=0), label='Stippling')

    plt.title(f'Extreme tail width trends in JRA-55 (>{i})', fontsize=18, fontweight='bold')
    save_figure(plt, f'JRA-55 tail width trend threshold ({i})', 'Maps')
    #plt.show()