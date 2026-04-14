"""
13/11/2025
----------
CANARI intercomparison
"""

import iris
import iris.analysis
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import glob
import cf_units
import matplotlib.cm as mpl_cm
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from matplotlib.colors import ListedColormap, BoundaryNorm
import matplotlib.colors as mcolors
from matplotlib import cm
from Lat_Land_mask_function import lat_land_mask
from Save_figure_function import save_figure

member = sys.argv[1]

# Load regridding template
template = iris.load_cube('template_cube_1deg.nc')

lat_minimum = -60
lat_maximum = 90

# Colour bar parameters
cmap = mcolors.ListedColormap([
    "#5AB4AC",  
    "#99D594",  
    "#E6F5B0",  
    "#F6B661",  
    "#F87234",  
])

cmap.set_under("#544B94")  # exactly 0
cmap.set_over("#9E0142")   # exactly 100
bounds = [0.1 ,5 ,25, 75, 95, 100]

# Load ERA5 agreement plot
ERA5 = np.load('/home/users/samuel.sancar/Documents/Temperature_extremes/Model_agreement/CANARI_ERA5_agreement.npy')
ERA5_area = np.load('/home/users/samuel.sancar/Documents/Temperature_extremes/Model_agreement/CANARI_ERA5_area.npy')
JRA55 = np.load('/home/users/samuel.sancar/Documents/Temperature_extremes/Model_agreement/CANARI_JRA55_agreement.npy')
outdir = '/data/users/samuel.sancar/Heat_extremes/Trends/CANARI_trends/'

# Load ERA5 trends
datadir = '/data/users/samuel.sancar/Heat_extremes/Trends/Reanalysis_trends/'
era5 = iris.load_cube(f"{datadir}era5_daily_maximum_2m_temperature_tailwidth_trends_1958-2023.nc")
ERA5_OHpcnt = np.load('/home/users/samuel.sancar/Documents/Temperature_extremes/Model_agreement/ERA5_OHpcnt_grids_CANARI.npy')
jradir = '/data/users/samuel.sancar/Heat_extremes/Trends/Reanalysis_trends/'
jra55 = iris.load_cube(f"{jradir}jra55__minmax_surf.015_tmax.reg_tl319._tailwidth_trends_1958-2023.nc")

# Get list of filenames
CANARI_filenames = sorted(glob.glob('/data/users/samuel.sancar/Heat_extremes/Trends/CANARI_trends/*_tailwidth_trend_1958-2023.nc'))

# Load cubes and pair with filenames
CANARI_files = [(iris.load_cube(fname), fname) for fname in CANARI_filenames]

# Apply mask to ERA5 data 
#era5_masked_cube, _ = lat_land_mask(era5.data, -60, 90)
era5_masked = era5.data
#jra55_masked_cube, _ = lat_land_mask(jra55.data, -60, 90)
jra55_masked = jra55.data

# Store threshold curves for all models
all_fraction_gridpoints = []
fraction_area_list = []
OHpcnt_list = []
labels_list = []

lon = template.coord('longitude').points
lat = template.coord('latitude').points
lon2d, lat2d = np.meshgrid(lon, lat)

if member == 'era5':
    p87p5_trend = iris.load_cube(f"{datadir}era5_daily_maximum_2m_temperature_p87p5_trends_1958-2023.nc")
    indices = [j for j in range(len(CANARI_files))]
    model_stack = np.concatenate(
        [np.array([CANARI_files[j][0].data for j in indices])],  axis=0) # model_stack = np.concatenate([np.array([CANARI_files[j][0].data for j in indices]),era5.data[np.newaxis, :, :]], axis=0), JRA55_regridded.data[np.newaxis, :, :]

    p87p5_masked, _ = lat_land_mask(p87p5_trend, lat_minimum, lat_maximum)
    neg_indices = np.where(p87p5_masked.data < 0)
    neg_lats = lat2d[neg_indices]
    neg_lons = lon2d[neg_indices]
    
    agreement = np.mean(model_stack <= era5.data, axis=0)
    agreement_land_cube, n_gridpoints = lat_land_mask(agreement, -60, 90)
    agreement_land = agreement_land_cube.data

    # Extract lat/lon
    lats = template.coord('latitude').points
    lons = template.coord('longitude').points
    lon2d, lat2d = np.meshgrid(lons, lats)
            
    # Plot individual map
    fig = plt.figure(figsize=(12, 6))
    ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=0))
    ax.set_extent([-180, 180, lat_minimum, lat_maximum], crs=ccrs.PlateCarree(central_longitude=0))

    norm = mcolors.BoundaryNorm(bounds, cmap.N, clip=False)
    
    for spine in ax.spines.values():
        spine.set_linewidth(2)  
        spine.set_edgecolor('black')

    c = ax.pcolormesh(lon2d, lat2d, agreement_land*100, cmap=cmap, norm=norm,
                    shading='auto', transform=ccrs.PlateCarree(central_longitude=0))
    
    gl = ax.gridlines(draw_labels=False, linewidth=0.5, color='gray', alpha=0.5, linestyle='--')
    gl.top_labels = False
    gl.right_labels = False
    
    # Overlay white dots where p87p5_trend is negative
    ax.scatter(neg_lons, neg_lats, color='white', s=0.9, alpha=0.9, transform=ccrs.PlateCarree(), label='p87p5 trend < 0')

    ax.coastlines(resolution='50m', linewidth=0.8)
    ax.add_feature(cfeature.LAND, facecolor='lightgray')
    ax.add_feature(cfeature.OCEAN, facecolor='aliceblue')

    #cbar = plt.colorbar(
        #c,
        #ax=ax,
        #orientation='horizontal',
        #pad=0.05,
        #aspect=50,
        #extend='both')

    #cbar.set_ticklabels(['0', '5', '25', '75', '95', '100'])

    #cbar.set_label(f'Percent of CANARI-LE ensemble member trends less than ERA5 [%]', fontsize=11)
    plt.title(f'ERA5', fontsize=22, weight='bold')
    plt.tight_layout()
    save_figure(plt, 'ERA5_percentile_map', 'Maps')
    plt.show()

else:
    # Loop through each model
    for i, (canari_member, fname) in enumerate(CANARI_files):
        
            # Get the base filename 
        filename = os.path.splitext(os.path.basename(fname))[0]  

        # Combine first and second parts
        model_name = filename.split('_')[0]
        
        canari_p87p5 = iris.load_cube(f'/data/users/samuel.sancar/Heat_extremes/Trends/CANARI_trends/{model_name}_p87p5_trend_1958-2023.nc')
        
        indices = [j for j in range(len(CANARI_files)) if j != i]
        model_stack = np.concatenate(
            [np.array([CANARI_files[j][0].data for j in indices]),
            era5.data[np.newaxis, :, :]],  axis=0) # model_stack = np.concatenate([np.array([CANARI_files[j][0].data for j in indices]),era5.data[np.newaxis, :, :]], axis=0), JRA55_regridded.data[np.newaxis, :, :]

        p87p5_masked, _ = lat_land_mask(canari_p87p5, lat_minimum, lat_maximum)
        neg_indices = np.where(p87p5_masked.data < 0)
        neg_lats = lat2d[neg_indices]
        neg_lons = lon2d[neg_indices]
        
        agreement = np.mean(model_stack <= canari_member.data, axis=0)
        agreement_land_cube, n_gridpoints = lat_land_mask(agreement, -60, 90)
        agreement_land = agreement_land_cube.data

        # Extract lat/lon
        lats = template.coord('latitude').points
        lons = template.coord('longitude').points
        lon2d, lat2d = np.meshgrid(lons, lats)
        
        print(model_name)
        
        if model_name == member:
            
            # Plot individual map
            fig = plt.figure(figsize=(12, 6))
            ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=0))
            ax.set_extent([-180, 180, lat_minimum, lat_maximum], crs=ccrs.PlateCarree(central_longitude=0))

            norm = mcolors.BoundaryNorm(bounds, cmap.N, clip=False)
            
            for spine in ax.spines.values():
                spine.set_linewidth(2)  
                spine.set_edgecolor('black')
                
            gl = ax.gridlines(draw_labels=False, linewidth=0.5, color='gray', alpha=0.5, linestyle='--')
            gl.top_labels = False
            gl.right_labels = False

            c = ax.pcolormesh(lon2d, lat2d, agreement_land*100, cmap=cmap, norm=norm,
                        shading='auto', transform=ccrs.PlateCarree(central_longitude=0))
            
            # Overlay white dots where p87p5_trend is negative
            ax.scatter(neg_lons, neg_lats, color='white', s=0.9, alpha=0.6, transform=ccrs.PlateCarree(), label='p87p5 trend < 0')

            ax.coastlines(resolution='50m', linewidth=0.8)
            ax.add_feature(cfeature.LAND, facecolor='lightgray')
            ax.add_feature(cfeature.OCEAN, facecolor='aliceblue')

                
            #cbar = plt.colorbar(
            #        c,
            #        ax=ax,
            #        orientation='horizontal',
            #        pad=0.05,
            #        aspect=50,
            #        extend='both'
            #    )

            #cbar.set_ticklabels(['0', '5', '25', '75', '95', '100'])
            #cbar.ax.tick_params(labelsize=16)
            #cbar.set_label(f'Percent of CANARI-LE trends less than member [%]', fontsize=18)
            plt.title(f'CANARI-LE member {model_name}', fontsize=22, weight='bold')
            plt.tight_layout()
            save_figure(plt, f'{model_name}_percentile_map', 'Maps')
            plt.show()

        

