"""

 Code to calculate and plot model agreement intercomparison.
 
 To run: python Model_intercomparison.
 
"""

import iris
import os
import numpy as np
import matplotlib.pyplot as plt
import glob
import matplotlib.cm as mpl_cm
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from Lat_Land_mask_function import lat_land_mask
from Save_figure_function import save_figure

# Load regridding template
template = iris.load_cube('template_cube_1deg.nc')

lat_minimum = -60
lat_maximum = 90

# Load ERA5 agreement plot
ERA5 = np.load('/home/users/samuel.sancar/Documents/Heat_Extremes/Model_agreement/ERA5_agreement.npy')
JRA55 = np.load('/home/users/samuel.sancar/Documents/Heat_Extremes/Model_agreement/JRA55_agreement.npy')

# Load ERA5 trends
datadir = '/data/users/samuel.sancar/Heat_extremes/Trends/Reanalysis_trends/'
era5 = iris.load_cube(f"{datadir}era5_daily_maximum_2m_temperature_tailwidth_trends_1958-2023.nc")
jradir = '/data/users/samuel.sancar/Heat_extremes/Trends/Reanalysis_trends/'
jra55 = iris.load_cube(f"{jradir}jra55__minmax_surf.015_tmax.reg_tl319._tailwidth_trends_1958-2023.nc")

# Get list of filenames
coupled_filenames = sorted(glob.glob('/data/users/samuel.sancar/Heat_extremes/Trends/Coupled_model_trends/*_tailwidth_trend.nc'))
sst_filenames = sorted(glob.glob('/data/users/samuel.sancar/Heat_extremes/Trends/SST_model_trends/*_tailwidth_trend.nc'))

# Load cubes and pair with filenames
coupled_files = [(iris.load_cube(fname), fname) for fname in coupled_filenames]
sst_files = [(iris.load_cube(fname), fname) for fname in sst_filenames]

# Regrid data
coupled_regridded = [(cube.regrid(template, iris.analysis.Linear()), fname) for cube, fname in coupled_files]
sst_regridded = [(cube.regrid(template, iris.analysis.Linear()), fname) for cube, fname in sst_files]
ERA5_regrided = era5.regrid(template, iris.analysis.Linear())

# Apply mask to ERA5 data (if needed)
era5_masked_cube, _ = lat_land_mask(ERA5_regrided.data, -60, 90)
era5_masked = era5_masked_cube.data
jra55_masked_cube, _ = lat_land_mask(jra55.data, -60, 90)
jra55_masked = jra55_masked_cube.data
# Store threshold curves for all models
all_fraction_gridpoints = []

# Loop through each model
for i, (coupled, fname) in enumerate(coupled_regridded):
    indices = [j for j in range(len(coupled_regridded)) if j != i]
    model_stack = np.concatenate(
        [np.array([coupled_regridded[j][0].data for j in indices]),
         np.array([cube[0].data for cube in sst_regridded]),
         ERA5_regrided.data[np.newaxis, :, :]], axis=0)

    agreement = np.mean(model_stack <= coupled.data, axis=0)
    agreement_land_cube, n_gridpoints = lat_land_mask(agreement, -60, 90)
    agreement_land = agreement_land_cube.data

    # Get the base filename without path or extension
    filename = os.path.splitext(os.path.basename(fname))[0]  # "ACCESS-CM2_sst"

    # Split by underscore
    parts = filename.split('_')

    # Combine first and second parts
    model_name = f"{parts[0]}_{parts[2]}" # "ACCESS-CM2_sst"

    # Extract lat/lon
    lats = template.coord('latitude').points
    lons = template.coord('longitude').points
    lon2d, lat2d = np.meshgrid(lons, lats)

    # Plot individual map
    #fig = plt.figure(figsize=(12, 6))
    #ax = plt.axes(projection=ccrs.Robinson(central_longitude=0))
    #ax.set_extent([-180, 180, lat_minimum, lat_maximum], crs=ccrs.PlateCarree(central_longitude=0))
    
    #cmap = mpl_cm.get_cmap('coolwarm', 20)
    #c = ax.pcolormesh(lon2d, lat2d, agreement_land, cmap=cmap, vmin=0, vmax=1, shading='auto', transform=ccrs.PlateCarree(central_longitude=0))

    #ax.coastlines(resolution='50m', linewidth=0.8)
    #ax.add_feature(cfeature.BORDERS.with_scale('50m'), linewidth=0.5)
    #ax.add_feature(cfeature.LAND, facecolor='lightgray')
    #ax.add_feature(cfeature.OCEAN, facecolor='aliceblue')

    #gl = ax.gridlines(draw_labels=True, linewidth=0.5, color='gray', alpha=0.5, linestyle='--')
    #gl.top_labels = False
    #gl.right_labels = False

    #cbar = plt.colorbar(c, ax=ax, orientation='horizontal', pad=0.05, aspect=50)
    #cbar.set_label(f'Fraction of models less then {model_name}', fontsize=10)

    #plt.title(f'Model intercomparison plot for coupled {model_name} (1958-2022)', fontsize=16, weight='bold')
    #plt.tight_layout()
    #plt.show()

    # Calculate threshold curve
    thresholds = np.arange(0, 1.01, 0.01)
    fraction_gridpoints = [
        100 * np.count_nonzero((agreement_land >= thresh) & ~np.isnan(agreement_land)) / n_gridpoints
        for thresh in thresholds
    ]
    all_fraction_gridpoints.append((fname, fraction_gridpoints))

# Plot all threshold curves together
plt.figure(figsize=(10, 6))
y = -np.arange(0, 101, 1) + 100
for fname, fraction_gridpoints in all_fraction_gridpoints:
    plt.plot(np.arange(0, 101, 1), fraction_gridpoints, marker='.', color = 'grey') # Add label if needed "label=fname"

plt.plot(np.arange(0, 101, 1), y, linestyle='--', color='black', label='1:1 Line')
plt.plot(np.arange(0, 101, 1), ERA5, marker='.', color='red', label='ERA5 agreement')
plt.plot(np.arange(0, 101, 1), JRA55, marker='.', color='blue', label='JRA-55 agreement')

plt.xlabel('X % of Models Less Than X at Grid Point')
plt.ylabel('Percentage of Land Grid Points (%)')
plt.title('Model Agreement Curves for ERA5 and Coupled Models')
plt.grid(True)
plt.legend(loc='best', fontsize='small')
plt.tight_layout()
#plt.show()

# Store threshold curves for all models
all_fraction_gridpoints_sst = []

plt.figure(figsize=(12, 8))

# Looping over SST model agreement
for j, (sst, fname) in enumerate(sst_regridded):
    indices = [i for i in range(len(sst_regridded)) if i != j]
    model_stack = np.concatenate(
        [np.array([cube[0].data for cube in coupled_regridded]),
         np.array([sst_regridded[i][0].data for i in indices]),
         ERA5_regrided.data[np.newaxis, :, :]], axis=0)


    # Apply latitude and land mask to all models and ERA5
    agreement = np.mean(model_stack <= sst.data, axis=0)
    agreement_land_cube, n_gridpoints = lat_land_mask(agreement, -60, 90)
    agreement_land = agreement_land_cube.data

    # Get the base filename without path or extension
    filename = os.path.splitext(os.path.basename(fname))[0]  # "ACCESS-CM2_sst"

    # Split by underscore
    parts = filename.split('_')

    # Combine first and second parts
    model_name = f"{parts[0]}_{parts[2]}" # "ACCESS-CM2_sst"

    # Extract lat/lon
    lats = template.coord('latitude').points
    lons = template.coord('longitude').points
    lon2d, lat2d = np.meshgrid(lons, lats)

    if model_name == 'HadGEM3-GC31-HM_r1i2p1f1':
        # Plot individual map
        fig = plt.figure(figsize=(12, 6))
        ax = plt.axes(projection=ccrs.Robinson(central_longitude=0))
        ax.set_extent([-180, 180, lat_minimum, lat_maximum], crs=ccrs.PlateCarree(central_longitude=0))

        cmap = mpl_cm.get_cmap('coolwarm', 20)
        c = ax.pcolormesh(lon2d, lat2d, agreement_land, cmap=cmap, vmin=0, vmax=1, shading='auto', transform=ccrs.PlateCarree(central_longitude=0))

        ax.coastlines(resolution='50m', linewidth=0.8)
        ax.add_feature(cfeature.BORDERS.with_scale('50m'), linewidth=0.5)
        ax.add_feature(cfeature.LAND, facecolor='lightgray')
        ax.add_feature(cfeature.OCEAN, facecolor='aliceblue')

        gl = ax.gridlines(draw_labels=True, linewidth=0.5, color='gray', alpha=0.5, linestyle='--')
        gl.top_labels = False
        gl.right_labels = False

        cbar = plt.colorbar(c, ax=ax, orientation='horizontal', pad=0.05, aspect=50)
        cbar.set_label(f'Fraction of models with trends less than a) ERA5, b) HadGEM3-GC31-HM highresSST-present r1i2p1f1', fontsize=11)
        plt.title(f'Fraction of ensemble members where trends are less than \n HadGEM3-GC31-HM highresSST-present r1i2p1f1 trends (per gridpoint)', fontsize=12, weight='bold')
        plt.tight_layout()
        save_figure(plt, 'HadGEM-GC31-HM_SST_trend_intercomparison', 'Maps')
        plt.show()
        
        


    # Calculate threshold curve
    thresholds = np.arange(0, 1.01, 0.01)
    fraction_gridpoints_sst = [
        100 * np.count_nonzero((agreement_land >= thresh) & ~np.isnan(agreement_land)) / n_gridpoints
        for thresh in thresholds
    ]
    all_fraction_gridpoints_sst.append((fname, fraction_gridpoints_sst))

# Plot all threshold curves together
plt.figure(figsize=(8, 8))
y = -np.arange(0, 101, 1) + 100

for fname, fraction_gridpoints_sst in all_fraction_gridpoints_sst:
    plt.plot(np.arange(0, 101, 1), fraction_gridpoints_sst, marker='.', color = 'grey') # Add label if needed "label=fname"

plt.plot(np.arange(0, 101, 1), y, linestyle='--', color='black', label='1:1 Line')
plt.plot(np.arange(0, 101, 1), ERA5, marker='.', color='red', label='ERA5 agreement')
plt.plot(np.arange(0, 101, 1), JRA55, marker='.', color='blue', label='JRA-55 agreement')

plt.xlabel('Percentage of models where trends are less X (where X is an ensemble member, ERA5 or JRA-55) per grid point')
plt.ylabel('Percentage of Land Grid Points (%)')
plt.title('Model Agreement Curves for ERA5 and SST Models')
plt.grid(True)
plt.legend(loc='best', fontsize='small')
plt.tight_layout()
plt.show()



