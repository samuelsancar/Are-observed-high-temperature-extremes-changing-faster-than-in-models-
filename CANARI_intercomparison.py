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
from Lat_Land_mask_function import lat_land_mask
from Save_figure_function import save_figure

member = sys.argv[1]

# Load regridding template
template = iris.load_cube('template_cube_1deg.nc')

lat_minimum = -60
lat_maximum = 90

# Load ERA5 agreement plot
ERA5 = np.load('/home/users/samuel.sancar/Documents/Heat_Extremes/Model_agreement/CANARI_ERA5_agreement.npy')
ERA5_area = np.load('/home/users/samuel.sancar/Documents/Heat_Extremes/Model_agreement/CANARI_ERA5_area.npy')
JRA55 = np.load('/home/users/samuel.sancar/Documents/Heat_Extremes/Model_agreement/CANARI_JRA55_agreement.npy')
outdir = '/data/users/samuel.sancar/Heat_extremes/Trends/CANARI_trends/'

# Load ERA5 trends
datadir = '/data/users/samuel.sancar/Heat_extremes/Trends/Reanalysis_trends/'
era5 = iris.load_cube(f"{datadir}era5_daily_maximum_2m_temperature_tailwidth_trends_1958-2023.nc")
jradir = '/data/users/samuel.sancar/Heat_extremes/Trends/Reanalysis_trends/'
jra55 = iris.load_cube(f"{jradir}jra55__minmax_surf.015_tmax.reg_tl319._tailwidth_trends_1958-2023.nc")

# Get list of filenames
CANARI_filenames = sorted(glob.glob('/data/users/samuel.sancar/Heat_extremes/Trends/CANARI_trends/*_tailwidth_trend.nc'))

# Load cubes and pair with filenames
CANARI_files = [(iris.load_cube(fname), fname) for fname in CANARI_filenames]

CANARI_regridded = [(cube.regrid(template, iris.analysis.Linear()), fname) for cube, fname in CANARI_files]
ERA5_regrided = era5.regrid(template, iris.analysis.Linear())
#JRA55_regridded = jra55.regrid(template, iris.analysis.Linear())

# Apply mask to ERA5 data 
era5_masked_cube, _ = lat_land_mask(ERA5_regrided.data, -60, 90)
era5_masked = era5_masked_cube.data
jra55_masked_cube, _ = lat_land_mask(jra55.data, -60, 90)
jra55_masked = jra55_masked_cube.data
# Store threshold curves for all models
all_fraction_gridpoints = []
fraction_area_list = []

# Loop through each model
for i, (canari_member, fname) in enumerate(CANARI_regridded):
    indices = [j for j in range(len(CANARI_regridded)) if j != i]
    model_stack = np.concatenate(
        [np.array([CANARI_regridded[j][0].data for j in indices]),
         ERA5_regrided.data[np.newaxis, :, :]],  axis=0) # model_stack = np.concatenate([np.array([CANARI_regridded[j][0].data for j in indices]),ERA5_regrided.data[np.newaxis, :, :]], axis=0), JRA55_regridded.data[np.newaxis, :, :]

    agreement = np.mean(model_stack <= canari_member.data, axis=0)
    agreement_land_cube, n_gridpoints = lat_land_mask(agreement, -60, 90)
    agreement_land = agreement_land_cube.data
    
    # Get the base filename 
    filename = os.path.splitext(os.path.basename(fname))[0]  

    # Combine first and second parts
    model_name = filename.split('_')[0]

    # Extract lat/lon
    lats = template.coord('latitude').points
    lons = template.coord('longitude').points
    lon2d, lat2d = np.meshgrid(lons, lats)
    
    print(model_name)
    
    if model_name == member:
        
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
        cbar.set_label(f'Fraction of models with trends less than {fname}', fontsize=11)
        plt.title(f'Fraction of ensemble members where trends are less than \n {fname} trends (per gridpoint)', fontsize=12, weight='bold')
        plt.tight_layout()
        save_figure(plt, f'{fname}', 'Maps')
        plt.show()
    
    np.save(outdir + f'{model_name}_agreement.npy', agreement_land)

    # Calculate threshold curve
    thresholds = np.arange(0, 1.01, 0.01)
    fraction_gridpoints = [
        100 * np.count_nonzero((agreement_land >= thresh) & ~np.isnan(agreement_land)) / n_gridpoints
        for thresh in thresholds
    ]
    
    full = 100 * np.count_nonzero((agreement_land >= 1) & ~np.isnan(agreement_land)) / n_gridpoints
    number = np.count_nonzero((agreement_land >= 1) & ~np.isnan(agreement_land))
    print(f"Model {model_name} has {full}% of the globe ({number} grid points) covered by grid points where it's trend is greater than 100% of models.")
    
    all_fraction_gridpoints.append((fname, fraction_gridpoints))
    
    fraction_area = []
    
    for thresh in thresholds:
        
        grid_points = (agreement_land >= thresh) & (~np.isnan(agreement_land))
        cube_masked = agreement_land_cube.copy()
        cube_masked.data = grid_points.astype(int)  

        if not cube_masked.coord('latitude').has_bounds():
            cube_masked.coord('latitude').guess_bounds()
        if not cube_masked.coord('longitude').has_bounds():
            cube_masked.coord('longitude').guess_bounds()
        
        # Fix units if missing or unknown
        if cube_masked.coord('latitude').units.is_unknown():
            cube_masked.coord('latitude').units = cf_units.Unit('degrees')
        if cube_masked.coord('longitude').units.is_unknown():
            cube_masked.coord('longitude').units = cf_units.Unit('degrees')

        area_weights = iris.analysis.cartography.area_weights(cube_masked)

        area = cube_masked.collapsed(['latitude', 'longitude'], iris.analysis.SUM, weights=area_weights).data
        
        # Convert to percentage of total area
        total_area = np.sum(area_weights)
        fraction_area.append(100 * area / total_area)
        
    fraction_area_list.append((fname, fraction_area))

# Plot all threshold curves together
plt.figure(figsize=(8, 8))
y = -np.arange(0, 101, 1) + 100

for fname, fraction_gridpoints in all_fraction_gridpoints:
    plt.plot(np.arange(0, 101, 1), fraction_gridpoints, marker='.', color = 'grey') 

plt.plot(np.arange(0, 101, 1), y, linestyle='--', color='black', label='1:1 Line')
plt.plot(np.arange(0, 101, 1), ERA5, marker='.', color='red', label='ERA5 agreement')
#plt.plot(np.arange(0, 101, 1), JRA55, marker='.', color='blue', label='JRA-55 agreement')

plt.xlabel('Percentage of models where trends are less X (where X is an ensemble member, ERA5 or JRA-55) per grid point')
plt.ylabel('Percentage of Land Grid Points (%)')
plt.title('Model Agreement Curves for ERA5 and SST Models')
plt.grid(True)
plt.legend(loc='best', fontsize='small')
plt.tight_layout()
plt.show()

# Plot all threshold curves together
plt.figure(figsize=(8, 8))
y = -np.arange(0, 101, 1) + 100

for fname, fraction_area in fraction_area_list:
    plt.plot(np.arange(0, 101, 1), fraction_area, marker='.', color = 'grey') 

plt.plot(np.arange(0, 101, 1), y, linestyle='--', color='black', label='1:1 Line')
plt.plot(np.arange(0, 101, 1), ERA5_area, marker='.', color='red', label='ERA5 agreement')
#plt.plot(np.arange(0, 101, 1), JRA55, marker='.', color='blue', label='JRA-55 agreement')

plt.xlabel('Percentage of models where trends are less X (where X is an ensemble member, ERA5 or JRA-55) per grid point')
plt.ylabel('Percentage of Land Grid Points (%)')
plt.title('Model Agreement Curves for ERA5 and SST Models')
plt.grid(True)
plt.legend(loc='best', fontsize='small')
plt.tight_layout()
plt.show()



