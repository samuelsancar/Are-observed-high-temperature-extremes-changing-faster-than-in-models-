"""
    In progress
"""

import iris
import xarray as xr
import os
import numpy as np
import matplotlib.pyplot as plt
import glob
import cartopy.crs as ccrs
import cartopy.feature as cfeature

template = iris.load_cube('template_cube_1deg.nc')

coupled_files = glob.glob('/data/users/samuel.sancar/Heat_extremes/trends/Coupled_model_trends/*')
sst_files = glob.glob('/data/users/samuel.sancar/Heat_extremes/trends/SST_model_trends/*')


coupled_regridded = [cube.regrid(template, iris.analysis.Linear()) for cube in coupled_files]
sst_regridded = [cube.regrid(template, iris.analysis.Linear()) for cube in sst_files]

trend_sum = 


for files in coupled_regridded:

    main = files[0].split('/')[-1].split('_')[0]


import iris
import numpy as np
import matplotlib.pyplot as plt

template = iris.load_cube('template_cube_1deg.nc')

coupled_files = iris.load('/data/users/samuel.sancar/Heat_extremes/trends/Coupled_model_trends/*')
sst_files = iris.load('/data/users/samuel.sancar/Heat_extremes/trends/SST_model_trends/*')

coupled_regridded = [cube.regrid(template, iris.analysis.Linear()) for cube in coupled_files]
sst_regridded = [cube.regrid(template, iris.analysis.Linear()) for cube in sst_files]

era5 = iris.load_cube('/data/users/samuel.sancar/Heat_extremes/trends/Reanalysis_trends/era5_daily_maximum_2m_temperature_tailwidth_trends_1958-2023.nc')
era5_regridded = era5.regrid(template, iris.analysis.Linear())

# Load and regrid land-sea mask
lsm = iris.load_cube('lsm_1279l4_0.1x0.1.grb_v4_unpack.nc')
lsm_regridded = lsm.regrid(template, iris.analysis.Nearest())

# Mask for land
land = np.squeeze(lsm_regridded.data) > 0.5
no_of_gridpoints = np.count_nonzero(land)

# Get latitude points
latitudes = template.coord('latitude').points

# Create a 2D mask for latitudes >= -60 (exclude Antarctica)
lat_mask = latitudes >= -60
lat_mask_2d = np.repeat(lat_mask[:, np.newaxis], land.shape[1], axis=1)

land_mask = land & lat_mask_2d

# Stack all models into a 3D array (n_models, lat, lon)
model_stack = np.concatenate([np.array([cube.data for cube in coupled_regridded]),np.array([cube.data for cube in sst_regridded])],axis=0)  

# for each grid point, the fraction of models where model trend <= ERA5
agreement = np.mean(model_stack <= era5_regridded.data, axis=0)  

# Only consider land points
agreement_land = np.where(land_mask, agreement, np.nan)

# Scan over thresholds (0 to 1, step 0.01)
thresholds = np.arange(0, 1.01, 0.01)
n_gridpoints = np.count_nonzero(~np.isnan(agreement_land))

fraction_gridpoints = []
for thresh in thresholds:
    n = np.count_nonzero((agreement_land >= thresh) & ~np.isnan(agreement_land))
    fraction_gridpoints.append(100 * n / n_gridpoints)

# 1:1 line
y = - np.arange(0, 101, 1) + 100


# Plotting
plt.figure(figsize=(10, 6))
plt.plot(thresholds * 100, fraction_gridpoints, marker='.')
plt.plot(thresholds * 100, y, linestyle='--', color='gray', label='1:1 Line')
plt.xlabel('Minimum % of Models Less Than ERA5 at Grid Point')
plt.ylabel('Percentage of Land Grid Points (%)')
plt.title('Grid points with at least X% of models less than ERA5')
plt.grid(True)
plt.tight_layout()
plt.show()

# Plot model agreement as a map
import cartopy.crs as ccrs
import cartopy.feature as cfeature

fig = plt.figure(figsize=(12, 6))
ax = plt.axes(projection=ccrs.PlateCarree())

# lats/lons from template cube
lats = template.coord('latitude').points
lons = template.coord('longitude').points
lon2d, lat2d = np.meshgrid(lons, lats)
agreement_plot = np.where(lat2d >= -60, agreement_land, np.nan)

c = ax.pcolormesh(lon2d, lat2d, np.squeeze(agreement_plot), cmap='viridis', vmin=0, vmax=1, shading='auto')

ax.coastlines(resolution='50m', linewidth=1)
ax.add_feature(cfeature.LAND, facecolor='none')
ax.add_feature(cfeature.OCEAN, facecolor='white')
plt.colorbar(c, ax=ax, orientation='vertical', label='Fraction of Models Less Than ERA5')
plt.title('Model Agreement: Fraction of Models Less Than ERA5 (per grid point)')
plt.tight_layout()
plt.show()
    
