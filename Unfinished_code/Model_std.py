"""
01/08/25
---------
Code to calculate standard deviation in model trands for daily tasmax tailwidths.

To run: python Model_std.py

"""

import xarray as xr
import numpy as np
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import iris
import cartopy.feature as cfeature


template = iris.load_cube('template_cube_1deg.nc')

coupled_files = iris.load('/data/users/samuel.sancar/Heat_extremes/trends/Coupled_model_trends/*')
sst_files = iris.load('/data/users/samuel.sancar/Heat_extremes/trends/SST_model_trends/*')


coupled_regridded = [cube.regrid(template, iris.analysis.Linear()) for cube in coupled_files]
sst_regridded = [cube.regrid(template, iris.analysis.Linear()) for cube in sst_files]

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
std = np.std(model_stack, axis=0)  

# Only consider land points
std_land = np.where(land_mask, std, np.nan)

fig = plt.figure(figsize=(12, 6))
ax = plt.axes(projection=ccrs.PlateCarree())

# lats/lons from template cube
lats = template.coord('latitude').points
lons = template.coord('longitude').points
lon2d, lat2d = np.meshgrid(lons, lats)
std_plot = np.where(lat2d >= -60, std_land, np.nan)

c = ax.pcolormesh(lon2d, lat2d, np.squeeze(std_plot), cmap='coolwarm', vmin=0, vmax=0.175, shading='auto')

ax.coastlines(resolution='50m', linewidth=1)
ax.add_feature(cfeature.LAND, facecolor='none')
ax.add_feature(cfeature.OCEAN, facecolor='white')
plt.colorbar(c, ax=ax, orientation='vertical', label='Standard deviation of model trends')
plt.title('Standard deviation of model trends (per grid point)')
plt.tight_layout()
plt.show()