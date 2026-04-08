"""
10/11/2025
----------
Code to generate analysis plots of models.

"""

import glob
import iris
import numpy as np
import matplotlib.pyplot as plt
import iris.plot as iplt
import matplotlib.cm as mpl_cm
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from Lat_Land_mask_function import lat_land_mask

cmap_rb = mpl_cm.get_cmap('coolwarm', 2)

t1 = 1979

# Load template and model data
template = iris.load_cube('template_cube_1deg.nc')

era5_cube = iris.load_cube(f'/data/users/samuel.sancar/Heat_extremes/Trends/Reanalysis_trends/era5_daily_maximum_2m_temperature_tailwidth_trends_{t1}-2023.nc')
jra_cube = iris.load_cube(f'/data/users/samuel.sancar/Heat_extremes/Trends/Reanalysis_trends/jra55__minmax_surf.015_tmax.reg_tl319._tailwidth_trends_{t1}-2023.nc')

era5 = era5_cube.data
jra = jra_cube.data

# Trend agreement logic
robust_trends = np.where((era5 > 0) & (jra > 0) | (era5 <= 0) & (jra <= 0), 1, np.where((era5 < 0) & (jra > 0) | (era5 > 0) & (jra < 0), -1, 0))

masked_result, _ = lat_land_mask(robust_trends, -60, 90)

fig = plt.figure(figsize=(14, 8))
ax = plt.axes(projection=ccrs.Robinson(central_longitude=0))

map = iplt.pcolormesh(masked_result, cmap=cmap_rb, vmin=-1, vmax=1)

ax.coastlines(resolution='50m', linewidth=1)
ax.add_feature(cfeature.LAND, facecolor='none')
ax.add_feature(cfeature.OCEAN, facecolor='white')

cbar = plt.colorbar(map, ax=ax, orientation='horizontal', pad=0.05, aspect=50)
cbar.set_ticks([-1,1])
cbar.ax.tick_params(labelsize=14)
cbar.set_label('Sign agreement (-1 for disagreement, 1 for agreement)', fontsize=15)
    
plt.title('ERA5 and JRA-55 Trend agreement', fontsize=14)
plt.tight_layout()
plt.show()



