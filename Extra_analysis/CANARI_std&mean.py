"""
07/11/25
---------

Code to calculate standard deviation in model trands for daily tasmax tailwidths.

To run: python Model_std.py
"""

import iris
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from Save_figure_function import save_figure

template = iris.load_cube('template_cube_1deg.nc')

CANARI_files = iris.load('/data/users/samuel.sancar/Heat_extremes/Trends/CANARI_trends/*_tailwidth_trend.nc')

CANARI_regridded = [cube.regrid(template, iris.analysis.Linear()) for cube in CANARI_files]

lsm = iris.load_cube('lsm_1279l4_0.1x0.1.grb_v4_unpack.nc')
lsm_regridded = lsm.regrid(template, iris.analysis.Nearest())

land = np.squeeze(lsm_regridded.data) > 0.5

latitudes = template.coord('latitude').points
lat_mask = latitudes >= -60
lat_mask_2d = np.repeat(lat_mask[:, np.newaxis], land.shape[1], axis=1)

land_mask = land & lat_mask_2d

def clean_data(cube):
    data = cube.data
    if hasattr(data, 'mask'):
        data = np.ma.filled(data, np.nan)
    data = np.where(np.abs(data) > 1e10, np.nan, data) 
    return np.nan_to_num(data, nan=0.0)

# Stack models into 3D array
CANARI_array = np.array([clean_data(cube) for cube in CANARI_regridded])

# Compute standard deviation and mean
std = np.nanstd(CANARI_array, axis=0)
mean = np.nanmean(CANARI_array, axis=0)

print(f"STD min: {np.nanmin(std):.4f}, max: {np.nanmax(std):.4f}")
print(f"Mean min: {np.nanmin(mean):.4f}, max: {np.nanmax(mean):.4f}")

# Land mask
std_land = np.where(land_mask, std, np.nan)
mean_land = np.where(land_mask, mean, np.nan)

def plot_map(data, title, cmap, vmin=None, vmax=None, label=''):
    fig = plt.figure(figsize=(12, 6))
    ax = plt.axes(projection=ccrs.Robinson(central_longitude=0))
    ax.set_extent([-180, 180, 90, -60], crs=ccrs.PlateCarree(central_longitude=0))
    lats = template.coord('latitude').points
    lons = template.coord('longitude').points
    lon2d, lat2d = np.meshgrid(lons, lats)

    c = ax.pcolormesh(lon2d, lat2d, data, transform=ccrs.PlateCarree(), cmap=cmap, vmin=vmin, vmax=vmax, shading='auto')
    ax.coastlines(resolution='50m', linewidth=1)
    ax.add_feature(cfeature.LAND, facecolor='none')
    ax.add_feature(cfeature.OCEAN, facecolor='white')
    plt.colorbar(c, ax=ax, orientation='horizontal', label=label, pad=0.05, aspect=50)
    plt.title(title)
    plt.tight_layout()
    save_figure(plt, title, 'Maps' )
    plt.show()

plot_map(std_land, 'CANARI tailwidth trends std', cmap='coolwarm', vmin=0, vmax=0.175, label='Std Dev')
plot_map(mean_land, 'CANARI mean tailwidth trends', cmap='RdBu_r', vmin=-0.2, vmax=0.2, label='Mean Trend')