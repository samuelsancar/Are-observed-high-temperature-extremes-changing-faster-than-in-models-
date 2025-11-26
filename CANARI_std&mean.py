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

# Load template and data cubes
template = iris.load_cube('template_cube_1deg.nc')

CANARI_files = iris.load('/data/users/samuel.sancar/Heat_extremes/Trends/CANARI_trends/*_tailwidth_trend.nc')

# Regrid all cubes to template
CANARI_regridded = [cube.regrid(template, iris.analysis.Linear()) for cube in CANARI_files]

# Load and regrid land-sea mask
lsm = iris.load_cube('lsm_1279l4_0.1x0.1.grb_v4_unpack.nc')
lsm_regridded = lsm.regrid(template, iris.analysis.Nearest())

land = np.squeeze(lsm_regridded.data) > 0.5

# Latitude mask (exclude Antarctica)
latitudes = template.coord('latitude').points
lat_mask = latitudes >= -60
lat_mask_2d = np.repeat(lat_mask[:, np.newaxis], land.shape[1], axis=1)

land_mask = land & lat_mask_2d

# Clean data function
def clean_data(cube):
    data = cube.data
    if hasattr(data, 'mask'):
        data = np.ma.filled(data, np.nan)
    data = np.where(np.abs(data) > 1e10, np.nan, data)  # remove absurd fill values
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

# Plotting functions
def plot_map(data, title, cmap, vmin=None, vmax=None, label=''):
    fig = plt.figure(figsize=(12, 6))
    ax = plt.axes(projection=ccrs.Robinson(central_longitude=0))
    lats = template.coord('latitude').points
    lons = template.coord('longitude').points
    lon2d, lat2d = np.meshgrid(lons, lats)

    c = ax.pcolormesh(lon2d, lat2d, data, transform=ccrs.PlateCarree(), cmap=cmap, vmin=vmin, vmax=vmax, shading='auto')
    ax.coastlines(resolution='50m', linewidth=1)
    ax.add_feature(cfeature.LAND, facecolor='none')
    ax.add_feature(cfeature.OCEAN, facecolor='white')
    plt.colorbar(c, ax=ax, orientation='horizontal', label=label)
    plt.title(title)
    plt.tight_layout()
    plt.show()

# Plot STD and Mean
plot_map(std_land, 'Standard deviation of model trends (per grid point)', cmap='coolwarm', vmin=0, vmax=0.175, label='Std Dev')
plot_map(mean_land, 'Mean of model trends (per grid point)', cmap='RdBu_r', vmin=-0.2, vmax=0.2, label='Mean Trend')