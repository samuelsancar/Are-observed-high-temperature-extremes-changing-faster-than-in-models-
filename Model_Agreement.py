import iris
import numpy as np
import matplotlib.pyplot as plt
import scipy.special
import os
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.cm as mpl_cm
from BnD_error_function import binomial_uncertainty
from Lat_Land_mask_function import lat_land_mask

# Load template
template = iris.load_cube('template_cube_1deg.nc')

# Load model trends
coupled_files = iris.load('/data/users/samuel.sancar/Heat_extremes/Trends/Coupled_model_trends/*')
sst_files = iris.load('/data/users/samuel.sancar/Heat_extremes/Trends/SST_model_trends/*')

# Regrid model trends
coupled_regridded = [cube.regrid(template, iris.analysis.Linear()) for cube in coupled_files]
sst_regridded = [cube.regrid(template, iris.analysis.Linear()) for cube in sst_files]

# Load and regrid ERA5 trends
era5 = iris.load_cube('/data/users/samuel.sancar/Heat_extremes/Trends/Reanalysis_trends/era5_daily_maximum_2m_temperature_tailwidth_trends_1958-2023.nc')
era5_regridded = era5.regrid(template, iris.analysis.Linear())

jra = iris.load_cube('/data/users/samuel.sancar/Heat_extremes/Trends/Reanalysis_trends/jra55__minmax_surf.015_tmax.reg_tl319._tailwidth_trends_1958-2023.nc')
#jra_regridded = jra.regrid(template, iris.analysis.Linear())

# Apply mask to ERA5 data 
era5_masked_cube, _ = lat_land_mask(era5_regridded.data, -60, 90)
era5_masked = era5_masked_cube.data

jra_masked_cube, _ = lat_land_mask(jra.data, -60, 90)
jra_masked = jra_masked_cube.data

# Stack all models into a 3D array (n_models, lat, lon)
model_stack = np.concatenate([np.array([cube.data for cube in coupled_regridded]),
                              np.array([cube.data for cube in sst_regridded])], axis=0)

# for each grid point, the fraction of models where model trend <= ERA5
agreement = np.mean(model_stack <= era5_masked, axis=0)
agreement_jra = np.mean(model_stack <= jra_masked, axis=0)

# Apply new land and latitude mask
agreement_masked_cube, n_gridpoints = lat_land_mask(agreement, -60, 90)
agreement_land = agreement_masked_cube.data
agreement_jra_masked_cube, n_gridpoints_jra = lat_land_mask(agreement_jra, -60, 90)
agreement_land_jra = agreement_jra_masked_cube.data

# Compute sampling error
p_low_list = []
p_high_list = []

for i in np.arange(0, 101, 1):
    k = (32 / 100) * i if i != 0 else np.inf
    p_low, p_high, _, _, _ = binomial_uncertainty(32, k)
    p_high_list.append(100 / 32 * p_high)
    p_low_list.append(100 / 32 * p_low)

# Scan over thresholds
thresholds = np.arange(0, 1.01, 0.01)
fraction_gridpoints = []
fraction_gridpoints_jra = []

for thresh in thresholds:
    grid_points = (agreement_land >= thresh) & (~np.isnan(agreement_land))
    grid_points_jra = (agreement_land_jra >= thresh) & (~np.isnan(agreement_land_jra))
    k = np.count_nonzero(grid_points)
    k_jra = np.count_nonzero(grid_points_jra)
    # Use n_gridpoints from ERA5 mask for consistency
    fraction_gridpoints.append(100 * k / n_gridpoints)
    fraction_gridpoints_jra.append(100 * k_jra / n_gridpoints_jra)

# 1:1 line
y = -np.arange(0, 101, 1) + 100

# Compute error margins
lower_errors = np.maximum(0, np.array(p_low_list))
upper_errors = np.maximum(0, np.array(p_high_list))
xerr = [lower_errors, upper_errors]

# Save model agreement
outpath = "/home/users/samuel.sancar/Documents/Heat_Extremes/Model_agreement/"
outpath_jra = "/home/users/samuel.sancar/Documents/Heat_Extremes/Model_agreement/"
filename = "ERA5_agreement.npy"
filename_jra = "JRA55_agreement.npy"    
full_path = os.path.join(outpath, filename)
full_path_jra = os.path.join(outpath_jra, filename_jra)
os.makedirs(outpath, exist_ok=True)
np.save(full_path, fraction_gridpoints)
np.save(full_path_jra, fraction_gridpoints_jra)

# Plotting threshold plot
plt.figure(figsize=(12, 7))
plt.errorbar(thresholds * 100, fraction_gridpoints, xerr=xerr, fmt='.',
             color='darkblue', ecolor='salmon', elinewidth=2, capsize=5,
             label='Model Agreement')
plt.errorbar(thresholds * 100, fraction_gridpoints_jra, xerr=xerr, fmt='.',
             color='darkgreen', ecolor='lightgreen', elinewidth=2, capsize=5,
             label='Model Agreement with JRA-55')
plt.plot(thresholds * 100, y, linestyle='--', color='gray', linewidth=2, label='1:1 Line')
plt.xlabel('Percentage of Models Less Than ERA5 and JRA-55 at Grid Point (%)', fontsize=14, weight='bold')
plt.ylabel('Percentage of Land Grid Points (%)', fontsize=14, weight='bold')
plt.title('Model Agreement with ERA5 Across Land Grid Points', fontsize=16, weight='bold')
plt.grid(True, linestyle='--', alpha=0.6)
plt.legend(fontsize=12)
plt.tight_layout()
plt.show()

# Extract latitude and longitude points
lats = template.coord('latitude').points
lons = template.coord('longitude').points
lon2d, lat2d = np.meshgrid(lons, lats)

# Create figure and axis
fig = plt.figure(figsize=(14, 7))
ax = plt.axes(projection=ccrs.Robinson(central_longitude=0))
ax.set_extent([-180, 180, -60, 90], crs=ccrs.PlateCarree(central_longitude=0))

# Use colormap
cmap = mpl_cm.get_cmap('coolwarm', 20)

# Plot the data
c = ax.pcolormesh(lon2d, lat2d, agreement_land, cmap=cmap, vmin=0, vmax=1, shading='auto', transform=ccrs.PlateCarree(central_longitude=0))

# Add refined map features
ax.coastlines(resolution='50m', linewidth=0.8)
ax.add_feature(cfeature.BORDERS.with_scale('50m'), linewidth=0.5)
ax.add_feature(cfeature.LAND, facecolor='lightgray')
ax.add_feature(cfeature.OCEAN, facecolor='aliceblue')

# Add gridlines with labels
gl = ax.gridlines(draw_labels=True, linewidth=0.5, color='gray', alpha=0.5, linestyle='--')
gl.top_labels = False
gl.right_labels = False

# Add colourbar
cbar = plt.colorbar(c, ax=ax, orientation='horizontal', pad=0.05, aspect=50)
cbar.set_label('Fraction of Models Less Than ERA5 trends', fontsize=12)

# Add title
plt.title('Fraction of ensemble members where trends are less than ERA5 trends (per grid point)', fontsize=16, weight='bold')

# Final layout adjustments
plt.tight_layout()
plt.savefig("/home/users/samuel.sancar/Documents/Heat_Extremes/CSSP_Images/Model_agreement.png", dpi = 500)
plt.show()

plt.figure(figsize=(14, 7))
ax = plt.axes(projection=ccrs.PlateCarree())
ax.set_extent([-180, 180, -60, 90], crs=ccrs.PlateCarree())
cmap = mpl_cm.get_cmap('coolwarm', 20)
c = ax.pcolormesh(lon2d, lat2d, agreement_land_jra, cmap=cmap, vmin=0, vmax=1, shading='auto')          
ax.coastlines(resolution='50m', linewidth=0.8)
ax.add_feature(cfeature.BORDERS.with_scale('50m'), linewidth=0.5)
ax.add_feature(cfeature.LAND, facecolor='lightgray')
ax.add_feature(cfeature.OCEAN, facecolor='aliceblue')
gl = ax.gridlines(draw_labels=True, linewidth=0.5, color='gray', alpha=0.5, linestyle='--')
gl.top_labels = False
gl.right_labels = False
cbar = plt.colorbar(c, ax=ax, orientation='vertical', pad=0.02, aspect=40)
cbar.set_label('Fraction of Models Less Than JRA-55', fontsize=12)
plt.title('Model Agreement: Fraction of Models Less Than JRA-55 (per grid point)', fontsize=16, weight='bold')
plt.tight_layout()
plt.show()



