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
coupled_files = iris.load('/data/users/samuel.sancar/Heat_extremes/trends/Coupled_model_trends/*')
sst_files = iris.load('/data/users/samuel.sancar/Heat_extremes/trends/SST_model_trends/*')

# Regrid model trends
coupled_regridded = [cube.regrid(template, iris.analysis.Linear()) for cube in coupled_files]
sst_regridded = [cube.regrid(template, iris.analysis.Linear()) for cube in sst_files]

# Load and regrid ERA5 trends
era5 = iris.load_cube('/data/users/samuel.sancar/Heat_extremes/trends/Reanalysis_trends/era5_daily_maximum_2m_temperature_tailwidth_trends_1958-2023.nc')
era5_regridded = era5.regrid(template, iris.analysis.Linear())

# Apply mask to ERA5 data 
era5_masked_cube, _ = lat_land_mask(era5_regridded.data, -60)
era5_masked = era5_masked_cube.data

# Stack all models into a 3D array (n_models, lat, lon)
model_stack = np.concatenate([np.array([cube.data for cube in coupled_regridded]),
                              np.array([cube.data for cube in sst_regridded])], axis=0)

# for each grid point, the fraction of models where model trend <= ERA5
agreement = np.mean(model_stack <= era5_masked, axis=0)

# Apply new land and latitude mask
agreement_masked_cube, n_gridpoints = lat_land_mask(agreement, -60)
agreement_land = agreement_masked_cube.data

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

for thresh in thresholds:
    grid_points = (agreement_land >= thresh) & (~np.isnan(agreement_land))
    k = np.count_nonzero(grid_points)
    fraction_gridpoints.append(100 * k / n_gridpoints)

# 1:1 line
y = -np.arange(0, 101, 1) + 100

# Compute error margins
lower_errors = np.maximum(0, np.array(p_low_list))
upper_errors = np.maximum(0, np.array(p_high_list))
xerr = [lower_errors, upper_errors]

# Save model agreement
outpath = "/home/users/samuel.sancar/Heat_Extremes/Model_agreement/"
filename = "ERA5_agreement.npy"
full_path = os.path.join(outpath, filename)
os.makedirs(outpath, exist_ok=True)
np.save(full_path, fraction_gridpoints)

# Plotting threshold plot
plt.figure(figsize=(12, 7))
plt.errorbar(thresholds * 100, fraction_gridpoints, xerr=xerr, fmt='.',
             color='darkblue', ecolor='salmon', elinewidth=2, capsize=5,
             label='Model Agreement')
plt.plot(thresholds * 100, y, linestyle='--', color='gray', linewidth=2, label='1:1 Line')
plt.xlabel('Percentage of Models Less Than ERA5 at Grid Point (%)', fontsize=14, weight='bold')
plt.ylabel('Percentage of Land Grid Points (%)', fontsize=14, weight='bold')
plt.title('Model Agreement with ERA5 Across Land Grid Points', fontsize=16, weight='bold')
plt.grid(True, linestyle='--', alpha=0.6)
plt.legend(fontsize=12)
plt.tight_layout()
plt.show()

# Plot model agreement as a map
lats = template.coord('latitude').points
lons = template.coord('longitude').points
lon2d, lat2d = np.meshgrid(lons, lats)
agreement_plot = np.where(lat2d >= -60, agreement_land, np.nan)

fig = plt.figure(figsize=(14, 7))
ax = plt.axes(projection=ccrs.PlateCarree())
cmap = mpl_cm.get_cmap('coolwarm', 18)

c = ax.pcolormesh(lon2d, lat2d, agreement_plot, cmap=cmap, vmin=0, vmax=1, shading='auto')
ax.coastlines(resolution='50m', linewidth=1)
ax.add_feature(cfeature.LAND, facecolor='none', edgecolor='black')
ax.add_feature(cfeature.OCEAN, facecolor='lightblue')
plt.colorbar(c, ax=ax, orientation='vertical', label='Fraction of Models Less Than ERA5')
plt.title('Model Agreement: Fraction of Models Less Than ERA5 (per grid point)', fontsize=16, weight='bold')
plt.tight_layout()
plt.show()

