"""
20/08/25
---------
Code to demonstrate weighted Model agreement with ERA5 for daily tasmax tailwidths and the associated sampling error. Output: 2 plots.

To run: python Model_Agreement.py

"""

import iris
import numpy as np
import matplotlib.pyplot as plt
import scipy.special
import os

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

# Average across models
agreement = np.mean(model_stack <= era5_regridded.data, axis=0)*era5_regridded.data#/(np.max(era5_regridded.data)-np.min(era5_regridded)) 

# Calculating sampling error using Binomial_uncertainty:
#from BnD_sampling_errors import binomial_uncertainty

# Make error lists
#p_low_list = []
#p_high_list = []

# Compute errors by looping over percentages
#for i in np.arange(0,101,1):
    #if i !=0:
        #k = (32/100)*i
    #else:
        #k = np.inf
    #p_low, p_high, _,_,_= binomial_uncertainty(32, k)

    # Error propagation and appending to list
    #p_high_list.append(100/32 * p_high)
    #p_low_list.append(100/32 * p_low)

# Only consider land points
agreement_land = np.where(land_mask, agreement_magnitude_weighted, np.nan)

# Scan over thresholds (0 to 1, step 0.01)
thresholds = np.arange(0, 1.01, 0.01)
n_gridpoints = np.count_nonzero(~np.isnan(agreement_land))

fraction_gridpoints = [] #Create list to store fractions

# Loop over thresholds to count number of gridpoints
for thresh in thresholds:
    grid_points = (agreement_land >= thresh) & (~np.isnan(agreement_land))
    k = np.count_nonzero(grid_points)
    fraction_gridpoints.append(100 * k / n_gridpoints)
 
# 1:1 line
y = - np.arange(0, 101, 1) + 100

# Compute error margins
#lower_errors = np.maximum(0, np.array(p_low_list))
#upper_errors = np.maximum(0, np.array(p_high_list))
#xerr = [lower_errors, upper_errors]

# Save model agreement
outpath = "/home/users/samuel.sancar/Heat_Extremes/Model_agreement/"
filename = "Weighted_ERA5_agreement.npy"
full_path = os.path.join(outpath, filename)

os.makedirs(outpath, exist_ok=True)

np.save(full_path, fraction_gridpoints)

# Plotting threshold plot
plt.figure(figsize=(10, 6))
plt.plot(thresholds * 100, fraction_gridpoints, marker='.', color = 'red')
plt.plot(thresholds * 100, y, linestyle='--', color='gray', label='1:1 Line')
plt.xlabel('X% of Models Less Than ERA5 at Grid Point')
plt.ylabel('Percentage of Land Grid Points (%)')
plt.title('Grid points with at least X% of models less than ERA5')
plt.grid(True)
plt.tight_layout()
plt.show()

