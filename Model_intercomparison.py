"""

 Code to calculate and plot model agreement intercomparison.
 
 To run: python Model_intercomparison.
 
"""

import iris
import os
import numpy as np
import matplotlib.pyplot as plt
import glob
from Lat_Land_mask_function import lat_land_mask
import matplotlib.cm as mpl_cm
import cartopy.crs as ccrs
import cartopy.feature as cfeature


# Load regridding template
template = iris.load_cube('template_cube_1deg.nc')

# Load ERA5 agreement plot
ERA5 = np.load('/home/users/samuel.sancar/Documents/Heat_Extremes/Model_agreement/ERA5_agreement.npy')

# Load ERA5 trends
datadir = '/data/users/samuel.sancar/Heat_extremes/trends/Reanalysis_trends/'
era5 = iris.load_cube(f"{datadir}era5_daily_maximum_2m_temperature_tailwidth_trends_1958-2023.nc")

# Get list of filenames
coupled_filenames = sorted(glob.glob('/data/users/samuel.sancar/Heat_extremes/trends/Coupled_model_trends/*'))
sst_filenames = sorted(glob.glob('/data/users/samuel.sancar/Heat_extremes/trends/SST_model_trends/*'))

# Load cubes and pair with filenames
coupled_files = [(iris.load_cube(fname), fname) for fname in coupled_filenames]
sst_files = [(iris.load_cube(fname), fname) for fname in sst_filenames]

# Regrid data
coupled_regridded = [(cube.regrid(template, iris.analysis.Linear()), fname) for cube, fname in coupled_files]
sst_regridded = [(cube.regrid(template, iris.analysis.Linear()), fname) for cube, fname in sst_files]
ERA5_regrided = era5.regrid(template, iris.analysis.Linear())

# Apply mask to ERA5 data (if needed)
era5_masked_cube, _ = lat_land_mask(ERA5_regrided.data, -60)
era5_masked = era5_masked_cube.data

plt.figure(figsize=(12, 8))

# Looping across Coupled model agreement
for i, (coupled, fname) in enumerate(coupled_regridded):
    indices = [j for j in range(len(coupled_regridded)) if j != i]
    model_stack = np.concatenate(
        [np.array([coupled_regridded[j][0].data for j in indices]),
         np.array([cube[0].data for cube in sst_regridded]),
         ERA5_regrided.data[np.newaxis, :, :]], axis=0)

    # Apply latitude and land mask to all models and ERA5
    agreement = np.mean(model_stack <= coupled.data, axis=0)
    agreement_land_cube, n_gridpoints = lat_land_mask(agreement, -60)
    agreement_land = agreement_land_cube.data

    # Define threshold interval
    thresholds = np.arange(0, 1.01, 0.01)
    fraction_gridpoints = []


    # Calculate fraction of gridpoints under threshold
    for thresh in thresholds:
        n = np.count_nonzero((agreement_land >= thresh) & ~np.isnan(agreement_land))
        fraction_gridpoints.append(100 * n / n_gridpoints)

    # Plot 1:1 line
    y = -np.arange(0, 101, 1) + 100
    plt.plot(thresholds * 100, fraction_gridpoints, marker='.', color='grey')  # Add label if needed

plt.plot(np.arange(0, 101, 1), y, linestyle='--', color='black', label='1:1 Line')
plt.plot(np.arange(0, 101, 1), ERA5, marker='.', color='red', label='ERA5 agreement')
plt.xlabel('X % of Models Less Than X at Grid Point')
plt.ylabel('Percentage of Land Grid Points (%)')
plt.title('Model Agreement Curves for ERA5 and Coupled models')
plt.grid(True)
plt.legend(loc='best', fontsize='small')
plt.tight_layout()
plt.show()

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
    agreement_land_cube, n_gridpoints = lat_land_mask(agreement, -60)
    agreement_land = agreement_land_cube.data

    # Define threshold interval
    thresholds = np.arange(0, 1.01, 0.01)
    fraction_gridpoints = []

    # Calculate fraction of gridpoints under threshold
    for thresh in thresholds:
        n = np.count_nonzero((agreement_land >= thresh) & ~np.isnan(agreement_land))
        fraction_gridpoints.append(100 * n / n_gridpoints)

    # Plot 1:1 line
    y = -np.arange(0, 101, 1) + 100
    plt.plot(thresholds * 100, fraction_gridpoints, marker='d', color='grey')  # Add label if needed

plt.plot(np.arange(0, 101, 1), y, linestyle='--', color='black', label='1:1 Line')
plt.plot(np.arange(0, 101, 1), ERA5, marker='.', color='red', label='ERA5 agreement')
plt.xlabel('X % of Models Less Than X at Grid Point')
plt.ylabel('Percentage of Land Grid Points (%)')
plt.title('Model Agreement Curves for ERA5 and SST Models')
plt.grid(True)
plt.legend(loc='best', fontsize='small')
plt.tight_layout()
plt.show()



