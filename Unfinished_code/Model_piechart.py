"""

 Code to calculate and plot model agreement intercomparison.
 
 To run: python Model_intercomparison.
 
"""

import iris
import numpy as np
import matplotlib.pyplot as plt
import scipy.special
import glob
import os
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.cm as mpl_cm
from Lat_Land_mask_function import lat_land_mask
from Save_figure_function import save_figure

# Load template
template = iris.load_cube('template_cube_1deg.nc')

# Load model trends
coupled_files = iris.load('/data/users/samuel.sancar/Heat_extremes/Trends/Coupled_model_trends/*_tailwidth_trend.nc')
sst_files = iris.load('/data/users/samuel.sancar/Heat_extremes/Trends/SST_model_trends/*_tailwidth_trend.nc')

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

# Scan over thresholds
thresh = 100
fraction_gridpoints = []
fraction_gridpoints_jra = []


grid_points = (agreement_land >= thresh) & (~np.isnan(agreement_land))
grid_points_jra = (agreement_land_jra >= thresh) & (~np.isnan(agreement_land_jra))
k = np.count_nonzero(grid_points)
k_jra = np.count_nonzero(grid_points_jra)
# Use n_gridpoints from ERA5 mask for consistency
fraction_gridpoints.append(100 * k / n_gridpoints)
fraction_gridpoints_jra.append(100 * k_jra / n_gridpoints_jra)

# Load cubes and pair with filenames
coupled = [(iris.load_cube(fname), fname) for fname in coupled_files]
sst = [(iris.load_cube(fname), fname) for fname in sst_files]

# Regrid data
coupled_regridded = [(cube.regrid(template, iris.analysis.Linear()), fname) for cube, fname in coupled]
sst_regridded = [(cube.regrid(template, iris.analysis.Linear()), fname) for cube, fname in sst]
ERA5_regrided = era5.regrid(template, iris.analysis.Linear())

# Combine all models (coupled + SST) into one list for comparison
all_models = coupled_regridded + sst_regridded
model_names = []

# Extract model names for labels
for _, fname in all_models:
    filename = os.path.splitext(os.path.basename(fname))[0]
    parts = filename.split('_')
    model_name = f"{parts[0]}_{parts[2]}" if len(parts) > 2 else filename
    model_names.append(model_name)

# Stack all model data into a single array: shape (n_models, lat, lon)
model_data_stack = np.array([cube.data for cube, _ in all_models])

# Apply land/latitude mask to all models
masked_models = []
for data in model_data_stack:
    masked_cube, _ = lat_land_mask(data, -60, 90)
    masked_models.append(masked_cube.data)

masked_models = np.array(masked_models)

# For each grid point, find index of model with maximum trend
max_indices = np.nanargmax(masked_models, axis=0)

# Count occurrences for each model
counts = [np.count_nonzero(max_indices == i) for i in range(len(model_names))]

# Convert to percentages
total_points = np.sum(counts)
percentages = [100 * c / total_points for c in counts]

# Plot pie chart
plt.figure(figsize=(10, 8))
plt.pie(percentages, labels=model_names, autopct='%1.1f%%', startangle=140)
plt.title('Fraction of Land Grid Points Dominated by Each Model (Greatest Trend)', fontsize=14, weight='bold')
plt.tight_layout()
plt.show()
