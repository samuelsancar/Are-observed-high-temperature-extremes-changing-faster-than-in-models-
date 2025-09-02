"""

26/08/25

Plot histogram for significant tail widths in Models and ERA5 trends. Land masks have already been applied to models.

To run: python plot_tailwidth_trend.py [1958] [2023]

"""

import iris
import sys
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
import glob
from Lat_Land_mask_function import lat_land_mask


variable = 'daily_maximum_2m_temperature'
reanalysis = 'era5'

# Load template cube
template_cube = iris.load_cube('template_cube_1deg.nc')

# Load and regrid ERA5 trend data
datadir = '/data/users/samuel.sancar/Heat_extremes/trends/Reanalysis_trends/'
era5_cube = iris.load_cube(f"{datadir}{reanalysis}_{variable}_tailwidth_trends_{1958}-{2023}.nc")
era5_regridded = era5_cube.regrid(template_cube, iris.analysis.Linear())

# Apply land and latitude mask to ERA5
era5_masked_cube, _ = lat_land_mask(era5_regridded, -60)
era5_data = era5_masked_cube.data.flatten()
era5_data = era5_data[~np.isnan(era5_data)]

# Load and regrid model cubes
coupled_filenames = sorted(glob.glob('/data/users/samuel.sancar/Heat_extremes/trends/Trend_significance/*_hist-1950_*'))
sst_filenames = sorted(glob.glob('/data/users/samuel.sancar/Heat_extremes/trends/Trend_significance/*_highresSST-present_*'))

coupled_data = [iris.load_cube(f).regrid(template_cube, iris.analysis.Linear()) for f in coupled_filenames]
sst_data = [iris.load_cube(f).regrid(template_cube, iris.analysis.Linear()) for f in sst_filenames]

significant_gridpoints = []

for f in coupled_filenames:
    sig_gridpoints = f.data > 0.05
    sum_sig_gridpoints = np.sum(sig_gridpoints)
    significant_gridpoints.append(sum_sig_gridpoints)

for f in sst_filenames:
    sig_gridpoints = f.data > 0.05
    sum_sig_gridpoints = np.sum(sig_gridpoints)
    significant_gridpoints.append(sum_sig_gridpoints)

era5_sig_gridpoints = era5_masked_cube.data > 0.05
sum_era5_sig_gridpoints = np.sum(era5_sig_gridpoints)
significant_gridpoints.append(sum_era5_sig_gridpoints)
    
plt.figure(figsize=(10, 6))
counts, bins, patches = plt.hist(significant_gridpoints, bins=33, color='blue', alpha=0.7)

# Find the bin that contains the value 60
bin_index = np.digitize(sum_era5_sig_gridpoints, bins) - 1
bin_center = (bins[bin_index] + bins[bin_index + 1]) / 2
frequency = counts[bin_index]

# Annotate the histogram at value 60
plt.text(bin_center, frequency, f'{frequency:.0f}', ha='center', va='bottom', color='red')
plt.axvline(x=60, color='red', linestyle='--', label='Value 60')
plt.legend()

plt.xlabel('Number of Significant Grid Points')
plt.ylabel('Number of Models')
plt.label()
plt.title(f'Histogram of Significant Grid Points in Tail Width Trends ({t1}–{t2})')
plt.grid(True)
plt.tight_layout()
plt.show()