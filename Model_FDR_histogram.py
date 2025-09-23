"""
23/09/2025
-----------

Histogram of FDR correction for tailwidth trends for all models, ERA5 and JRA-55 reanalysis data.

To run: python Model_FDR_histogram.py

"""

import iris
import numpy as np
import sys
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
import glob
import os
from Lat_Land_mask_function import lat_land_mask

variable = 'daily_maximum_2m_temperature'
reanalysis = 'era5'

# Load template cube
template_cube = iris.load_cube('template_cube_1deg.nc')

# Load and regrid ERA5 trend data
datadir = '/data/users/samuel.sancar/Heat_extremes/Trends/Reanalysis_trends/'
era5_cube = np.load(f"{datadir}{reanalysis}_{variable}_tailwidth_fdr_1958-2023.npy")

datadir_jra = '/data/users/samuel.sancar/Heat_extremes/Trends/Reanalysis_trends/'
jra_cube = np.load(f"{datadir_jra}jra55__minmax_surf.015_tmax.reg_tl319._tailwidth_fdr_1958-2023.npy")

# Prepare data fro ERA5
era5_masked_cube = era5_cube
era5_data = era5_masked_cube[~np.isnan(era5_masked_cube)]

# Prepare data for JRA-55
jra_masked_cube = jra_cube
jra_data = jra_masked_cube[~np.isnan(jra_masked_cube)]

# Load and regrid model cubes
coupled_filenames = sorted(glob.glob('/data/users/samuel.sancar/Heat_extremes/Trends/Trend_significance/FDR_corrections/*_hist-1950_*_fdr_1958-2023.npy'))
sst_filenames = sorted(glob.glob('/data/users/samuel.sancar/Heat_extremes/Trends/Trend_significance/FDR_corrections/*_highresSST-present_*_fdr_1958-2023.npy'))

coupled_data = [np.load(f) for f in coupled_filenames]
sst_data = [np.load(f) for f in sst_filenames]

significant_gridpoints = []

for filename in coupled_filenames:
    cube = np.load(filename)
    model_name = os.path.basename(filename)
    #sig_gridpoints = cube.data < 0.05
    sum_sig_gridpoints = np.count_nonzero(cube)
    print(f"Significant grid points in {model_name}: {sum_sig_gridpoints}")
    significant_gridpoints.append(sum_sig_gridpoints)

for filename in sst_filenames:
    cube = np.load(filename)
    model_name = os.path.basename(filename)
    #sig_gridpoints = cube.data < 0.05
    sum_sig_gridpoints = np.count_nonzero(cube)
    print(f"Significant grid points in {model_name}: {sum_sig_gridpoints}")
    significant_gridpoints.append(sum_sig_gridpoints)


#era5_sig_gridpoints = era5_data.data < 0.05
sum_era5_sig_gridpoints = np.count_nonzero(era5_data)
significant_gridpoints.append(sum_era5_sig_gridpoints)

#jra_sig_gridpoints = jra_cube.data < 0.05
sum_jra_sig_gridpoints = np.count_nonzero(jra_data)
significant_gridpoints.append(sum_jra_sig_gridpoints)

plt.figure(figsize=(12, 7))
counts, bins, patches = plt.hist(significant_gridpoints,bins=30,range=(0, 1200),histtype='bar',color='#4682B4',  # Steel bluealpha=0.85,edgecolor='white',
    linewidth=2)

# Label histogram at ERA5 value
plt.axvline(x=sum_era5_sig_gridpoints, color='#D7263D', linestyle='--', linewidth=2.5, label='ERA5 Value')
plt.text(
    sum_era5_sig_gridpoints + 50, max(counts) * 0.92, 'ERA5',
    color='#D7263D', fontsize=13, fontweight='bold',
    bbox=dict(facecolor='white', edgecolor='#D7263D', boxstyle='round,pad=0.3', alpha=0.8)
)

# Label histogram at JRA-55 value
plt.axvline(x=sum_jra_sig_gridpoints, color='#1B998B', linestyle='--', linewidth=2.5, label='JRA-55 Value')
plt.text(
    sum_jra_sig_gridpoints + 50, max(counts) * 0.82, 'JRA-55',
    color='#1B998B', fontsize=13, fontweight='bold',
    bbox=dict(facecolor='white', edgecolor='#1B998B', boxstyle='round,pad=0.3', alpha=0.8)
)

# title
plt.xlabel('Number of FDR corrected significant grid points', fontsize=15, fontweight='bold')
plt.ylabel('Number of Models', fontsize=15, fontweight='bold')
plt.title('Histogram of Significant Grid Points in Tail Width Trends (1958–2022)', fontsize=18, fontweight='bold')

# Ticks and grid
plt.xticks(fontsize=13)
plt.yticks(fontsize=13)
plt.grid(True, linestyle='--', alpha=0.5, linewidth=1.2)

# Legend
plt.legend(fontsize=13, loc='upper right', frameon=False)
plt.tight_layout()
plt.show()



