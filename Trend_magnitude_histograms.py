"""
17/07/25
---------
Plot histogram in tail widths trends of ERA5 data.
To run: python plot_tailwidth_trend.py [1958] [2023]

"""

import iris
import sys
import matplotlib.pyplot as plt
import numpy as np
import glob
from Lat_Land_mask_function import lat_land_mask
from Save_figure_function import save_figure

# Read command-line arguments
reanalysis = 'era5'
t1 = int(sys.argv[1])
t2 = int(sys.argv[2])
variable = 'daily_maximum_2m_temperature'

jra = 'jra55'
variable_jra = '_minmax_surf.015_tmax.reg_tl319.'

# Load ERA5 trend data
datadir = '/data/users/samuel.sancar/Heat_extremes/Trends/Reanalysis_trends/'
trend = iris.load_cube(f"{datadir}{reanalysis}_{variable}_tailwidth_trends_{t1}-{t2}.nc")
masked_trend, _ = lat_land_mask(trend, -60, 90)
era5_data = masked_trend.data.flatten()
era5_data = era5_data[~np.isnan(era5_data)]

# Load JRA-55 trend data
datadir_jra = '/data/users/samuel.sancar/Heat_extremes/Trends/Reanalysis_trends/'
trend_jra = iris.load_cube(f"{datadir_jra}{jra}_{variable_jra}_tailwidth_trends_{t1}-{t2}.nc")
masked_trend_jra, _ = lat_land_mask(trend_jra, -60, 90)
jra_data = masked_trend_jra.data.flatten()
jra_data = jra_data[~np.isnan(jra_data)]

# Load model data
coupled_files = sorted(glob.glob('/data/users/samuel.sancar/Heat_extremes/Trends/Coupled_model_trends/*_tailwidth_trend.nc'))
sst_files = sorted(glob.glob('/data/users/samuel.sancar/Heat_extremes/Trends/SST_model_trends/*_tailwidth_trend.nc'))
coupled_data = [iris.load_cube(f) for f in coupled_files]
sst_data = [iris.load_cube(f) for f in sst_files]

# Define bins
bins = np.arange(-0.5, 1.5, 0.1)

# --- Coupled Models Histogram ---
plt.figure(figsize=(14, 5))
for cube in coupled_data:
    masked_model, _ = lat_land_mask(cube, -60, 90)
    model_data = masked_model.data.flatten()
    model_data = model_data[~np.isnan(model_data)]
    plt.hist(model_data, bins=bins, histtype='stepfilled', alpha=0.2, color='skyblue', edgecolor='black')

plt.hist(era5_data, bins=bins, histtype='step', color='red', linewidth=2.5, linestyle='--', label='ERA5')
plt.hist(jra_data, bins=bins, histtype='step', color='blue', linewidth=2.5, linestyle='-.', label='JRA-55')
plt.xlabel('Trend Magnitude (°C per Decade)', fontsize=13, fontweight='bold')
plt.ylabel('Number of Grid Points', fontsize=13, fontweight='bold')
plt.title('Tail Width Trend Histogram (Coupled Models vs ERA5, 1958–2022)', fontsize=15, fontweight='bold')
plt.legend(loc='upper right', fontsize=12)
plt.grid(True, linestyle='--', alpha=0.6)
plt.tight_layout()
plt.show()

# --- SST Models Histogram ---
plt.figure(figsize=(14, 12))
for cube in sst_data:
    masked_model, _ = lat_land_mask(cube, -60, 90)
    model_data = masked_model.data.flatten()
    model_data = model_data[~np.isnan(model_data)]
    plt.hist(model_data, bins=bins, range = (-0.5, 1.5), histtype='stepfilled', alpha=0.2, color='lightgreen', edgecolor='black')
for cube in coupled_data:
    masked_model, _ = lat_land_mask(cube, -60, 90)
    model_data = masked_model.data.flatten()
    model_data = model_data[~np.isnan(model_data)]
    plt.hist(model_data, bins=bins, range = (-0.5, 1.5),histtype='stepfilled', alpha=0.2, color='skyblue', edgecolor='black')
plt.hist(era5_data, bins=bins, range = (-0.5, 1.5), histtype='step', color='red', linewidth=2.5, linestyle='--', label='ERA5')
plt.hist(jra_data, bins=bins, range = (-0.5, 1.5), histtype='step', color='blue', linewidth=2.5, linestyle='-.', label='JRA-55')
plt.xlabel('Trend Magnitude (°C per Decade)', fontsize=20, fontweight='bold')
plt.xticks(fontsize = 15)
plt.ylabel('Number of Grid Points', fontsize=15, fontweight='bold')
plt.yticks(fontsize = 15)

# Set x-axis limit
plt.xlim(-0.5, 1.5)

plt.title('Tail Width Trend Histogram for Models, ERA5 and JRA-55', fontsize=20, fontweight='bold')
plt.legend(loc='upper right', fontsize=15)
plt.grid(True, linestyle='--', alpha=0.6)
#plt.tight_layout()
#save_figure(plt, 'Trend_magnitude_histogram', 'Graphs')
plt.show()



