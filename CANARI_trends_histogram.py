"""
04/11/25
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
CANARI_files = sorted(glob.glob('/data/users/samuel.sancar/Heat_extremes/Trends/CANARI_trends/*_tailwidth_trend.nc'))
CANARI_data = [iris.load_cube(f) for f in CANARI_files]

# Define bins
bins = np.arange(-0.5, 1.5, 0.1)

plt.figure(figsize=(14, 12))
for cube in CANARI_data:
    masked_model, _ = lat_land_mask(cube, -60, 90)
    model_data = masked_model.data.flatten()
    model_data = model_data[~np.isnan(model_data)]
    plt.hist(model_data, bins=bins, histtype='stepfilled', alpha=0.2, color='skyblue', edgecolor='black')

plt.hist(era5_data, bins=bins, histtype='step', color='red', linewidth=2.5, linestyle='--', label='ERA5')
plt.hist(jra_data, bins=bins, histtype='step', color='blue', linewidth=2.5, linestyle='-.', label='JRA-55')
plt.xlabel('Trend Magnitude (°C per Decade)', fontsize=10, fontweight='bold')
plt.ylabel('Number of Grid Points', fontsize=10, fontweight='bold')
plt.title('Histogram of Tail-Width Trends (CANARI LE, ERA5 & JRA-55, 1958–2022)', fontsize=15, fontweight='bold')
plt.legend(loc='upper right', fontsize=10)
plt.xlim(-0.5, 1.5)
plt.grid(True, linestyle='--', alpha=0.6)
plt.tight_layout()
plt.show()

