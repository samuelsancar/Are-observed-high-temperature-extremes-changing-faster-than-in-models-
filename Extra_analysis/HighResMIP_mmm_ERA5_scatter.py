"""
14/05/26
--------

Code to plot scatter plot figure of ERA5 tail width trends vs HighResMIP multi-model mean tail width trends.
Demonstrates that ERA5 heatwave hotspots are not fully coincident with a forced signal in models.

To run: HighResMIP_mmm_ERA5_scatter.py
"""

import iris
import sys
import matplotlib.pyplot as plt
import matplotlib.cm as mpl_cm
import numpy as np
from Lat_Land_mask_function import lat_land_mask
from Save_figure_function import save_figure


cmap_rb = mpl_cm.get_cmap('coolwarm', 15)

t1 = 1958
t2= 2023

alpha = 0.05

lat_minimum = -60
lat_maximum = 60

datadir = f'/data/users/samuel.sancar/Heat_extremes/Trends/'

reanalysis = 'era5'
variable = 'daily_maximum_2m_temperature'
era5_datadir = '/data/users/samuel.sancar/Heat_extremes/Trends/Reanalysis_trends/'

template = iris.load_cube('template_cube_1deg.nc')
era5_trend = iris.load_cube(f"{era5_datadir}{reanalysis}_{variable}_tailwidth_trends_{t1}-{t2-1}.nc")
era5_p87p5_trend = iris.load_cube(f"{era5_datadir}{reanalysis}_{variable}_p87p5_trends_{t1}-{t2-1}.nc")

filename = f"HighResMIP_MMM_tailwidth_trend_1958-2022.nc"
cube = iris.load_cube(datadir + filename)
p87p5_trend = iris.load_cube(datadir + filename.replace('tailwidth','p87p5'))

# Guess bounds
lon_coord = template.coord('longitude')
if not lon_coord.has_bounds():
    lon_coord.guess_bounds()

lat_coord = template.coord('latitude')
if not lat_coord.has_bounds():
    lat_coord.guess_bounds()

# Apply land mask
masked_cube, _ = lat_land_mask(cube, template, lat_minimum, lat_maximum)
p87p5_masked, _ = lat_land_mask(p87p5_trend, template, lat_minimum, lat_maximum)

# Apply latitude and land mask using function
era5_model_mask, _ = lat_land_mask(era5_trend, template, lat_minimum, lat_maximum)
era5_p87p5_masked, _ = lat_land_mask(era5_p87p5_trend, template, lat_minimum, lat_maximum)

print(masked_cube.coord('latitude'))
print(era5_model_mask.coord('latitude'))

print(masked_cube.coord('longitude'))
print(era5_model_mask.coord('longitude'))

x = masked_cube.data
y = era5_model_mask.data

x = masked_cube.data.filled(np.nan)
y = era5_model_mask.data.filled(np.nan)

# Filter out unphysical numbers arising from masking issues.
mask = (
    np.isfinite(x) &
    np.isfinite(y) &
    (np.abs(x) < 5) &
    (np.abs(y) < 5)
)

x_valid = x[mask]
y_valid = y[mask]

# Correlation
r = np.corrcoef(x_valid, y_valid)[0,1]
print(f"Correlation: {r:.2f}")

# Linear regression (slope + intercept)
slope, intercept = np.polyfit(x_valid, y_valid, 1)
print(f"Slope: {slope:.2f}, Intercept: {intercept:.2f}")

plt.figure(figsize=(12, 6))
ax = plt.gca()
for spine in ax.spines.values():
    spine.set_linewidth(2)
    
plt.scatter(x_valid, y_valid, s=5)

# Line of best fit
x_fit = np.linspace(x_valid.min(), x_valid.max(), 100)
y_fit = slope * x_fit + intercept

plt.plot(x_fit, y_fit, color='black', linewidth=2, label=f'Fit: y={slope:.2f}x + {intercept:.2f}')

plt.axhline(0, color='grey', linewidth=1)
plt.axvline(0, color='grey', linewidth=1)
plt.scatter(0, 0, color='red', s=40, zorder=5, label='(0,0)')

plt.title('c) HRMIP MMM and ERA5 tail width trends', fontsize=22, fontweight='bold')
plt.xlabel('HighResMIP MMM tail width trend (°C/decade)', fontsize=21)
plt.ylabel('ERA5 tail width trend (°C/decade)', fontsize=21)

plt.xticks(fontsize=19)
plt.yticks(fontsize=19)
plt.grid(True, linestyle='--', alpha=0.5, linewidth=1.2)
plt.legend(fontsize=15)
save_figure(plt, 'c) HighResMIP_MMM_ERA5_scatter_plot', 'Graphs')
plt.show()