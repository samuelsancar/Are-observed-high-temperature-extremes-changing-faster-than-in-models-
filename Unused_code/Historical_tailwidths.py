"""
25/07/25

Example code to calculate tail width (p99 - p87.5) in daily tasmax,
for each model year 1958–2014 (for models with only historical runs).

To run: python Hist_tailwidths.py [HadGEM3-GC31-??] [hist-1950] [r??] 
"""

import iris
import numpy as np
import sys
import iris.analysis.stats
import matplotlib.pyplot as plt
import iris.plot as iplt
import iris.quickplot as qplt
from Trend_calculation_functions import calculate_trend

# Command-line arguments
model = sys.argv[1]
experiment = sys.argv[2]
member = sys.argv[3]

# Define directories
datadir = f'/data/users/samuel.sancar/Climate_models/hist/{Group}/'
outdir = '/data/users/samuel.sancar/Heat_extremes/trends/'

# Load template grid for regridding
template = iris.load_cube('template_cube_1deg.nc')

# CubeList to store tail widths
tail_width_list = iris.cube.CubeList([])

# Load and process each year's data
subdir = f'/tasmax_day_{model}_{experiment}_{member}/'
print(subdir)

for y in range(1950, 2015):
    fn = f'tasmax_day_{model}_{experiment}_{member}_*_{y}0101-{y}12*.nc'
    cube = iris.load_cube(datadir + subdir + fn)

    # Calculate tail width (p99 - p87.5)
    p99 = cube.collapsed('time', iris.analysis.PERCENTILE, percent=99)
    p87p5 = cube.collapsed('time', iris.analysis.PERCENTILE, percent=87.5)
    tail_width = p99 - p87p5

    # Regrid to template
    tail_width_rg = tail_width.regrid(template, iris.analysis.Linear())
    tail_width_list.append(tail_width_rg)

# Merge into a single cube
iris.util.equalise_attributes(tail_width_list)
yearly_tail_widths = tail_width_list.merge()[0]

# Create a 'year' cube with the same shape
year_cube = yearly_tail_widths.copy()
year_cube_data = np.empty_like(year_cube.data)

n_years = year_cube_data.shape[0]
start_year = 1950

for i in range(n_years):
    year_cube_data[i, :, :] = start_year + i

year_cube.data = year_cube_data

# Calculate Pearson correlation and trend
trend = calculate_trend(year_cube, yearly_tail_widths)

iris.save(trend, f'{outdir}{model}_{experiment}_{member}_tailwidth_trend.nc')
