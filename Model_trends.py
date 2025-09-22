"""
28/07/25
---------

Example code to calculate trends in tail widths (p99-p87.5) for daily tasmax (using pearson r coefficient), for 
each model year 1958-2022 (Coupled models).

To run: python Model_trends.py [model] [experiment] [member]

"""

import iris
import numpy as np
import sys
import iris.analysis.stats
import matplotlib.pyplot as plt
import iris.plot as iplt
import iris.quickplot as qplt
from Trend_calc_function import calculate_trend

#Set model, experiment, ensemble member:
Mode = sys.argv[1]  # 'Coupled' or 'SST'
model = sys.argv[2]
experiment = sys.argv[3]
member = sys.argv[4]

# Define data file paths
datadir = '/data/users/samuel.sancar/Heat_extremes/Tail_widths/'
outdir = f'/data/users/samuel.sancar/Heat_extremes/Trends/{Mode}_model_trends/'

fn = model+'_'+experiment+'_'+member+'_tailwidths.nc'
yearly_tail_widths = iris.load_cube(datadir+fn)

#Creat 'time' cube:
year_cube = yearly_tail_widths.copy()
year_cube_data = np.empty(shape=year_cube.data.shape)
for i,y in enumerate(np.arange(1958,2023)):
    #print(y)
    year_cube_data[i,:,:] = y
year_cube.data = year_cube_data

# Extract the time coordinate
time_coord = yearly_tail_widths.coord('time')

# Print the raw time points
print("Raw time points:", time_coord.points)

# Convert to datetime objects
time_datetimes = time_coord.units.num2date(time_coord.points)
print("Datetime values:")
for dt in time_datetimes:
    print(dt)

# Calculate trend using function
trend = calculate_trend(year_cube, yearly_tail_widths)

iris.save(trend,outdir+model+'_'+experiment+'_'+member+'_tailwidth_trend.nc')
