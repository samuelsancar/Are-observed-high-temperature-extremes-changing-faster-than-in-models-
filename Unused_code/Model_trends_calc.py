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
from Trend_calculation_functions import calculate_trend

#Set model, experiment, ensemble member:
Mode = sys.argv[1]  # 'Coupled' or 'SST'
model = sys.argv[2]
experiment = sys.argv[3]
member = sys.argv[4]

# Define data file paths
datadir = '/data/users/samuel.sancar/Heat_extremes/Tail_widths/'
p87p5 = '/data/users/samuel.sancar/Heat_extremes/Tail_widths/'
outdir = f'/data/users/samuel.sancar/Heat_extremes/Trends/{Mode}_model_trends/'

fn = model+'_'+experiment+'_'+member+'_tailwidths.nc'
yearly_tail_widths = iris.load_cube(datadir+fn)
yearly_p87p5 = iris.load_cube(p87p5+fn.replace('tailwidths','p87.5'))

#Creat 'time' cube:
year_cube = yearly_tail_widths.copy()
year_cube_data = np.empty(shape=year_cube.data.shape)
for i,y in enumerate(np.arange(1958,2023)):
    #print(y)
    year_cube_data[i,:,:] = y
year_cube.data = year_cube_data

year_p87p5 = yearly_p87p5.copy()
year_p87p5_cube = np.empty(shape=year_p87p5.data.shape)
for i,y in enumerate(np.arange(1958,2023)):
    #print(y)
    year_p87p5_cube[i,:,:] = y
year_p87p5.data = year_p87p5_cube

# Extract the time coordinate
time_coord = yearly_tail_widths.coord('time')
p87p5_time_coord = yearly_p87p5.coord('time')

# Print the raw time points
print("Raw time points:", time_coord.points)

# Convert to datetime objects
time_datetimes = time_coord.units.num2date(time_coord.points)
print("Datetime values:")
for dt in time_datetimes:
    print(dt)
    
# Convert to datetime objects
p87p5_time_datetimes = p87p5_time_coord.units.num2date(p87p5_time_coord.points)
print("Datetime values:")
for dt in p87p5_time_datetimes:
    print(dt)

# Calculate trend using function
trend = calculate_trend(year_cube, yearly_tail_widths)
p87p5_trend = calculate_trend(year_p87p5, yearly_p87p5)

iris.save(trend,outdir+model+'_'+experiment+'_'+member+'_tailwidth_trend.nc')
iris.save(trend,outdir+model+'_'+experiment+'_'+member+'_p87p5_trend.nc')
