"""
17/07/25
---------
Code to merge SST yearly 'daily tasmax tailwidth' files and merge for years 1958-2022.

To run: python SST_merge&trend.py

"""

import iris
import numpy as np
import sys
import iris.analysis.stats
import iris.plot as iplt
import iris.quickplot as qplt
import os
import cf_units
from Trend_calc_function import calculate_trend

# Define data file paths
datadir = '/data/users/samuel.sancar/Heat_extremes/tailwidthslist/'
outdir = '/data/users/samuel.sancar/Heat_extremes/Trends/SST_model_trends/'

# Set data variables
Group = sys.argv[1]
model = sys.argv[2]
experiment = sys.argv[3]
member = sys.argv[4]
t1 = 1958
t2 = 2023

os.makedirs(outdir, exist_ok=True)

# Define timeline
years = np.arange(t1, t2)
time_unit = cf_units.Unit('days since 1970-01-01', calendar='gregorian')

# Create a time coordinate
time_coord = iris.coords.DimCoord(
    years,
    standard_name='time',
    long_name='valid_time',
    units=time_unit,
)

# Create list of year cubes
yearly_tail_widths = iris.cube.CubeList([])

# Combine Year cubes into list
for y in np.arange(t1, t2):
    fn = Group+'/'+model+'_'+experiment+'_'+member+'_'+'tailwidth'+'_'+str(y)+'.nc'
    year_cube = iris.load_cube(datadir+fn)

    # Remove singleton time dimension
    year_cube = year_cube.slices(['latitude', 'longitude']).next()
    yearly_tail_widths.append(year_cube)


# Use spatial coordinates from one of the cubes
lat_coord = yearly_tail_widths[0].coord('latitude')
lon_coord = yearly_tail_widths[0].coord('longitude')

# Equalise attributes
iris.util.equalise_attributes(yearly_tail_widths)
data_stack = np.stack([cube.data for cube in yearly_tail_widths], axis=0)

# Create the 3D cube
yearly_tail_widths = iris.cube.Cube(
    data_stack,
    dim_coords_and_dims=[(time_coord, 0), (lat_coord, 1), (lon_coord, 2)],
    long_name='tail_width_over_years',
    units=yearly_tail_widths[0].units,
    attributes=yearly_tail_widths[0].attributes
)

yearly = yearly_tail_widths.copy()
year_cube_data = np.empty(shape=yearly.data.shape)

# Define each year
for i,y in enumerate(np.arange(t1, t2)):
    #print(y)
    year_cube_data[i, :, :] = y

yearly.data = year_cube_data

# Calculate trend using function
trend = calculate_trend(yearly, yearly_tail_widths)

out_file_trend = os.path.join(outdir, f"{model}_{experiment}_{member}_tailwidth_trend.nc")
iris.save(trend,out_file_trend)