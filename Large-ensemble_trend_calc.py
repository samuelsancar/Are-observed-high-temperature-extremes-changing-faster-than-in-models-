"""
    Code to merge tailwidth trends and calculate tailwidth trends for MASS large ensemble
   
"""

import iris
import numpy as np
import sys
import os
import cf_units
from iris.cube import CubeList
from Trend_calculation_function import calculate_trend
from ERA5_trend_calc import maskp87p5

# file paths
datadir = '/data/users/samuel.sancar/Heat_extremes/tailwidthslist/MASS/'
outtrend = '/data/users/samuel.sancar/Heat_extremes/Trends/MASS_model_trends/'
outtailwidth = '/data/users/samuel.sancar/Heat_extremes/Tail_widths/MASS/'

# data variables
member = sys.argv[1]
t1 = 1958
t2 = 2023

os.makedirs(outtrend, exist_ok=True)

# Initialise lists
yearly_tail_widths = iris.cube.CubeList([])
yearly_p87p5 = iris.cube.CubeList([])
valid_years_tail = []
valid_years_p87p5 = []

# Load tail width cubes
for y in np.arange(t1, t2):
    tail_fn = f"LE_r0{member}i1p00000_tailwidth_{y}.nc"
    tail_path = os.path.join(datadir, tail_fn)
    try:
        year_cube = iris.load_cube(tail_path)
        year_cube = year_cube.slices(['latitude', 'longitude']).next()
        yearly_tail_widths.append(year_cube)
        valid_years_tail.append(y)
    except Exception as e:
        print(f"Skipping year {y} for tail width: {e}")

# Load p87p5 cubes
for y in np.arange(t1, t2):
    p87p5_fn = f"LE_r0{member}i1p00000_p87p5_{y}.nc"
    p87p5_path = os.path.join(datadir, p87p5_fn)
    try:
        p87p5_year_cube = iris.load_cube(p87p5_path)
        p87p5_year_cube = p87p5_year_cube.slices(['latitude', 'longitude']).next()
        yearly_p87p5.append(p87p5_year_cube)
        valid_years_p87p5.append(y)
    except Exception as e:
        print(f"Skipping year {y} for p87.5 percentile: {e}")

# Equalise attributes
iris.util.equalise_attributes(yearly_tail_widths)
iris.util.equalise_attributes(yearly_p87p5)

# Stack data
data_stack_tail = np.stack([cube.data for cube in yearly_tail_widths], axis=0)
data_stack_p87p5 = np.stack([cube.data for cube in yearly_p87p5], axis=0)

# Create time coordinates
time_unit = cf_units.Unit('days since 1970-01-01', calendar='gregorian')
time_coord_tail = iris.coords.DimCoord(valid_years_tail, standard_name='time', long_name='valid_time', units=time_unit)
time_coord_p87p5 = iris.coords.DimCoord(valid_years_p87p5, standard_name='time', long_name='valid_time', units=time_unit)

# Get spatial coordinates
lat_coord_tail = yearly_tail_widths[0].coord('latitude')
lon_coord_tail = yearly_tail_widths[0].coord('longitude')
lat_coord_p87p5 = yearly_p87p5[0].coord('latitude')
lon_coord_p87p5 = yearly_p87p5[0].coord('longitude')

# Create cubes
cube_tail = iris.cube.Cube(
    data_stack_tail,
    dim_coords_and_dims=[(time_coord_tail, 0), (lat_coord_tail, 1), (lon_coord_tail, 2)],
    long_name='tail_width_over_years',
    units=yearly_tail_widths[0].units,
    attributes=yearly_tail_widths[0].attributes
)

cube_p87p5 = iris.cube.Cube(
    data_stack_p87p5,
    dim_coords_and_dims=[(time_coord_p87p5, 0), (lat_coord_p87p5, 1), (lon_coord_p87p5, 2)],
    long_name='p87p5_over_years',
    units=yearly_p87p5[0].units,
    attributes=yearly_p87p5[0].attributes
)

# Save tail width cube
out_file_tailwidth = os.path.join(outtailwidth, f"{member}_tailwidths.nc")
iris.save(cube_tail, out_file_tailwidth)

# Create year cubes for trend calculation
year_cube_tail = cube_tail.copy()
year_cube_tail.data = np.tile(np.array(valid_years_tail)[:, None, None], (1, cube_tail.shape[1], cube_tail.shape[2]))

year_cube_p87p5 = cube_p87p5.copy()
year_cube_p87p5.data = np.tile(np.array(valid_years_p87p5)[:, None, None], (1, cube_p87p5.shape[1], cube_p87p5.shape[2]))

# Calculate trends
trend_tail = calculate_trend(year_cube_tail, cube_tail)
trend_p87p5 = calculate_trend(year_cube_p87p5, cube_p87p5)
masked_trend = maskp87p5(trend_p87p5, trend_tail)

out_file_trend_tail = os.path.join(outtrend, f"{member}_tailwidth_trend.nc")
iris.save(masked_trend, out_file_trend_tail)

out_file_trend_p87p5 = os.path.join(outtrend, f"{member}_p87p5_trend.nc")
iris.save(trend_p87p5, out_file_trend_p87p5)
