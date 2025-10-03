"""
17/09/2025
----------

Code to download JRA-55 data, aggregate for daily maximum 2m temperature data, and find annual tail-width.

To run: python JRA55_tasmax.py

"""

import iris
import iris.analysis
import iris.coord_categorisation
import os

# Define directories and parameters
datadir = '/data/users/samuel.sancar/Climate_models/JRA-55/'
outdir = '/data/users/samuel.sancar/Heat_extremes/JRA-55_tailwidths/daily_maximum__2m_temperature/'
template = iris.load_cube('template_cube_1deg.nc')
reanalysis = 'jra55'
variable = 'minmax_surf.015_tmax.reg_tl319.'

os.makedirs(outdir, exist_ok=True)

# Process each year
for y in range(1958, 2023):

    # Load monthly data files for the year
    file_pattern = f"{datadir}{variable}{y}*"
    cube_list = iris.load(file_pattern)

    iris.util.equalise_attributes(cube_list)
    iris.util.unify_time_units(cube_list)

    merged_cube = cube_list.concatenate_cube()

    # Reshape to daily blocks (8 time steps per day)
    no_of_days = merged_cube.coord('time').shape[0] // 8
    merged_cube = merged_cube[:no_of_days * 8]

    # Add 'day' categorisation for aggregation
    iris.coord_categorisation.add_day_of_year(merged_cube, 'time', name='day')

    # Aggregate to daily maximum
    daily_cube = merged_cube.aggregated_by('day', iris.analysis.MAX)

    # Create a new time coordinate based on the aggregated cube's shape
    time_points = daily_cube.coord('day').points
    time_units = merged_cube.coord('time').units

    # Ensure time_points are strictly monotonic
    unique_sorted_time_points = sorted(set(time_points))
    new_time_coord = iris.coords.DimCoord(unique_sorted_time_points, standard_name='time', units=time_units)

    # Replace time coordinate in daily_cube
    daily_cube.remove_coord('time')
    daily_cube.add_dim_coord(new_time_coord, 0)

    # Remove problematic attributes before saving
    daily_cube.attributes.pop('GRIB_PARAM', None)

    output_filename = f"{outdir}{reanalysis}_{variable}_{y}.nc"
    iris.save(daily_cube, output_filename)

print("Processing complete.")



 