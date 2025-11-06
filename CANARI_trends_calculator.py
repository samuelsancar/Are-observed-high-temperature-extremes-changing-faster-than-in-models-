"""
03/11/2025
----------
Code to caluclate trends in tail widths and 87.5th percentiles for CANARI model data.
 
"""

import iris
import numpy as np
import sys
from iris.cube import CubeList
from iris.analysis import stats
import iris.quickplot as qplt
import os
import cf_units
from Trend_calculation_function import calculate_trend
from ERA5_trend_calc import maskp87p5

t1 = 1958
t2 = 2023

outtrend = '/data/users/samuel.sancar/Heat_extremes/Trends/CANARI_trends/'
outtailwidth = '/data/users/samuel.sancar/Heat_extremes/Tail_widths/canari/'

for member in range(24,41):
    
    datadir = f'/data/users/samuel.sancar/Heat_extremes/Tail_widths/canari/member{member}/'
    
    os.makedirs(outtrend, exist_ok=True)

    # Slice safely
    def slice_spatial(cube, year):
        coords_available = [coord.name() for coord in cube.coords()]
        print(f"Year {year} coords: {coords_available}")  

        if cube.coords('latitude') and cube.coords('longitude'):
            return cube.slices(['latitude', 'longitude']).next()
        elif cube.coords('grid_latitude') and cube.coords('grid_longitude'):
            return cube.slices(['grid_latitude', 'grid_longitude']).next()
        else:
            print(f"⚠ Skipping year {year}: missing longitude or latitude")
            return None

    # Collect cubes and track missing years
    yearly_tail_widths = iris.cube.CubeList([])
    missing_tailwidth_years = []
    processed_tailwidth_years = []

    for y in range(t1, t2):
        fn = f"tailwidth_member{member}_{y}.nc"
        year_cube = iris.load_cube(datadir + fn)

        sliced = slice_spatial(year_cube, y)
        if sliced is not None:
            yearly_tail_widths.append(sliced)
            processed_tailwidth_years.append(y)
        else:
            missing_tailwidth_years.append(y)

    yearly_p87p5 = iris.cube.CubeList([])
    missing_p87p5_years = []
    processed_p87p5_years = []

    for y in range(t1, t2):
        fn = f"p875_member{member}_{y}.nc"
        p87p5_year_cube = iris.load_cube(datadir + fn)

        sliced = slice_spatial(p87p5_year_cube, y)
        if sliced is not None:
            yearly_p87p5.append(sliced)
            processed_p87p5_years.append(y)
        else:
            missing_p87p5_years.append(y)

    # Report missing years
    if missing_tailwidth_years:
        print(f"\n Missing coords in tailwidth files for years: {missing_tailwidth_years}")
    if missing_p87p5_years:
        print(f"Missing coords in p87p5 files for years: {missing_p87p5_years}\n")

    # Build time coordinates from processed years
    time_unit = cf_units.Unit('days since 1970-01-01', calendar='gregorian')
    time_coord_tailwidth = iris.coords.DimCoord(processed_tailwidth_years, standard_name='time', units=time_unit)
    time_coord_p87p5 = iris.coords.DimCoord(processed_p87p5_years, standard_name='time', units=time_unit)

    # Use spatial coordinates from first valid cube
    lat_coord = yearly_tail_widths[0].coord('latitude') if yearly_tail_widths[0].coords('latitude') else yearly_tail_widths[0].coord('grid_latitude')
    lon_coord = yearly_tail_widths[0].coord('longitude') if yearly_tail_widths[0].coords('longitude') else yearly_tail_widths[0].coord('grid_longitude')

    p87p5_lat_coord = yearly_p87p5[0].coord('latitude') if yearly_p87p5[0].coords('latitude') else yearly_p87p5[0].coord('grid_latitude')
    p87p5_lon_coord = yearly_p87p5[0].coord('longitude') if yearly_p87p5[0].coords('longitude') else yearly_p87p5[0].coord('grid_longitude')

    # Equalise attributes and stack data
    iris.util.equalise_attributes(yearly_tail_widths)
    data_stack = np.stack([cube.data for cube in yearly_tail_widths], axis=0)
    iris.util.equalise_attributes(yearly_p87p5)
    p87p5_data_stack = np.stack([cube.data for cube in yearly_p87p5], axis=0)

    # Create merged cubes
    yearly_tail_widths = iris.cube.Cube(
        data_stack,
        dim_coords_and_dims=[(time_coord_tailwidth, 0), (lat_coord, 1), (lon_coord, 2)],
        long_name='tail_width_over_years',
        units=yearly_tail_widths[0].units,
        attributes=yearly_tail_widths[0].attributes
    )

    yearly_p87p5 = iris.cube.Cube(
        p87p5_data_stack,
        dim_coords_and_dims=[(time_coord_p87p5, 0), (p87p5_lat_coord, 1), (p87p5_lon_coord, 2)],
        long_name='p87p5_over_years',
        units=yearly_p87p5[0].units,
        attributes=yearly_p87p5[0].attributes
    )

    # Prepare yearly cubes for trend calculation 
    yearly = yearly_tail_widths.copy()
    year_cube_data = np.repeat(np.array(processed_tailwidth_years)[:, None, None],
                                yearly.shape[1] * yearly.shape[2]).reshape(yearly.shape)
    yearly.data = year_cube_data

    p87p5_yearly = yearly_p87p5.copy()
    p87p5_cube_data = np.repeat(np.array(processed_p87p5_years)[:, None, None],
                                p87p5_yearly.shape[1] * p87p5_yearly.shape[2]).reshape(p87p5_yearly.shape)
    p87p5_yearly.data = p87p5_cube_data

    # Save tailwidth cube
    out_file_tailwidth = os.path.join(outtailwidth, f"{member}_tailwidths.nc")
    iris.save(yearly_tail_widths, out_file_tailwidth)

    # Calculate trend
    trend = calculate_trend(yearly, yearly_tail_widths)
    p87p5_trend = calculate_trend(p87p5_yearly, yearly_p87p5)
    masked_trend = maskp87p5(p87p5_trend, trend)

    # Long name
    masked_trend.long_name = 'Tailwidth trends (°C/decade)'
    p87p5_trend.long_name = 'p87p5 trends (°C/decade)'  

    out_file_trend = os.path.join(outtrend, f"{member}_tailwidth_trend.nc")
    iris.save(masked_trend, out_file_trend)

    p87p5_out_file_trend = os.path.join(outtrend, f"{member}_p87p5_trend.nc")
    iris.save(p87p5_trend, p87p5_out_file_trend)