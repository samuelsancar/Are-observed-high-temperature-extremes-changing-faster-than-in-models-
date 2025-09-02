"""
17/07/25
---------
Code to calculate trends in ERA5 daily tasmax tailwidths.

To run: python ERA5_trends_calc [1958] [2023]

"""

import os
import sys
import numpy as np
import iris
import iris.analysis
import iris.util
import cf_units

def load_yearly_cubes(datadir, reanalysis, variable, years):
    cubes = iris.cube.CubeList()
    for year in years:
        filename = f"{reanalysis}_{variable}_{year}_tail_width.nc"
        cube = iris.load_cube(os.path.join(datadir, filename))
        cubes.append(cube)
    return cubes

def create_time_coord(years):
    time_unit = cf_units.Unit('days since 1970-01-01', calendar='gregorian')
    return iris.coords.DimCoord(years, standard_name='time', long_name='valid_time', units=time_unit)

def stack_cubes(cubes, time_coord):
    iris.util.equalise_attributes(cubes)
    data_stack = np.stack([cube.data for cube in cubes], axis=0)
    lat_coord = cubes[0].coord('latitude')
    lon_coord = cubes[0].coord('longitude')
    return iris.cube.Cube(
        data_stack,
        dim_coords_and_dims=[(time_coord, 0), (lat_coord, 1), (lon_coord, 2)],
        long_name='tail_width_over_years',
        units=cubes[0].units,
        attributes=cubes[0].attributes
    )

def create_year_cube(data_cube, years):
    year_data = np.empty_like(data_cube.data)
    for i, year in enumerate(years):
        year_data[i, :, :] = year
    year_cube = data_cube.copy()
    year_cube.data = year_data
    return year_cube

def calculate_trend(year_cube, data_cube):
    r = iris.analysis.stats.pearsonr(year_cube, data_cube, corr_coords='time')
    std_data = data_cube.collapsed('time', iris.analysis.STD_DEV)
    std_year = year_cube.collapsed('time', iris.analysis.STD_DEV)
    return r * std_data / std_year * 10  # decadal trend

def main():
    datadir = '/data/users/samuel.sancar/Heat_extremes/tailwidthslist/era5/'
    outdir = '/data/users/samuel.sancar/Heat_extremes/trends/Reanalysis_trends/'
    os.makedirs(outdir, exist_ok=True)

    reanalysis = 'era5'
    variable = 'daily_maximum_2m_temperature'
    t1 = int(sys.argv[1])
    t2 = int(sys.argv[2])
    years = np.arange(t1, t2)

    time_coord = create_time_coord(years)
    cubes = load_yearly_cubes(datadir, reanalysis, variable, years)
    data_cube = stack_cubes(cubes, time_coord)
    year_cube = create_year_cube(data_cube, years)
    trend = calculate_trend(year_cube, data_cube)

    out_file = os.path.join(outdir, f"{reanalysis}_{variable}_tailwidth_trends_{t1}-{t2}.nc")
    iris.save(trend, out_file)

if __name__ == "__main__":
    main()


