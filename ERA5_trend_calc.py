"""
17/07/25
---------
Code to calculate trends in ERA5 daily tasmax tailwidths.

To run: python ERA5_trend_calc [1958] [2023]

"""

import os
import sys
import numpy as np
import iris
from iris.analysis import stats
import iris.util
from iris.cube import CubeList
import cf_units

def load_yearly_cubes(datadir, p87p5, reanalysis, variable, years):
    
    cubes = iris.cube.CubeList()
    for year in years:
        filename = f"{reanalysis}_{variable}_{year}_tail_width.nc"
        cube = iris.load_cube(os.path.join(datadir, filename))
        cubes.append(cube)
        
    p87p5_cubes = iris.cube.CubeList()
    for year in years:
        p87p5_filename = f"{reanalysis}_{variable}_{year}_87p5.nc"
        p87p5_cube = iris.load_cube(os.path.join(p87p5, p87p5_filename))
        p87p5_cubes.append(p87p5_cube)
        
    return cubes, p87p5_cubes
    

def create_time_coord(years):
    time_unit = cf_units.Unit('days since 1970-01-01', calendar='gregorian')
    return iris.coords.DimCoord(years, standard_name='time', long_name='valid_time', units=time_unit)

def stack_cubes(cubes, p87p5_cubes, time_coord):
    iris.util.equalise_attributes(cubes)
    data_stack = np.stack([cube.data for cube in cubes], axis=0)
    lat_coord = cubes[0].coord('latitude')
    lon_coord = cubes[0].coord('longitude')
    
    iris.util.equalise_attributes(p87p5_cubes)
    p87p5_data_stack = np.stack([p87p5_cube.data for p87p5_cube in p87p5_cubes], axis=0)
    p87p5_lat_coord = p87p5_cubes[0].coord('latitude')
    p87p5_lon_coord = p87p5_cubes[0].coord('longitude')
    
    return iris.cube.Cube(
        data_stack,
        dim_coords_and_dims=[(time_coord, 0), (lat_coord, 1), (lon_coord, 2)],
        long_name='tail_width_over_years',
        units=cubes[0].units,
        attributes=cubes[0].attributes
    ), iris.cube.Cube(
        p87p5_data_stack,
        dim_coords_and_dims=[(time_coord, 0), (p87p5_lat_coord, 1), (p87p5_lon_coord, 2)],
        long_name='tail_width_over_years',
        units=p87p5_cubes[0].units,
        attributes=p87p5_cubes[0].attributes)
    

def maskp87p5(p87p5_cube, data_cube):
    masked_data = np.where(p87p5_cube.data < 0, np.nan, data_cube.data)
    masked_cube = data_cube.copy()
    masked_cube.data = masked_data
    return masked_cube

def create_year_cube(data_cube, p87p5_data_cube, years):
    year_data = np.empty_like(data_cube.data)
    for i, year in enumerate(years):
        year_data[i, :, :] = year
    year_cube = data_cube.copy()
    year_cube.data = year_data

    p87p5_year_data = np.empty_like(p87p5_data_cube.data)
    for i, year in enumerate(years):
        p87p5_year_data[i, :, :] = year
    p87p5_year_cube = p87p5_data_cube.copy()
    p87p5_year_cube.data = p87p5_year_data  

    return year_cube, p87p5_year_cube

def calculate_trend(year_cube, data_cube, p87p5_year_cube, p87p5_data_cube):
    r = iris.analysis.stats.pearsonr(year_cube, data_cube, corr_coords='time')
    std_data = data_cube.collapsed('time', iris.analysis.STD_DEV)
    std_year = year_cube.collapsed('time', iris.analysis.STD_DEV)

    p87p5_r = iris.analysis.stats.pearsonr(p87p5_year_cube, p87p5_data_cube, corr_coords='time')
    p87p5_std_data = p87p5_data_cube.collapsed('time', iris.analysis.STD_DEV)
    p87p5_std_year = p87p5_year_cube.collapsed('time', iris.analysis.STD_DEV)
    return r * std_data / std_year * 10, p87p5_r * p87p5_std_data / p87p5_std_year * 10  # decadal trend

def main():
    datadir = '/data/users/samuel.sancar/Heat_extremes/ERA5_tailwidths/'
    outdir = '/data/users/samuel.sancar/Heat_extremes/Trends/Reanalysis_trends/'
    p87p5 = '/data/users/samuel.sancar/Heat_extremes/ERA5_tailwidths/'
    os.makedirs(outdir, exist_ok=True)

    reanalysis = 'era5'
    variable = 'daily_maximum_2m_temperature'
    t1 = int(sys.argv[1])
    t2 = int(sys.argv[2])
    years = np.arange(t1, t2)

    time_coord = create_time_coord(years)
    cubes, p87p5_cubes = load_yearly_cubes(datadir, p87p5 , reanalysis, variable, years)
    data_cube, p87p5_data_cube = stack_cubes(cubes, p87p5_cubes, time_coord)
    
     # Ensure time coordinate is attached
    if not data_cube.coords('time'):
        data_cube.add_dim_coord(time_coord, 0)
    if not p87p5_data_cube.coords('time'):
        p87p5_data_cube.add_dim_coord(time_coord, 0)
        
    year_cube, p87p5_year_cube = create_year_cube(data_cube, p87p5_data_cube, years)
    trend, p87p5_trend = calculate_trend(year_cube, data_cube, p87p5_year_cube, p87p5_data_cube)
    masked_neg_trends = maskp87p5(p87p5_trend, trend)

    out_file = os.path.join(outdir, f"{reanalysis}_{variable}_tailwidth_trends_{t1}-{t2}.nc")
    iris.save(masked_neg_trends, out_file)
    
    p87p5_out_file = os.path.join(outdir, f"{reanalysis}_{variable}_p87p5_trends_{t1}-{t2}.nc")
    iris.save(p87p5_trend, p87p5_out_file)

if __name__ == "__main__":
    main()


    
    



