"""
    Functions to prepare iris cubes and calculate 
    linear trends using pearson r coefficient.

"""

import os
import numpy as np
import iris
from iris.analysis import stats
import iris.util
from iris.cube import CubeList
from iris.coords import AuxCoord
import cf_units
from scipy.stats import pearsonr

def add_year_coord(cube, first_year: int):
    '''Add a year coordinate to the cube based on the time dimension.'''
    
    n = cube.shape[0]
    sim_years = np.arange(first_year, first_year + n, dtype=int)

    # Attach as an AuxCoord along the time dimension
    time_dim = cube.coord_dims('time')[0]
    sim_year_aux = AuxCoord(sim_years, long_name='sim_year', units='year')
    cube.add_aux_coord(sim_year_aux, data_dims=(time_dim,))
    return cube

def calculate_trend(year_cube, data_cube):
    '''Calculate linear trend using Pearson correlation coefficient at each grid point.'''
    
    r = iris.analysis.stats.pearsonr(year_cube, data_cube, corr_coords='time')
    std_data = data_cube.collapsed('time', iris.analysis.STD_DEV)
    std_year = year_cube.collapsed('time', iris.analysis.STD_DEV)
    return r * std_data / std_year * 10   # decadal trend

def regional_trends(regional_series):
    '''Calculate linear trend using Pearson correlation coefficient for a regional series.'''
    
    years = regional_series.coord('year').points
    values = regional_series.data
    
    # Pearson correlation
    r, p_value = pearsonr(years, values)

    return r * (values.std() / years.std()) * 10  # decadal trend

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
    '''Create a time coordinate with units of days since 1970-01-01.'''
    
    time_unit = cf_units.Unit('days since 1970-01-01', calendar='gregorian')
    
    return iris.coords.DimCoord(years, standard_name='time', long_name='valid_time', units=time_unit)

def stack_cubes(cubes, p87p5_cubes, time_coord):
    '''Stack yearly cubes into a single cube with time as the first dimension.'''
    
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
    '''Mask tail width values where p87.5 is negative.'''
    
    masked_data = np.where(p87p5_cube.data < 0, np.nan, data_cube.data)
    masked_cube = data_cube.copy()
    masked_cube.data = masked_data
    return masked_cube

def create_year_cube(data_cube, p87p5_data_cube, years):
    '''Create cubes with the same shape as data_cube but with year values for each time step.'''
    
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

