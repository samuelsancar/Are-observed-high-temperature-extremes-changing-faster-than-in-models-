"""
Code to calculate MERRA tailwidth trends from 1958 to 2022.
"""

import iris
import iris.util
from iris.analysis import stats
import os
import numpy as np
import sys
from iris.cube import CubeList
import cf_units
from Trend_calculation_functions import (create_time_coord, create_year_cube, maskp87p5, stack_cubes, calculate_trend)

def load_yearly_cubes(datadir, p87p5, reanalysis, variable, years):
    
    cubes = iris.cube.CubeList()
    for year in years:
        filename = f"{reanalysis}_{variable}_{year}_tail_width.nc"
        cube = iris.load_cube(os.path.join(datadir, filename))
        cubes.append(cube)
        
    p87p5_cubes = iris.cube.CubeList()
    for year in years:
        p87p5_filename = f"{variable}_{year}_87p5.nc"
        p87p5_cube = iris.load_cube(os.path.join(p87p5, p87p5_filename))
        p87p5_cubes.append(p87p5_cube)
        
    return cubes, p87p5_cubes

def main():

    t1 = 1980 #int(sys.argv[1])
    t2 = 2023 # int(sys.argv[2])
    years = np.arange(t1, t2)
    
    reanalysis = 'MERRA'
    variable = 'daily_maximum_2m_temperature'
    
    datadir = '/data/users/samuel.sancar/Heat_extremes/MERRA_tailwidths/'
    outdir = '/data/users/samuel.sancar/Heat_extremes/Trends/Reanalysis_trends/'
    p87p5 = '/data/users/samuel.sancar/Heat_extremes/MERRA_tailwidths/'
    
    os.makedirs(outdir, exist_ok=True)

    time_coord = create_time_coord(years)
    cubes, p87p5_cubes = load_yearly_cubes(datadir, p87p5 , reanalysis, variable, years)
    data_cube, p87p5_data_cube = stack_cubes(cubes, p87p5_cubes, time_coord)
    
     # Ensure time coordinate is attached
    if not data_cube.coords('time'):
        data_cube.add_dim_coord(time_coord, 0)
    if not p87p5_data_cube.coords('time'):
        p87p5_data_cube.add_dim_coord(time_coord, 0)
        
    year_cube, p87p5_year_cube = create_year_cube(data_cube, p87p5_data_cube, years)
    trend = calculate_trend(year_cube, data_cube)
    p87p5_trend = calculate_trend(p87p5_year_cube, p87p5_data_cube)
    masked_neg_trends = maskp87p5(p87p5_trend, trend)
    
    masked_neg_trends.rename('Tail width trends (1958-2022)')
    p87p5_trend.rename('p87p5 trends (1958-2022)')

    out_file = os.path.join(outdir, f"{reanalysis}_{variable}_tailwidth_trends_{t1}-{t2}.nc")
    iris.save(masked_neg_trends, out_file)
    
    p87p5_out_file = os.path.join(outdir, f"{reanalysis}_{variable}_p87p5_trends_{t1}-{t2}.nc")
    iris.save(p87p5_trend, p87p5_out_file)

if __name__ == "__main__":
    main()