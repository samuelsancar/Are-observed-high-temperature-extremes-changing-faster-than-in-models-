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
from iris.coords import AuxCoord
import os
import cf_units
from Trend_calculation_functions import (create_time_coord, create_year_cube, maskp87p5, stack_cubes, calculate_trend)
from ERA5_trend_calc import maskp87p5

t1 = 1980
t2 = 2023
years = np.arange(t1, t2)

outtrend = '/data/users/samuel.sancar/Heat_extremes/Trends/CANARI_trends/'
outtailwidth = '/data/users/samuel.sancar/Heat_extremes/Tail_widths/canari/'

def add_year_coord(cube, first_year: int):
    n = cube.shape[0]
    sim_years = np.arange(first_year, first_year + n, dtype=int)

    # Attach as an AuxCoord along the time dimension
    time_dim = cube.coord_dims('time')[0]
    sim_year_aux = AuxCoord(sim_years, long_name='sim_year', units='year')
    cube.add_aux_coord(sim_year_aux, data_dims=(time_dim,))
    return cube
    
def load_yearly_cubes(datadir, year, member):
        
    cubes = iris.cube.CubeList()
    p87p5_cubes = iris.cube.CubeList()
    
    for year in years:
        
        fn = f"tailwidth_member{member}_{year}.nc"
        fn_87p5 = f"p875_member{member}_{year}.nc"
            
        year_cube = iris.load_cube(datadir + fn)
        p87p5_year_cube = iris.load_cube(datadir + fn_87p5)
        
        cubes.append(year_cube)
        p87p5_cubes.append(p87p5_year_cube)
            
    return cubes, p87p5_cubes

for member in range(1, 41):

    datadir = f'/data/users/samuel.sancar/Heat_extremes/Tail_widths/canari/member{member}/'
        
    # Stack yealy tailwidths into sinle file.
    cubes, p87p5_cubes = load_yearly_cubes(datadir, years, member)
    time_coord = create_time_coord(years)
    data_cube, p87p5_cube = stack_cubes(cubes, p87p5_cubes, time_coord)
    
    data_cube = add_year_coord(data_cube, t1)
    p87p5_cube = add_year_coord(p87p5_cube, t1)
    
    data_cube.long_name = f'Yealy tail widths ({t1}-2022)'
    p87p5_cube.long_name = f'Yearly p87p5 ({t1}-2022)'

    # Save stacked tail tailwidth cube
    out_file_tailwidth = os.path.join(outtailwidth, f"{member}_tailwidths_{t1}-{t2}.nc")
    iris.save(data_cube, out_file_tailwidth)

    # Calculate trend
    year_cube, p87p5_year_cube = create_year_cube(data_cube, p87p5_cube, years)
    trend = calculate_trend(year_cube, data_cube)
    p87p5_trend = calculate_trend(p87p5_year_cube, p87p5_cube)
    masked_trend = maskp87p5(p87p5_trend, trend)

    # Long name
    masked_trend.rename(f'Tail width trends ({t1}-{t2})')
    p87p5_trend.rename(f'p87p5 trends ({t1}-{t2})') 

    out_file_trend = os.path.join(outtrend, f"{member}_tailwidth_trend_{t1}-{t2}.nc")
    iris.save(masked_trend, out_file_trend)

    p87p5_out_file_trend = os.path.join(outtrend, f"{member}_p87p5_trend_{t1}-{t2}.nc")
    iris.save(p87p5_trend, p87p5_out_file_trend)