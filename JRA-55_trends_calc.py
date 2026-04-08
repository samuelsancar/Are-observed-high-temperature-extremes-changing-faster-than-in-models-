"""
    Code to calulate trends in tail widths of JRA-55 reanalysis data (p99 - p87.5 of daily maximum 2m temperature)
    over 1958-2023.

"""

import os
import sys
import numpy as np
import iris
from iris.analysis import stats
import iris.util
from iris.cube import CubeList
from iris.coords import AuxCoord
import cf_units
from Trend_calculation_functions import (create_time_coord, create_year_cube, maskp87p5, load_yearly_cubes, stack_cubes, calculate_trend)

def add_year_coord(cube, first_year: int):
    n = cube.shape[0]
    sim_years = np.arange(first_year, first_year + n, dtype=int)

    # Attach as an AuxCoord along the time dimension
    time_dim = cube.coord_dims('time')[0]
    sim_year_aux = AuxCoord(sim_years, long_name='sim_year', units='year')
    cube.add_aux_coord(sim_year_aux, data_dims=(time_dim,))
    return cube

def main():
    datadir = '/data/users/samuel.sancar/Heat_extremes/JRA-55_tailwidths/'
    outdir = '/data/users/samuel.sancar/Heat_extremes/Trends/Reanalysis_trends/'
    p87p5 = '/data/users/samuel.sancar/Heat_extremes/JRA-55_tailwidths/'
    os.makedirs(outdir, exist_ok=True)

    reanalysis = 'jra55_'
    variable = 'minmax_surf.015_tmax.reg_tl319.'
    t1 = 1958
    t2 = 2023
    years = np.arange(t1, t2)

    time_coord = create_time_coord(years)
    cubes, p87p5_cubes = load_yearly_cubes(datadir, p87p5 , reanalysis, variable, years)
    data_cube, p87p5_data_cube = stack_cubes(cubes, p87p5_cubes, time_coord)
     
     # Ensure time coordinate is attached
    if not data_cube.coords('time'):
        data_cube.add_dim_coord(time_coord, 0)
    if not p87p5_data_cube.coords('time'):
        p87p5_data_cube.add_dim_coord(time_coord, 0)
        
    data_cube = add_year_coord(data_cube, t1)
    p87p5_data_cube = add_year_coord(p87p5_data_cube, t1)
    
    data_cube.long_name = f'Yealy tail widths ({t1}-2022)'
    p87p5_data_cube.long_name = f'Yearly p87p5 ({t1}-2022)'
        
    year_cube, p87p5_year_cube = create_year_cube(data_cube, p87p5_data_cube, years)
    trend = calculate_trend(year_cube, data_cube)
    p87p5_trend = calculate_trend(p87p5_year_cube, p87p5_data_cube)
    masked_neg_trends = maskp87p5(p87p5_trend, trend)
    
    masked_neg_trends.rename(f'Tail width trends ({t1}-2022)')
    p87p5_trend.rename(f'p87p5 trends ({t1}-2022)')

    out_file = os.path.join(outdir, f"{reanalysis}_{variable}_tailwidth_trends_{t1}-{t2}.nc")
    iris.save(masked_neg_trends, out_file)
    
    p87p5_out_file = os.path.join(outdir, f"{reanalysis}_{variable}_p87p5_trends_{t1}-{t2}.nc")
    iris.save(p87p5_trend, p87p5_out_file)
    
if __name__ == "__main__":
    main()