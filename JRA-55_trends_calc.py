"""
    Code to calulate trends in tail widths of JRA-55 reanalysis data
    (p99 - p87.5 of daily maximum 2m temperature)
    over 1958-2023.

    To run: python3 JRA-55_trends.py

"""

import os
import sys
import numpy as np
import iris
from iris.analysis import stats
import iris.util
from iris.cube import CubeList
import cf_units
from ERA5_trend_calc import load_yearly_cubes, stack_cubes, create_year_cube, calculate_trend, maskp87p5
import matplotlib.pyplot as plt

def create_time_coord(years):
    time_unit = cf_units.Unit('days since 1970-01-01', calendar='gregorian')
    return iris.coords.DimCoord(years, standard_name='time', long_name='time', units=time_unit)

def main():
    datadir = '/data/users/samuel.sancar/Heat_extremes/JRA-55_tailwidths/Yearly_tailwidths/'
    outdir = '/data/users/samuel.sancar/Heat_extremes/Trends/Reanalysis_trends/'
    p87p5 = '/data/users/samuel.sancar/Heat_extremes/JRA-55_tailwidths/Yearly_tailwidths/'
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
        
    year_cube, p87p5_year_cube = create_year_cube(data_cube, p87p5_data_cube, years)
    trend, p87p5_trend = calculate_trend(year_cube, data_cube, p87p5_year_cube, p87p5_data_cube)
    masked_neg_trends = maskp87p5(p87p5_trend, trend)

    out_file = os.path.join(outdir, f"{reanalysis}_{variable}_tailwidth_trends_{t1}-{t2}.nc")
    iris.save(masked_neg_trends, out_file)
    
    p87p5_out_file = os.path.join(outdir, f"{reanalysis}_{variable}_p87p5_trends_{t1}-{t2}.nc")
    iris.save(p87p5_trend, p87p5_out_file)
    
if __name__ == "__main__":
    main()