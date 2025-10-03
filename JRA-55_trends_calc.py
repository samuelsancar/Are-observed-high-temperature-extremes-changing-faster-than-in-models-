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
from ERA5_trend_calc import load_yearly_cubes, stack_cubes, create_year_cube, calculate_trend
import matplotlib.pyplot as plt


def load_yearly_cubes(datadir, reanalysis, variable, years):
    cubes = iris.cube.CubeList()
    for year in years:
        filename = f"{reanalysis}{variable}_{year}_tail_width.nc"
        cube = iris.load_cube(os.path.join(datadir, filename))
        cubes.append(cube)
    return cubes

def create_time_coord(years):
    time_unit = cf_units.Unit('days since 1970-01-01', calendar='gregorian')
    return iris.coords.DimCoord(years, standard_name='time', long_name='time', units=time_unit)

def main():
    datadir = '/data/users/samuel.sancar/Heat_extremes/JRA-55_tailwidths/Yearly_tailwidths/'
    outdir = '/data/users/samuel.sancar/Heat_extremes/Trends/Reanalysis_trends/'
    os.makedirs(outdir, exist_ok=True)

    reanalysis = 'jra55_'
    variable = '_minmax_surf.015_tmax.reg_tl319.'
    t1 = 1958
    t2 = 2023
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