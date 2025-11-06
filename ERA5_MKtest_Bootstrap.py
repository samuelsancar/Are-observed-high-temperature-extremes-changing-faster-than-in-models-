"""

24/07/25
---------

Code that applies MannKendall test and Bootstrapping on ERA5 tailwidth trends. Output is KendallTau and p-values (for all trends and positive only trends).

To run: python ERA5_MK_tests.py [1958] [2023]

"""

import sys
import numpy as np
import xarray as xr
import iris
import os
from MK_functions import (load_and_mask_yearly_data,create_time_coord,stack_cubes,create_year_cube,compute_kendall_tau, save_results, save_positive_p_values)
from Bootstrap_function import bootstrap_tailwidth_trends, save_results_b_mask

def significance_plot(lat_minimum, lat_maximum):      # Function is called when using ERA5_trend_plotter.py
    datadir = '/data/users/samuel.sancar/Heat_extremes/ERA5_tailwidths/'
    p87p5 = '/data/users/samuel.sancar/Heat_extremes/Trends/Reanalysis_trends/'
    p87p5_cube = iris.load_cube(f"{p87p5}era5_daily_maximum_2m_temperature_p87p5_trends_1958-2023.nc")
    
    reanalysis = 'era5'
    variable = 'daily_maximum_2m_temperature'
    t1 = 1958
    t2 = 2023
    lat_min = lat_minimum
    lat_max = lat_maximum
    iterations = 1000
    years = np.arange(t1, t2)

    time_coord = create_time_coord(years)
    cubes = load_and_mask_yearly_data(datadir, reanalysis, variable, years, lat_min, lat_max)
    data_cube = stack_cubes(cubes, time_coord)
    year_cube = create_year_cube(data_cube, years)

    masked_data = np.where(p87p5_cube.data < 0, np.nan, data_cube.data)
    data_cube_masked = data_cube.copy()
    data_cube_masked.data = masked_data

    tau_matrix, p_matrix, lat_points, lon_points = compute_kendall_tau(year_cube, data_cube_masked)

    return p_matrix, lat_points, lon_points

def main():               # 'Main' example when applied to ERA5 tailwidth data. (Standalone script to run MK and Bootstrap)
    datadir = '/data/users/samuel.sancar/Heat_extremes/ERA5_tailwidths/'
    outdir = '/data/users/samuel.sancar/Heat_extremes/Trends/Reanalysis_trends'
    outdir_pos = outdir
    reanalysis = 'era5'
    variable = 'daily_maximum_2m_temperature'
    t1 = 1958
    t2 = 2023
    lat_min = -60
    lat_max = int(sys.argv[1])
    iterations = int(sys.argv[2])
    ci = float(sys.argv[3])
    years = np.arange(t1, t2)

    time_coord = create_time_coord(years)
    cubes = load_and_mask_yearly_data(datadir, reanalysis, variable, years, lat_min, lat_max)
    data_cube = stack_cubes(cubes, time_coord)
    year_cube = create_year_cube(data_cube, years)
    tau_matrix, p_matrix, lat_points, lon_points = compute_kendall_tau(year_cube, data_cube)
    b_matrix, _, _ = bootstrap_tailwidth_trends(data_cube, iterations, ci)

    save_results(tau_matrix, p_matrix, lat_points, lon_points, outdir, reanalysis, variable, t1, t2)
    save_positive_p_values(tau_matrix, p_matrix, lat_points, lon_points, outdir_pos, reanalysis, variable, t1, t2)
    save_results_b_mask(b_matrix)


if __name__ == "__main__":
    main()
