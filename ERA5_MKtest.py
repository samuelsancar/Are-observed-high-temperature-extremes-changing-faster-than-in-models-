"""

24/07/25
---------
Code that applies MannKendall test on ERA5 tailwidth trends. Output is KendallTau and p-values (for all trends and positive only trends).

To run: python ERA5_MK_tests.py [1958] [2023]

"""

import sys
import numpy as np
import xarray as xr
import os
from MK_functions import (load_and_mask_yearly_data,create_time_coord,stack_cubes,create_year_cube,compute_kendall_tau, save_results, save_positive_p_values)

def significance_plot(lat_minimum, lat_maximum):
    datadir = '/data/users/samuel.sancar/Heat_extremes/ERA5_tailwidths/'
    reanalysis = 'era5'
    variable = 'daily_maximum_2m_temperature'
    t1 = 1958
    t2 = 2023
    lat_min = lat_minimum
    lat_max = lat_maximum
    years = np.arange(t1, t2)
    years = np.arange(t1, t2)

    time_coord = create_time_coord(years)
    cubes = load_and_mask_yearly_data(datadir, reanalysis, variable, years, lat_min, lat_max)
    data_cube = stack_cubes(cubes, time_coord)
    year_cube = create_year_cube(data_cube, years)
    tau_matrix, p_matrix, lat_points, lon_points = compute_kendall_tau(year_cube, data_cube)
    return p_matrix, lat_points, lon_points

def main():
    datadir = '/data/users/samuel.sancar/Heat_extremes/ERA5_tailwidths/'
    outdir = '/data/users/samuel.sancar/Heat_extremes/Trends/Reanalysis_trends'
    outdir_pos = outdir
    reanalysis = 'era5'
    variable = 'daily_maximum_2m_temperature'
    t1 = 1958
    t2 = 2023
    lat_min = -60
    lat_max = int(sys.argv[1])
    years = np.arange(t1, t2)

    time_coord = create_time_coord(years)
    cubes = load_and_mask_yearly_data(datadir, reanalysis, variable, years, lat_min, lat_max)
    data_cube = stack_cubes(cubes, time_coord)
    year_cube = create_year_cube(data_cube, years)
    tau_matrix, p_matrix, lat_points, lon_points = compute_kendall_tau(year_cube, data_cube)

    save_results(tau_matrix, p_matrix, lat_points, lon_points, outdir, reanalysis, variable, t1, t2)
    save_positive_p_values(tau_matrix, p_matrix, lat_points, lon_points, outdir_pos, reanalysis, variable, t1, t2)

if __name__ == "__main__":
    main()
