"""
 Script performs MannKendall trend analysis on JRA-55 reanalysis data
 and plots the results on a global map.

    To run: python JRA-55_trends_plotter.py <start_year> <end_year>

"""
import numpy as np
import sys
import iris
from iris.analysis import stats
import iris.util
from iris.cube import CubeList
import cf_units
import os
from MK_functions import create_time_coord, stack_cubes, create_year_cube, compute_kendall_tau, save_results
from Lat_Land_mask_function import lat_land_mask

def load_and_mask_yearly_data(datadir, reanalysis, variable, years, lat_min=-60, lat_max=60):
    cubes = iris.cube.CubeList()
    for year in years:
        filename = f"{reanalysis}{variable}_{year}_tailwidth.nc"
        cube = iris.load_cube(os.path.join(datadir, filename))
        masked_cube, _ = lat_land_mask(cube, lat_min, lat_max)
        cubes.append(masked_cube)
    return cubes

def main():
    datadir = '/data/users/samuel.sancar/Heat_extremes/JRA-55_tailwidths/Yearly_tailwidths/'
    outdir = '/data/users/samuel.sancar/Heat_extremes/Trends/Reanalysis_trends'
    reanalysis = 'jra55'
    variable = '_minmax_surf.015_tmax.reg_tl319.'
    t1 = int(sys.argv[1])
    t2 = int(sys.argv[2])
    years = np.arange(t1, t2)

    time_coord = create_time_coord(years)
    cubes = load_and_mask_yearly_data(datadir, reanalysis, variable, years)
    data_cube = stack_cubes(cubes, time_coord)
    year_cube = create_year_cube(data_cube, years)
    tau_matrix, p_matrix, lat_points, lon_points = compute_kendall_tau(year_cube, data_cube)
    save_results(tau_matrix, p_matrix, lat_points, lon_points, outdir, reanalysis, variable, t1, t2)

if __name__ == "__main__":
    main()