"""
 Script performs MannKendall trend analysis on JRA-55 reanalysis data
 and plots the results on a global map.

To run: python JRA-55_trends_plotter.py <start_year> <end_year>

"""
import numpy as np
import sys
from MK_functions import load_and_mask_yearly_data, create_time_coord, stack_cubes, create_year_cube, compute_kendall_tau, save_results, save_positive_p_values
    
def significance_plot(lat_minimum, lat_maximum):
    datadir = '/data/users/samuel.sancar/Heat_extremes/JRA-55_tailwidths/Yearly_tailwidths/'
    reanalysis = 'jra55'
    variable = '_minmax_surf.015_tmax.reg_tl319.'
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
    datadir = '/data/users/samuel.sancar/Heat_extremes/JRA-55_tailwidths/Yearly_tailwidths/'
    outdir = '/data/users/samuel.sancar/Heat_extremes/Trends/Reanalysis_trends'
    outdir_pos_pvals = '/data/users/samuel.sancar/Heat_extremes/Trends/Reanalysis_trends/'  # <-- For positive trends
    reanalysis = 'jra55'
    variable = '_minmax_surf.015_tmax.reg_tl319.'
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
    save_positive_p_values(tau_matrix, p_matrix, lat_points, lon_points, outdir_pos_pvals, reanalysis, variable, t1, t2)  # <-- NEW

if __name__ == "__main__":
    main()
