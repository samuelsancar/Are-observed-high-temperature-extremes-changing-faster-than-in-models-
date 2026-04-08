"""
 Script performs MannKendall trend analysis on JRA-55 reanalysis data
 and plots the results on a global map.

To run: python JRA-55_trends_plotter.py <start_year> <end_year>

"""
import numpy as np
import sys
import iris
from MK_functions import load_and_mask_yearly_data, create_time_coord, stack_cubes, create_year_cube, compute_kendall_tau, save_results, save_positive_p_values

def main():
    
    t1 = 1958
    t2 = 2023
    lat_min = -60 
    lat_max = 90
    years = np.arange(t1, t2)
    
    reanalysis = 'jra55'
    variable = '_minmax_surf.015_tmax.reg_tl319.'
    datadir = '/data/users/samuel.sancar/Heat_extremes/JRA-55_tailwidths/'
    #correlation = iris.load_cube('/data/users/samuel.sancar/Heat_extremes/Trends/Reanalysis_trends/ERA5_JRA-55_correlation.nc')
    outdir = '/data/users/samuel.sancar/Heat_extremes/Trends/Reanalysis_trends'
    p87p5 = '/data/users/samuel.sancar/Heat_extremes/Trends/Reanalysis_trends/'
    p87p5_cube = iris.load_cube(f"{p87p5}jra55__minmax_surf.015_tmax.reg_tl319._p87p5_trends_{t1}-2023.nc")

    time_coord = create_time_coord(years)
    cubes = load_and_mask_yearly_data(datadir, reanalysis, variable, years, lat_min, lat_max)
    data_cube = stack_cubes(cubes, time_coord)
    year_cube = create_year_cube(data_cube, years)
    
    # Mask regions with negative trends in p87p5 percentile.
    mask = (p87p5_cube.data < 0) #| (correlation.data < 0.5)
    masked_data = np.where(mask, np.nan, data_cube.data)
    data_cube_masked = data_cube.copy()
    data_cube_masked.data = masked_data
    
    tau_matrix, p_matrix, lat_points, lon_points = compute_kendall_tau(year_cube, data_cube_masked)
    save_results(tau_matrix, p_matrix, lat_points, lon_points, outdir, reanalysis, variable, t1, t2)
    save_positive_p_values(tau_matrix, p_matrix, lat_points, lon_points, outdir, reanalysis, variable, t1, t2)

if __name__ == "__main__":
    main()
   