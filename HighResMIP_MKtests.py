"""
    Code to calculate KendallTau and p-values in model trends:

    To run: python Model_sig_trends
"""

import iris
import glob
import numpy as np
import os
import sys
import datetime as dt
import xarray as xr
import matplotlib.pyplot as plt
import iris.coord_categorisation as icc
from iris.coords import AuxCoord
from cftime import DatetimeGregorian as CFDT
from Lat_Land_mask_function import ThreeD_lat_land_mask
from MK_functions import compute_kendall_tau_models

lat_minimum = -60
lat_maximum = 90
t1 = sys.argv[1]
years = np.arange(int(t1), 2023)

template_cube = iris.load_cube('template_cube_1deg.nc')
Model_dir = '/data/users/samuel.sancar/Heat_extremes/Tail_widths/*_tailwidths_1958-2022.nc'
correlation = iris.load_cube('/data/users/samuel.sancar/Heat_extremes/Trends/Reanalysis_trends/ERA5_JRA-55_correlation.nc')
outdir = '/data/users/samuel.sancar/Heat_extremes/Trends/Trend_significance/P-vals/'
outdir_pos_pvals = '/data/users/samuel.sancar/Heat_extremes/Trends/Trend_significance/Positive_Pvals/'  # Positive trend p-vals

# Get filenames
Model_filenames = sorted(glob.glob(Model_dir, recursive=True))
print(f"Found {len(Model_filenames)} files to process.")

def add_year_coord(cube, first_year: int):
    n = cube.shape[0]
    sim_years = np.arange(first_year, first_year + n, dtype=int)

    # Attach as an AuxCoord along the time dimension
    time_dim = cube.coord_dims('time')[0]
    sim_year_aux = AuxCoord(sim_years, long_name='sim_year', units='year')
    cube.add_aux_coord(sim_year_aux, data_dims=(time_dim,))
    return cube

def save_results(tau_matrix, p_matrix, lat_points, lon_points, filename):

    base_name = os.path.splitext(os.path.basename(filename))[0]
    tau_save = xr.DataArray(tau_matrix, coords={"latitude": lat_points, "longitude": lon_points}, dims=["latitude", "longitude"], name = "Tailwidth_mktau")
    p_save = xr.DataArray(p_matrix, coords={"latitude": lat_points, "longitude": lon_points}, dims=["latitude", "longitude"], name = "Tailwidth_trend_p-values")
    tau_save.to_netcdf(os.path.join(outdir, f"{base_name}_mktau_{t1}-2022.nc"))
    p_save.to_netcdf(os.path.join(outdir, f"{base_name}_p_{t1}.nc"))
    
    print(f"Saved results for: {filename}")

def save_positive_p_values(tau_matrix, p_matrix, lat_points, lon_points, filename):
    
    pos_p_matrix = np.where(tau_matrix > 0, p_matrix, np.nan)
    pos_p_save = xr.DataArray(pos_p_matrix, coords={"latitude": lat_points, "longitude": lon_points}, dims=["latitude", "longitude"], name = 'Positive_trend_p-values')
    base_name = os.path.splitext(os.path.basename(filename))[0]
    pos_p_save.to_netcdf(os.path.join(outdir_pos_pvals, f"{base_name}_pos_p_{t1}.nc"))
    
    print(f"Saved positive trend p-values for: {filename}")

def load_valid_cube(filepath):
 
    cubes = iris.load(filepath)

    # Determine model
    Mode = 'SST' if filepath.split('/')[-1].split('_')[1] == 'highresSST-present' else 'Coupled'
    Model = filepath.split('/')[-1].split('_')[0]
    experiment = filepath.split('/')[-1].split('_')[1]
    member = filepath.split('/')[-1].split('_')[2]

    p87p5_dir = f'/data/users/samuel.sancar/Heat_extremes/Trends/{Mode}_model_trends/'
    p87p5_cube = iris.load_cube(p87p5_dir + f'{Model}_{experiment}_{member}_p87p5_trend_{t1}-2022.nc')

    # Require time dimension and year coordinate
    for cube in cubes:
        if cube.ndim == 3 and cube.coord('time'):  
                 
            cube.remove_coord('sim_year')
            first_year = 1958
            cube = add_year_coord(cube, first_year)
            
            if t1=='1980':
            
                print('Constraining')
                year_con = iris.Constraint(sim_year=lambda y: 1980 <= y <= 2022)
                cube = cube.extract(year_con)
                
                print(cube)
            
            else:
                continue
           
            return cube, p87p5_cube

def process_model_files(filenames):
    
    processed = 0
    
    for f in filenames:
        
        print(f"\nProcessing file: {f}")
        cube, p87p5_cube = load_valid_cube(f)
            
        cube_regridded = cube.regrid(template_cube, iris.analysis.Linear())
        data_cube, _ = ThreeD_lat_land_mask(cube_regridded, lat_minimum, lat_maximum)
        mask = (p87p5_cube.data < 0) 
        masked_data = np.where(mask, np.nan, data_cube.data)
            
        data_cube_masked = data_cube.copy()
        data_cube_masked.data = masked_data
        
        tau_matrix, p_matrix, lat_points, lon_points = compute_kendall_tau_models(years, data_cube_masked)         
        save_results(tau_matrix, p_matrix, lat_points, lon_points, f)
        save_positive_p_values(tau_matrix, p_matrix, lat_points, lon_points, f)  
            
        processed += 1

    print(f"Finished processing. Total files successfully processed: {processed}/{len(filenames)}")

# Run the processing
process_model_files(Model_filenames)







