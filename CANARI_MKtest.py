"""
    Code to calculate KendallTau and p-values in CANARI trends:

    To run: python Model_sig_trends

"""

import iris
import glob
import numpy as np
import os
import sys
import xarray as xr
import matplotlib.pyplot as plt
from Lat_Land_mask_function import ThreeD_lat_land_mask
from MK_functions import compute_kendall_tau_models

# Define the analysis period
years = np.arange(1958, 2023)
lat_minimum = -60
lat_maximum = 90

# Load the template cube for regridding
template_cube = iris.load_cube('template_cube_1deg.nc')

# Define input and output directories
Model_dir = '/data/users/samuel.sancar/Heat_extremes/Tail_widths/canari/*.nc'
correlation = iris.load_cube('/data/users/samuel.sancar/Heat_extremes/Trends/Reanalysis_trends/ERA5_JRA-55_correlation.nc')
outdirp = '/data/users/samuel.sancar/Heat_extremes/Trends/Trend_significance/CANARI/P-vals/'
outdirMK = '/data/users/samuel.sancar/Heat_extremes/Trends/Trend_significance/CANARI/P-vals/'
outdir_pos_pvals = '/data/users/samuel.sancar/Heat_extremes/Trends/Trend_significance/CANARI/Positive_Pvals/'  # Positive trend p-vals

# Get filenames
Model_filenames = sorted(glob.glob(Model_dir, recursive=True))
print(f"Found {len(Model_filenames)} files to process.")

def save_results(tau_matrix, p_matrix, lat_points, lon_points, filename):

    base_name = os.path.splitext(os.path.basename(filename))[0]
    tau_save = xr.DataArray(tau_matrix, coords={"latitude": lat_points, "longitude": lon_points}, dims=["latitude", "longitude"], name = "Tailwidth_mktau")
    p_save = xr.DataArray(p_matrix, coords={"latitude": lat_points, "longitude": lon_points}, dims=["latitude", "longitude"], name = "Tailwidth_trend_p-values")
    tau_save.to_netcdf(os.path.join(outdirMK, f"{base_name}_mktau_1958-2022.nc"))
    p_save.to_netcdf(os.path.join(outdirp, f"{base_name}_p_1958.nc"))
    
    print(f"Saved results for: {filename}") 

def save_positive_p_values(tau_matrix, p_matrix, lat_points, lon_points, filename):

    pos_p_matrix = np.where(tau_matrix > 0, p_matrix, np.nan)
    pos_p_save = xr.DataArray(pos_p_matrix, coords={"latitude": lat_points, "longitude": lon_points}, dims=["latitude", "longitude"], name =  'Positive_trend_p-values')
    base_name = os.path.splitext(os.path.basename(filename))[0]
    pos_p_save.to_netcdf(os.path.join(outdir_pos_pvals, f"{base_name}_pos_p_1958.nc"))
    
    print(f"Saved positive trend p-values for: {filename}")


def load_valid_cube(filepath):

    cubes = iris.load(filepath)
        

    member = filepath.split('/')[-1].split('_')[0]

    p87p5_dir = f'/data/users/samuel.sancar/Heat_extremes/Trends/CANARI_trends/'
    p87p5_cube = iris.load_cube(p87p5_dir + f'{member}_p87p5_trend.nc')

    # Relaxed check: only require time dimension and 3D shape
    for cube in cubes:
        if cube.ndim == 3 and cube.coord('time'):
            return cube, p87p5_cube

def process_model_files(filenames):
    processed = 0
    for f in filenames:
        print(f"\nProcessing file: {f}")

        cube, p87p5_cube = load_valid_cube(f)

        cube_regridded = cube.regrid(template_cube, iris.analysis.Linear())
        print(f"Shape of cube {cube_regridded.shape}")
        data_cube, _ = ThreeD_lat_land_mask(cube_regridded, lat_minimum, lat_maximum)
        mask = (p87p5_cube.data < 0) #| (correlation.data < 0.5)
        masked_data = np.where(mask, np.nan, data_cube.data)
            
        data_cube_masked = data_cube.copy()
        data_cube_masked.data = masked_data
            
        tau_matrix, p_matrix, lat_points, lon_points = compute_kendall_tau_models(years, data_cube_masked)         
        save_results(tau_matrix, p_matrix, lat_points, lon_points, f)
        save_positive_p_values(tau_matrix, p_matrix, lat_points, lon_points, f)  # <-- Saving positive p-values
            
        processed += 1

    print(f"Finished processing. Total files successfully processed: {processed}/{len(filenames)}")

# Run the processing
process_model_files(Model_filenames)
