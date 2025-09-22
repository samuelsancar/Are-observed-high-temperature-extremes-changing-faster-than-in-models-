"""
    Code to calculate KendallTau and p-values in model trends:

    To run: python Model_sig_trends

"""

import iris
import glob
import numpy as np
import os
import xarray as xr
import scipy.stats as stats
from Lat_Land_mask_function import ThreeD_lat_land_mask

# Define the analysis period
years = np.arange(1958, 2023)

# Load the template cube for regridding
template_cube = iris.load_cube('template_cube_1deg.nc')

# Define input and output directories
coupled_dir = '/data/users/samuel.sancar/Heat_extremes/Tail_widths/*.nc'
outdir = '/data/users/samuel.sancar/Heat_extremes/trends/Trend_significance/'
os.makedirs(outdir, exist_ok=True)
print(f"Output directory is set to: {outdir}")

# Get filenames
coupled_filenames = sorted(glob.glob(coupled_dir, recursive=True))
print(f"Found {len(coupled_filenames)} files to process.")

def compute_kendall_tau(years, data_cube):
    if data_cube.ndim != 3:
        raise ValueError("Expected a 3D cube with dimensions (time, lat, lon)")

    lat_points = data_cube.coord('latitude').points
    lon_points = data_cube.coord('longitude').points
    tau_matrix = np.empty((len(lat_points), len(lon_points)))
    p_matrix = np.empty((len(lat_points), len(lon_points)))

    for i in range(len(lat_points)):
        for j in range(len(lon_points)):
            data_series = data_cube.data[:, i, j]
            if np.all(np.isnan(data_series)):
                tau_matrix[i, j] = np.nan
                p_matrix[i, j] = np.nan
            else:
                tau, p = stats.kendalltau(years, data_series)
                tau_matrix[i, j] = tau
                p_matrix[i, j] = p

    return tau_matrix, p_matrix, lat_points, lon_points

def save_results(tau_matrix, p_matrix, lat_points, lon_points, outdir, filename):
    print(f"Saving results for: {filename}")
    tau_da = xr.DataArray(tau_matrix, coords={"latitude": lat_points, "longitude": lon_points}, dims=["latitude", "longitude"])
    p_da = xr.DataArray(p_matrix, coords={"latitude": lat_points, "longitude": lon_points}, dims=["latitude", "longitude"])
    base_name = os.path.splitext(os.path.basename(filename))[0]
    tau_path = os.path.join(outdir, f"{base_name}_mktau_1958-2022.nc")
    p_path = os.path.join(outdir, f"{base_name}_p_1958.nc")
    tau_da.to_netcdf(tau_path)
    p_da.to_netcdf(p_path)
    print(f"Saved: {tau_path}")
    print(f"Saved: {p_path}")

def load_valid_cube(filepath):
    try:
        cubes = iris.load(filepath)
        for cube in cubes:
            if cube.ndim == 3 and cube.shape[0] == len(years):
                return cube
        raise ValueError("No valid 3D cube with expected time dimension found.")
    except Exception as e:
        raise RuntimeError(f"Failed to load cube from {filepath}: {e}")

def process_model_files(filenames):
    processed = 0
    for f in filenames:
        print(f"\nProcessing file: {f}")
        try:
            cube = load_valid_cube(f)
        except Exception as e:
            print(e)
            continue

        try:
            cube_regridded = cube.regrid(template_cube, iris.analysis.Linear())
            masked_cube, _ = ThreeD_lat_land_mask(cube_regridded, -60, 90)
            tau_matrix, p_matrix, lat_points, lon_points = compute_kendall_tau(years, masked_cube)
            save_results(tau_matrix, p_matrix, lat_points, lon_points, outdir, f)
            processed += 1
        except Exception as e:
            print(f"Error during processing of {f}: {e}")
            continue

    print(f"Finished processing. Total files successfully processed: {processed}/{len(filenames)}")

# Run the processing
process_model_files(coupled_filenames)






