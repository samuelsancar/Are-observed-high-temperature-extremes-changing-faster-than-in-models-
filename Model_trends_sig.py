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
from Lat_Land_mask_function import lat_land_mask
import sys

# Define the analysis period
years = np.arange(1958, 2023)

# Load the template cube for regridding
template_cube = iris.load_cube('template_cube_1deg.nc')

# Define output directory
outdir = '/data/users/samuel.sancar/Heat_extremes/trends/Trend_significance/'
os.makedirs(outdir, exist_ok=True)
print(f"Output directory is set to: {outdir}")

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
    tau_path = os.path.join(outdir, f"{base_name}_tailwidth_mktau_1958-2022.nc")
    p_path = os.path.join(outdir, f"{base_name}_tailwidth_p_1958.nc")
    tau_da.to_netcdf(tau_path)
    p_da.to_netcdf(p_path)
    print(f"Saved: {tau_path}")
    print(f"Saved: {p_path}")

def process_single_model_file(filename):
    print(f"\nProcessing file: {filename}")
    try:
        cube = iris.load_cube(filename)
    except Exception as e:
        print(f"Failed to load cube from {filename}: {e}")
        return

    if cube.ndim != 3:
        print(f"Skipping {filename}: not a 3D cube")
        return

    try:

        # Regrid
        cube_regridded = cube.regrid(template_cube, iris.analysis.Linear())

        # Fix coordinate dimension assignments if needed
        if cube_regridded.coord_dims('latitude') != (1,):
            cube_regridded.remove_coord('latitude')
            cube_regridded.add_dim_coord(template_cube.coord('latitude'), 1)

        if cube_regridded.coord_dims('longitude') != (2,):
            cube_regridded.remove_coord('longitude')
            cube_regridded.add_dim_coord(template_cube.coord('longitude'), 2)

        masked_cube, _ = lat_land_mask(cube_regridded, -60)
        tau_matrix, p_matrix, lat_points, lon_points = compute_kendall_tau(years, masked_cube)
        save_results(tau_matrix, p_matrix, lat_points, lon_points, outdir, filename)
        print("Finished processing the file successfully.")
    except Exception as e:
        print(f"Error during processing of {filename}: {e}")

# Entry point for script execution
if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python Model_sig_trends.py <model_file_path>")
    else:
        model_file = f"/data/users/samuel.sancar/Heat_extremes/Tail_widths/{sys.argv[1]}"
        process_single_model_file(model_file)




