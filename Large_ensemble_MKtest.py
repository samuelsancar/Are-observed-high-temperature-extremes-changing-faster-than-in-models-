"""
    Code to calculate Mann-Kendall trend tests and significant trends of MASS large ensemble.
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

# Analysis period with missing years removed
years_full = np.arange(1958, 2023)
years = np.array([y for y in years_full if y < 2010 or y > 2014])
lat_minimum = -60
lat_maximum = sys.argv[1]  # e.g., 90

# Load the template cube for regridding
template_cube = iris.load_cube('template_cube_1deg.nc')

# Input and output directories
Model_dir = '/data/users/samuel.sancar/Heat_extremes/Tail_widths/MASS/*_tailwidths.nc'
outdir = '/data/users/samuel.sancar/Heat_extremes/Trends/Trend_significance/MASS/'

# Get filenames
Model_filenames = sorted(glob.glob(Model_dir, recursive=True))
print(f"Found {len(Model_filenames)} files to process.")

def save_results(tau_matrix, p_matrix, lat_points, lon_points, filename):
    print(f"Saving results for: {filename}")
    tau_da = xr.DataArray(tau_matrix, coords={"latitude": lat_points, "longitude": lon_points}, dims=["latitude", "longitude"])
    p_da = xr.DataArray(p_matrix, coords={"latitude": lat_points, "longitude": lon_points}, dims=["latitude", "longitude"])
    base_name = os.path.splitext(os.path.basename(filename))[0]
    tau_path = os.path.join(outdir, f"{base_name}_mktau_1958-2022_filtered.nc")
    p_path = os.path.join(outdir, f"{base_name}_p_1958_filtered.nc")
    tau_da.to_netcdf(tau_path)
    p_da.to_netcdf(p_path)
    print(f"Saved: {tau_path}")
    print(f"Saved: {p_path}")

def save_positive_p_values(tau_matrix, p_matrix, lat_points, lon_points, filename):
    print(f"Saving positive trend p-values for: {filename}")
    pos_p_matrix = np.where(tau_matrix > 0, p_matrix, np.nan)
    pos_p_da = xr.DataArray(pos_p_matrix, coords={"latitude": lat_points, "longitude": lon_points}, dims=["latitude", "longitude"])
    base_name = os.path.splitext(os.path.basename(filename))[0]
    pos_p_path = os.path.join(outdir, f"{base_name}_pos_p_1958_filtered.nc")
    pos_p_da.to_netcdf(pos_p_path)
    print(f"Saved: {pos_p_path}")

def load_valid_cube(filepath):
    iris.FUTURE.date_microseconds = True

    path = os.path.basename(filepath).replace('tailwidths', 'p87p5_trend')
    member = path.split('_')[0]
    p87p5_dir = f'/data/users/samuel.sancar/Heat_extremes/Trends/MASS_model_trends/{member}_p87p5_trend.nc'

    # Load tailwidths cubes
    cubes = iris.load(filepath)
    valid_cube = None
    for cube in cubes:

        if cube.coords('time'):
            time_coord = cube.coord('time')
            time_years = np.array([t.point.year for t in time_coord.cells()])
            valid_indices = np.where((time_years < 2010) | (time_years > 2014))[0]

            if len(valid_indices) > 0:
                cube_filtered = cube[valid_indices, :, :]
                valid_cube = cube_filtered
                break

    if valid_cube is None:
        raise ValueError("No valid cube with time coordinate found.")

    # Load p87p5 trend cube
    p87p5_cube = iris.load_cube(p87p5_dir)

    p87p5_data = np.array(p87p5_cube.data)
    if p87p5_data.dtype.kind in {'U', 'S'}:
        p87p5_data = p87p5_data.astype(float)

    cube_data = np.array(valid_cube.data)
    if cube_data.dtype.kind in {'U', 'S'}:
        cube_data = cube_data.astype(float)

    return valid_cube, p87p5_data


def plot_p_values(p_matrix, lat_points, lon_points):
    plt.figure(figsize=(10, 6))
    plt.contourf(lon_points, lat_points, p_matrix, levels=100, cmap='viridis')
    plt.colorbar(label='Kendall tau p-value')
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
    plt.title('Kendall tau p-value map')
    plt.grid(True)
    plt.tight_layout()
    plt.show()

def process_model_files(filenames):
    processed = 0
    for f in filenames:
        print(f"\nProcessing file: {f}")

        cube, p87p5_data = load_valid_cube(f)

        cube_regridded = cube.regrid(template_cube, iris.analysis.Linear())
        data_cube, _ = ThreeD_lat_land_mask(cube_regridded, lat_minimum, lat_maximum)
        masked_data = np.where(p87p5_data > 0, data_cube.data, np.nan)
        data_cube_masked = data_cube.copy()
        data_cube_masked.data = masked_data

        tau_matrix, p_matrix, lat_points, lon_points = compute_kendall_tau_models(years, data_cube_masked)
        plot_p_values(p_matrix, lat_points, lon_points)
        save_results(tau_matrix, p_matrix, lat_points, lon_points, f)
        save_positive_p_values(tau_matrix, p_matrix, lat_points, lon_points, f)
        processed += 1

    print(f"Finished processing. Total files successfully processed: {processed}/{len(filenames)}")

# Run the processing
process_model_files(Model_filenames)
