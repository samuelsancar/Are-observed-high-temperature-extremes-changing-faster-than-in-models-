"""
    Code to calculate KendallTau and p-values in model trends:

    To run: python Model_sig_trends

"""

import iris
import glob
import numpy as np
import matplotlib.pyplot as plt
import scipy.special
import os
import cf_units
import xarray as xr
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.cm as mpl_cm
import scipy.stats as stats
from Lat_Land_mask_function import lat_land_mask 

# Load and regrid model cubes
coupled_filenames = sorted(glob.glob('/data/users/samuel.sancar/Heat_extremes/trends/Coupled_model_trends/*'))
sst_filenames = sorted(glob.glob('/data/users/samuel.sancar/Heat_extremes/trends/SST_model_trends/*'))

# Load template grid for regridding
template_cube = iris.load_cube('template_cube_1deg.nc')
                       
coupled_data = [iris.load_cube(f).regrid(template_cube, iris.analysis.Linear()) for f in coupled_filenames]
sst_data = [iris.load_cube(f).regrid(template_cube, iris.analysis.Linear()) for f in sst_filenames]

years = np.arange(1958, 2023)

def create_time_coord(years):
    time_unit = cf_units.Unit('days since 1970-01-01', calendar='gregorian')
    return iris.coords.DimCoord(years, standard_name='time', long_name='valid_time', units=time_unit)

def create_year_cube(cube, years):
    lat_len = cube.shape[0]
    lon_len = cube.shape[1]
    year_data = np.empty((len(years), lat_len, lon_len), dtype=np.float32)
    for i, year in enumerate(years):
        year_data[i, :, :] = year
    year_cube = iris.cube.Cube(
        year_data,
        dim_coords_and_dims=[
            (create_time_coord(years), 0),
            (cube.coord('latitude'), 1),
            (cube.coord('longitude'), 2)
        ],
        long_name='year_cube',
        units='1'
    )
    return year_cube

def compute_kendall_tau(data_cube):
    lat_points = year_cube.coord('latitude').points
    lon_points = year_cube.coord('longitude').points
    tau_matrix = np.empty((len(lat_points), len(lon_points)))
    p_matrix = np.empty((len(lat_points), len(lon_points)))

    for i in range(len(lat_points)):
        for j in range(len(lon_points)):
            time_series = year_cube.data[:, i, j]
            data_series = data_cube.data[:, i, j]
            tau, p = stats.kendalltau(time_series, data_series)
            tau_matrix[i, j] = tau
            p_matrix[i, j] = p

    return tau_matrix, p_matrix, lat_points, lon_points

def save_results(tau_matrix, p_matrix, lat_points, lon_points, outdir, f):
    tau_da = xr.DataArray(tau_matrix, coords={"latitude": lat_points, "longitude": lon_points}, dims=["latitude", "longitude"])
    p_da = xr.DataArray(p_matrix, coords={"latitude": lat_points, "longitude": lon_points}, dims=["latitude", "longitude"])
    os.makedirs(outdir, exist_ok=True)
    tau_da.to_netcdf(os.path.join(outdir, f"{f}_tailwidth_mktau_1958-2022.nc"))
    p_da.to_netcdf(os.path.join(outdir, f"{f}_tailwidth_p_1958.nc"))


outdir = '/data/users/samuel.sancar/Heat_extremes/trends/Trend_significance/'

for f in coupled_filenames:
    cube = iris.load_cube(f).regrid(template_cube, iris.analysis.Linear())
    masked_cube, _ = lat_land_mask(cube, -60)
    time_coord = create_time_coord(years)
    year_cube = create_year_cube(masked_cube, years)
    tau_matrix, p_matrix, lat_points, lon_points = compute_kendall_tau(year_cube)
    save_results(tau_matrix, p_matrix, lat_points, lon_points, outdir, os.path.basename(f))

for f in sst_filenames:
    cube = iris.load_cube(f).regrid(template_cube, iris.analysis.Linear())
    masked_cube, _ = lat_land_mask(cube, -60)
    time_coord = create_time_coord(years)
    year_cube = create_year_cube(masked_cube, years)
    tau_matrix, p_matrix, lat_points, lon_points = compute_kendall_tau(year_cube)
    save_results(tau_matrix, p_matrix, lat_points, lon_points, outdir, os.path.basename(f))

