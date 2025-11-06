"""
24/07/25
---------
List of functions used to perform and calculate MannKendall Test data. Apply functions suitably for the nature of the data.

Functions present: load_and_mask_yearly_data, create_time_coord, stack_cubes, create_year_cube, compute_kendall_tau, save_results.

Returns: tau_matrix, p_matrix, lat_points, lon_points

"""

import os
import sys
import numpy as np
import iris
import cf_units
import scipy.stats as stats
import xarray as xr
from Lat_Land_mask_function import lat_land_mask

def load_and_mask_yearly_data(datadir, reanalysis, variable, years, lat_min, lat_max):
    cubes = iris.cube.CubeList()
    for year in years:
        filename = f"{reanalysis}_{variable}_{year}_tail_width.nc"
        cube = iris.load_cube(os.path.join(datadir, filename))
        masked_cube, _ = lat_land_mask(cube, lat_min, lat_max)
        cubes.append(masked_cube)
    return cubes

def create_time_coord(years):
    time_unit = cf_units.Unit('days since 1970-01-01', calendar='gregorian')
    return iris.coords.DimCoord(years, standard_name='time', long_name='valid_time', units=time_unit)

def stack_cubes(cubes, time_coord):
    iris.util.equalise_attributes(cubes)
    data_stack = np.stack([cube.data for cube in cubes], axis=0)
    lat_coord = cubes[0].coord('latitude')
    lon_coord = cubes[0].coord('longitude')
    return iris.cube.Cube(
        data_stack,
        dim_coords_and_dims=[(time_coord, 0), (lat_coord, 1), (lon_coord, 2)],
        long_name='tail_width_over_years',
        units=cubes[0].units,
        attributes=cubes[0].attributes
    )

def create_year_cube(cube, years):
    year_data = np.empty_like(cube.data)
    for i, year in enumerate(years):
        year_data[i, :, :] = year
    year_cube = cube.copy()
    year_cube.data = year_data
    return year_cube

def compute_kendall_tau(year_cube, data_cube):  # Use for reanalysis data
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

def compute_kendall_tau_models(years, data_cube):  # Use for model data
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

def save_results(tau_matrix, p_matrix, lat_points, lon_points, outdir, reanalysis, variable, t1, t2):
    p_da = xr.DataArray(p_matrix, coords={"latitude": lat_points, "longitude": lon_points}, dims=["latitude", "longitude"])
    os.makedirs(outdir, exist_ok=True)
    p_da.to_netcdf(os.path.join(outdir, f"{reanalysis}_{variable}_tailwidth_p_{t1}-{t2}.nc"))
    print(f"Saved p-values")
    
def save_positive_p_values(tau_matrix, p_matrix, lat_points, lon_points, outdir_pos_pvals, reanalysis, variable, t1, t2):
    pos_p_matrix = np.where(tau_matrix > 0, p_matrix, np.nan)
    pos_p_da = xr.DataArray(pos_p_matrix, coords={"latitude": lat_points, "longitude": lon_points}, dims=["latitude", "longitude"])
    filename = f"{reanalysis}{variable}_pos_p_{t1}-{t2}.nc"
    pos_p_path = os.path.join(outdir_pos_pvals, filename)
    pos_p_da.to_netcdf(pos_p_path)
    print(f"Saved positive p-values to: {pos_p_path}")


"""
'Main' example when applied to ERA5 tailwidth data.

To run: if __name__ == "__main__":
            main()

"""

def main():
    datadir = '/data/users/samuel.sancar/Heat_extremes/tailwidthslist/era5/'
    outdir = '/data/users/samuel.sancar/Heat_extremes/trends/Reanalysis_trends'
    reanalysis = 'era5'
    variable = 'daily_maximum_2m_temperature'
    t1 = int(sys.argv[1])  # need to import sys
    t2 = int(sys.argv[2])
    years = np.arange(t1, t2)

    time_coord = create_time_coord(years)
    cubes = load_and_mask_yearly_data(datadir, reanalysis, variable, years)
    data_cube = stack_cubes(cubes, time_coord)
    year_cube = create_year_cube(data_cube, years)
    tau_matrix, p_matrix, lat_points, lon_points = compute_kendall_tau(year_cube, data_cube)
    save_results(tau_matrix, p_matrix, lat_points, lon_points, outdir, reanalysis, variable, t1, t2)

