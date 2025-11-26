"""
    
    Perform bootstrapped linear trend estimation of tailwidth over time at each grid point.

    Parameters:
        data_cube (iris.cube.Cube): A 3D cube with dimensions (time, latitude, longitude).
        iterations (int): Number of bootstrap iterations.
        alpha int: Significance level for trend significance testing.

    Returns:
        slopes_matrix (np.ndarray): A 2D array of mean bootstrapped slopes for each grid point.
        quantiles_matrix (np.ndarray): A 3D array of quantiles [2.5%, 5%, 95%, 97.5%] for each grid point.
    
"""

import iris
import iris.coord_categorisation
import numpy as np
import os
import xarray as xr

def bootstrap_tailwidth_trends(data_cube, iterations=1000, alpha=0.05):
    
    if data_cube.ndim != 3:
        raise ValueError("Expected a 3D cube with dimensions (time, lat, lon)")

    # Add 'year' coordinate from 'time'
    iris.coord_categorisation.add_year(data_cube, 'time', name='year')
    years = data_cube.coord('year').points

    lat_points = data_cube.coord('latitude').points
    lon_points = data_cube.coord('longitude').points

    # Output arrays
    slopes_matrix = np.empty((len(lat_points), len(lon_points)))
    quantiles_matrix = np.empty((len(lat_points), len(lon_points), 4))  # For 4 quantiles
    b_matrix = np.empty((len(lat_points), len(lon_points)))

    rng = np.random.default_rng()

    for i in range(len(lat_points)):
        for j in range(len(lon_points)):
            tailwidths = data_cube.data[:, i, j]

            # Skip grid points with missing or constant data
            if np.ma.is_masked(tailwidths) or np.all(tailwidths == tailwidths[0]):
                slopes_matrix[i, j] = np.nan
                quantiles_matrix[i, j, :] = np.nan
                continue

            slopes = []
            for _ in range(iterations):
                idxs = rng.choice(len(years), len(years), replace=True)
                sampled_years = years[idxs]
                sampled_tailwidths = tailwidths[idxs]
                coeffs = np.polyfit(sampled_years, sampled_tailwidths, 1)
                slopes.append(coeffs[0])

            # Store mean slope and quantiles
            slopes_matrix[i, j] = np.mean(slopes)
            quantiles_matrix[i, j, :] = np.quantile(slopes, [0.025, 0.05, 0.95, 0.975])
            b_matrix[i, j] = ((np.any((quantiles_matrix > (alpha/2))<0) and np.any(slopes_matrix > 0)) or np.any((quantiles_matrix < (1-alpha/2))>0) and np.any(slopes_matrix < 0))

                

    return b_matrix, slopes_matrix, quantiles_matrix, lat_points, lon_points

def save_results_b_mask(b_matrix, lat_points, lon_points):      
    p_da = xr.DataArray(b_matrix, coords={"latitude": lat_points, "longitude": lon_points}, dims=["latitude", "longitude"])
    outdir = '/home/users/samuel.sancar/Heat_extremes/Trends/Reanalysis_trends/'
    os.makedirs(outdir, exist_ok=True)
    p_da.to_netcdf(os.path.join(outdir, f"era5_daily_maximum_2m_temperature_tailwidth_b_mask_1958-2023.nc"))





    
    