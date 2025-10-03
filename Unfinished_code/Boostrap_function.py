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

def bootstrap_tailwidth_trends(data_cube, iterations=1000, alpha=0.05):
    
    if data_cube.ndim != 3:
        raise ValueError("Expected a 3D cube with dimensions (time, lat, lon)")

    # Add 'year' coordinate from 'time'
    iris.coord_categorisation.add_year(data_cube, 'time', name='year')
    years = data_cube.coord('year').points

    lat_points = data_cube.coord('latitude').points
    lon_points = data_cube.coord('longitude').points

    # Prepare output arrays
    slopes_matrix = np.empty((len(lat_points), len(lon_points)))
    quantiles_matrix = np.empty((len(lat_points), len(lon_points), 4))  # For 4 quantiles
    sig_matrix = np.empty((len(lat_points), len(lon_points)))

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
            sig_matrix[i, j] = ((quantiles_matrix[i, j] > (alpha/2) and slopes_matrix[i, j] > 0) or (quantiles_matrix[i, j] < (1-alpha/2) and slopes_matrix[i, j] < 0))
                

    return sig_matrix, slopes_matrix, quantiles_matrix



    
    