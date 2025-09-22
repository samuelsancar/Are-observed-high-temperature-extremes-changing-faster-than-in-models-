"""
    Function to apply land and latitude mask to iris cubes. Latitude mask can be chosen.

    Returns masked_cube, no_of_gridpoints

"""

import numpy as np
import iris
import iris.coords
import iris.cube

def lat_land_mask(model, latitude_threshold):
    # Load template and land-sea mask
    template = iris.load_cube('template_cube_1deg.nc')
    lsm = iris.load_cube('lsm_1279l4_0.1x0.1.grb_v4_unpack.nc')
    lsm_regridded = lsm.regrid(template, iris.analysis.Nearest())

    # Create land mask
    land_mask = np.squeeze(lsm_regridded.data) > 0.5

    # Create latitude mask
    latitudes = template.coord('latitude').points
    lat_mask = latitudes >= latitude_threshold
    lat_mask_2d = np.repeat(lat_mask[:, np.newaxis], land_mask.shape[1], axis=1)

    # Combine masks
    combined_mask = land_mask & lat_mask_2d
    no_of_gridpoints = np.count_nonzero(combined_mask)

    # Handle model input: Iris cube or NumPy array
    if isinstance(model, iris.cube.Cube):
        model_data = model.data
        units = model.units
        attributes = model.attributes
    elif isinstance(model, np.ndarray):
        model_data = model
        units = '1'  # Default unit for agreement fraction
        attributes = {}
    else:
        raise TypeError("Input 'model' must be an Iris cube or a NumPy array.")

    # Apply mask
    masked_data = np.where(combined_mask, model_data, np.nan)

    # Wrap masked data into an Iris cube
    lat_coord = template.coord('latitude')
    lon_coord = template.coord('longitude')
    masked_cube = iris.cube.Cube(
        masked_data.astype(np.float32),
        dim_coords_and_dims=[(lat_coord, 0), (lon_coord, 1)],
        long_name='tail_width_trend_masked',
        units=units,
        attributes=attributes
    )

    return masked_cube, no_of_gridpoints

def lat_land_mask(model, latitude_threshold):
    # Load template and land-sea mask
    template = iris.load_cube('template_cube_1deg.nc')
    lsm = iris.load_cube('lsm_1279l4_0.1x0.1.grb_v4_unpack.nc')
    lsm_regridded = lsm.regrid(template, iris.analysis.Nearest())

    # Create land mask
    land_mask = np.squeeze(lsm_regridded.data) > 0.5

    # Create latitude mask
    latitudes = template.coord('latitude').points
    lat_mask = latitudes >= latitude_threshold
    lat_mask_2d = np.repeat(lat_mask[:, np.newaxis], land_mask.shape[1], axis=1)

    # Combine masks
    combined_mask = land_mask & lat_mask_2d
    no_of_gridpoints = np.count_nonzero(combined_mask)

    # Handle model input: Iris cube or NumPy array
    if isinstance(model, iris.cube.Cube):
        model_data = model.data
        units = model.units
        attributes = model.attributes
    elif isinstance(model, np.ndarray):
        model_data = model
        units = '1'  # Default unit for agreement fraction
        attributes = {}
    else:
        raise TypeError("Input 'model' must be an Iris cube or a NumPy array.")

    # Apply mask
    masked_data = np.where(combined_mask, model_data, np.nan)

    # Wrap masked data into an Iris cube
    lat_coord = template.coord('latitude')
    lon_coord = template.coord('longitude')
    masked_cube = iris.cube.Cube(
        masked_data.astype(np.float32),
        dim_coords_and_dims=[(lat_coord, 0), (lon_coord, 1)],
        long_name='tail_width_trend_masked',
        units=units,
        attributes=attributes
    )

    return masked_cube, no_of_gridpoints

