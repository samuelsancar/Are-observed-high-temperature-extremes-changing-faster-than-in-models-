
import iris
import numpy as np
import iris.analysis as ia
import iris.analysis.cartography as iac
import iris.coord_categorisation as icc
from iris.coord_categorisation import add_year
from cf_units import Unit
from Lat_Land_mask_function import lat_land_mask, ThreeD_lat_land_mask

def get_area_weights(cube):
    for coord in ['latitude', 'longitude']:
        if not cube.coord(coord).has_bounds():
            cube.coord(coord).guess_bounds()
    return iac.area_weights(cube)

def drop_problem_coords(cube):
    """Remove forecast-related coords that can block operations."""
    for bad in ('valid_time', 'forecast_reference_time'):
        if cube.coords(bad):
            cube.remove_coord(bad)
    return cube

def relabel_time_to_year(cube):
    years = cube.coord('time').points.astype(int)
    year_dim = iris.coords.DimCoord(years, long_name='year', units='1')
    out = iris.cube.Cube(
        np.ma.masked_invalid(np.asarray(cube.data)),
        long_name=cube.long_name or cube.name() or 'regional annual mean tail width',
        units=cube.units,
        attributes=cube.attributes
    )
    out.add_dim_coord(year_dim, 0)
    return out

def areaMean(cube):
    cube.data = np.ma.masked_invalid(cube.data)
    collapse_dims = [dim for dim in ['latitude', 'longitude'] if dim in [coord.name() for coord in cube.coords()]]
    if collapse_dims:
        weights = get_area_weights(cube)
        if np.ma.is_masked(cube.data):
            weights = np.ma.masked_where(cube.data.mask, weights)
        cube = cube.collapsed(collapse_dims, ia.MEAN, weights=weights)
    return cube

def area_mean_over_latlon_per_timestep(cube):

    # mask invalids to avoid propagating NaNs
    cube.data = np.ma.masked_invalid(cube.data)

    # Only collapse if both lat & lon are present
    coord_names = [c.name() for c in cube.coords()]
    collapse_dims = [dim for dim in ("latitude", "longitude") if dim in coord_names]
    if not collapse_dims:
        return cube

    weights = get_area_weights(cube)

    # Align weights to masked data so masked grid cells don't contribute
    if np.ma.is_masked(cube.data):
        weights = np.ma.masked_where(cube.data.mask, weights)

    # Collapses only spatial dims; 'time' remains
    return cube.collapsed(collapse_dims, ia.MEAN, weights=weights)

def regional_tailwidths(member_name, reglats, reglons):
    period = list(range(1958, 2023))
    annual_values = []
    units_out = None

    for year in period:
        # Load cube for this year
        if member_name == 'era5':
            cube = iris.load_cube(f'/data/users/samuel.sancar/Heat_extremes/ERA5_tailwidths/era5_daily_maximum_2m_temperature_{year}_tail_width.nc')
            cube, _ = lat_land_mask(cube, -60, 90)
        elif member_name == 'jra55':
            cube = iris.load_cube(f'/data/users/samuel.sancar/Heat_extremes/JRA-55_tailwidths/jra55__minmax_surf.015_tmax.reg_tl319._{year}_tail_width.nc')
            cube, _ = lat_land_mask(cube,-60, 90)
        else:
            cube = iris.load_cube(f'/data/users/samuel.sancar/Heat_extremes/Tail_widths/canari/member{member_name}/tailwidth_member{member_name}_{year}.nc')
            cube, _ = lat_land_mask(cube, -60, 90)    
            
        # Drop problematic coords
        cube  = drop_problem_coords(cube)

        if units_out is None:
            units_out = cube.units

        # Subset region and compute area mean
        cube_reg = cube.intersection(latitude=reglats, longitude=reglons)
        cube_reg_mean = areaMean(cube_reg)

        # Collapse any remaining dimensions (e.g., level)
        while cube_reg_mean.ndim > 0:
            cube_reg_mean = cube_reg_mean.collapsed(cube_reg_mean.coords()[0].name(), ia.MEAN)

        val = cube_reg_mean.data.item()
        annual_values.append(val)

    # Build 1-D cube
    data = np.ma.masked_invalid(np.asarray(annual_values))
    out = iris.cube.Cube(data, long_name=f'{member_name} regional annual mean tail width', units=units_out)
    year_dim = iris.coords.DimCoord(period, long_name='year', units='year')
    out.add_dim_coord(year_dim, 0)

    return out

def MIP_regional_tailwidths(member_name, reglats, reglons,
                            path_template='/data/users/samuel.sancar/Heat_extremes/Tail_widths/{member}_tailwidths.nc'):

    print(f"Processing member: {member_name}")
    cube = iris.load_cube(path_template.format(member=member_name))
    print("Input shape:", cube.shape)

    cube = drop_problem_coords(cube)

    cube_reg = cube.intersection(latitude=reglats, longitude=reglons)
    
    print("Original cube shape:", cube.shape)
    print("Subset region shape:", cube_reg.shape)

    reg_series = area_mean_over_latlon_per_timestep(cube_reg)
    print("After spatial collapse:", reg_series.shape)
    print("Coords:", [c.name() for c in reg_series.coords(dim_coords=True)])
    print(reg_series.coord('time').points[:10])

    yearly = relabel_time_to_year(reg_series)
    yearly.long_name = f'{member_name} regional annual mean tail width'
    print("Final yearly shape (expected (~65,)):", yearly.shape)

    return yearly










