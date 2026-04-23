"""
Code to calculate correlation between ERA5 and MERRA
"""

import iris
import os
from iris.analysis import stats
import iris.util
import numpy as np
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.cm as mpl_cm
from Lat_Land_mask_function import lat_land_mask
from Map_features import map_features
import matplotlib.pyplot as plt
import iris.plot as iplt
from iris.cube import CubeList
from Trend_calculation_functions import create_time_coord, load_yearly_cubes, stack_cubes

datadir = ['/data/users/samuel.sancar/Heat_extremes/MERRA_tailwidths/', '/data/users/samuel.sancar/Heat_extremes/ERA5_tailwidths/']
p87p5 = ['/data/users/samuel.sancar/Heat_extremes/MERRA_tailwidths/', '/data/users/samuel.sancar/Heat_extremes/ERA5_tailwidths/']
outdir = ['/data/users/samuel.sancar/Heat_extremes/Trends/Reanalysis_trends/']

reanalysis = ['MERRA', 'era5']
variable = ['daily_maximum_2m_temperature', 'daily_maximum_2m_temperature']
t1 = 1980
t2 = 2023
years = np.arange(t1, t2)

def MERRA_load_yearly_cubes(datadir, p87p5, reanalysis, variable, years):
    
    cubes = iris.cube.CubeList()
    for year in years:
        filename = f"{reanalysis}_{variable}_{year}_tail_width.nc"
        cube = iris.load_cube(os.path.join(datadir, filename))
        cubes.append(cube)
        
    p87p5_cubes = iris.cube.CubeList()
    for year in years:
        p87p5_filename = f"{variable}_{year}_87p5.nc"
        p87p5_cube = iris.load_cube(os.path.join(p87p5, p87p5_filename))
        p87p5_cubes.append(p87p5_cube)
        
    return cubes, p87p5_cubes

def MERRA_prepare(reanalysis, variable, datadir, p87p5):
    
    time_coord = create_time_coord(years)
    cubes, p87p5_cubes = MERRA_load_yearly_cubes(datadir, p87p5 , reanalysis, variable, years)
    data_cube, p87p5_data_cube = stack_cubes(cubes, p87p5_cubes, time_coord)
     
     # Ensure time coordinate is attached
    if not data_cube.coords('time'):
        data_cube.add_dim_coord(time_coord, 0)
    if not p87p5_data_cube.coords('time'):
        p87p5_data_cube.add_dim_coord(time_coord, 0)
        
    return data_cube

def prepare_reanlysis(reanalysis, variable, datadir, p87p5):
    time_coord = create_time_coord(years)
    cubes, p87p5_cubes = load_yearly_cubes(datadir, p87p5 , reanalysis, variable, years)
    data_cube, p87p5_data_cube = stack_cubes(cubes, p87p5_cubes, time_coord)
     
     # Ensure time coordinate is attached
    if not data_cube.coords('time'):
        data_cube.add_dim_coord(time_coord, 0)
    if not p87p5_data_cube.coords('time'):
        p87p5_data_cube.add_dim_coord(time_coord, 0)
        
    return data_cube
        
merra = MERRA_prepare(reanalysis[0], variable[0], datadir[0], p87p5[0])
era5 = prepare_reanlysis(reanalysis[1], variable[1], datadir[1], p87p5[1])

def sort_by_time(cube):
    """Sort cube along its 'time' dimension using NumPy (works across Iris versions)."""
    t = cube.coord('time')
    dims = cube.coord_dims('time')
    if len(dims) != 1:
        raise ValueError("Ambiguous 'time' dims")
    time_axis = dims[0]
    order = np.argsort(t.points)
    if np.all(order == np.arange(order.size)):
        return cube  # already sorted
    index = [slice(None)] * cube.ndim
    index[time_axis] = order
    return cube[tuple(index)]

def ensure_time_is_dim(c):
    """Promote 'time' to a dimension coordinate if needed."""
    dims = c.coord_dims('time')
    assert len(dims) == 1, "Ambiguous 'time' dims"
    time_coord = c.coord('time')
    if time_coord not in c.dim_coords:
        axis = dims[0]
        c2 = c.copy()
        c2.remove_coord('time')
        c2.add_dim_coord(time_coord.copy(), axis)
        return c2
    return c

def adopt_dim_coord(target, template, name, axis):
    """
    Force the dimension coordinate object on `target` to be an exact copy of `template`'s `name` coord,
    placed at the given axis. This makes coord identity match (not just values).
    """
    tgt = target.copy()
    if tgt.coords(name):
        tgt.remove_coord(name)
    # copy preserves units, bounds, coord_system, etc.
    tgt.add_dim_coord(template.coord(name).copy(), axis)
    return tgt

def corr_along_time_np(a, b):
    """
    Pearson r along axis 0 (time) with NaN handling.
    a, b shapes: (T, Y, X)
    Returns (Y, X).
    """
    # joint finite mask
    mask = np.isfinite(a) & np.isfinite(b)
    a = np.where(mask, a, np.nan)
    b = np.where(mask, b, np.nan)

    a_mean = np.nanmean(a, axis=0, keepdims=True)
    b_mean = np.nanmean(b, axis=0, keepdims=True)
    A = a - a_mean
    B = b - b_mean

    num = np.nansum(A * B, axis=0)
    den = np.sqrt(np.nansum(A * A, axis=0) * np.nansum(B * B, axis=0))
    with np.errstate(invalid='ignore', divide='ignore'):
        r = num / den
    return r

iris.util.unify_time_units([era5, merra])

cubes = CubeList([era5, merra])
iris.util.equalise_attributes(cubes)
era5, merra = cubes

try:
    if str(merra.units).lower() in ("unknown", "1", "unknown unit", "unknown units"):
        merra.units = era5.units
    elif merra.units != era5.units:
        merra.convert_units(era5.units)
except Exception:
    merra.units = era5.units

era5 = ensure_time_is_dim(era5)
merra  = ensure_time_is_dim(merra)

era5 = sort_by_time(era5)
merra  = sort_by_time(merra)

merra = adopt_dim_coord(merra, era5, 'time',      axis=0)
merra = adopt_dim_coord(merra, era5, 'latitude',  axis=1)
merra = adopt_dim_coord(merra, era5, 'longitude', axis=2)

joint_mask = np.isfinite(era5.data) & np.isfinite(merra.data)
era5.data = np.where(joint_mask, era5.data, np.nan)
merra.data  = np.where(joint_mask, merra.data,  np.nan)

r = iris.analysis.stats.pearsonr(era5, merra, corr_coords='time')

directory = '/data/users/samuel.sancar/Heat_extremes/Trends/Reanalysis_trends/'
title = 'ERA5_MERRA_correlation.nc'
outdir = os.path.join(directory, title)
iris.save(r, outdir)

# Apply latitude and land mask using function
model_mask, _ = lat_land_mask(r, -60, 90)

# Plot trends
fig = plt.figure(figsize=(14, 8))
ax = plt.axes(projection=ccrs.Robinson(central_longitude=0))
cmap_rb = mpl_cm.get_cmap('coolwarm', 8)

map = iplt.pcolormesh(model_mask, cmap=cmap_rb, vmin=-1, vmax=1)

map_features('Rob', ax, -60, 90)

cbar = plt.colorbar(map, ax=ax, orientation='horizontal', pad=0.05, aspect=50)
cbar.ax.tick_params(labelsize=14)
cbar.set_label('Correlations', fontsize=15)

plt.title(f"MERRA and ERA5 tailwidth correlations ({t1}-{t2-1})", fontsize=18, weight='bold')

plt.tight_layout()
plt.show()