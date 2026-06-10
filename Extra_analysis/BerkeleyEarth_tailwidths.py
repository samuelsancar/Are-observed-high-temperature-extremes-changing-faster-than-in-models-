""" 
10/05/2026
----------

Code to calculate BerkleyEarth gridded daily maximum temperature tail widths.

To run: python BerkeleyEarth_tailwidths.py
"""

import iris
import numpy as np
import sys
import iris.analysis.stats
import matplotlib.pyplot as plt
import iris.plot as iplt
import iris.quickplot as qplt
import os

#Beginning & end year+1
t1 = 1958
t2 = 2023

datadir = '/data/users/samuel.sancar/BerkeleyEarth/'
outdir = '/data/users/samuel.sancar/Heat_extremes/BerkeleyEarth_tailwidths/'

# Load template grid for regridding
template = iris.load_cube('template_cube_1deg.nc')

os.makedirs(outdir, exist_ok=True)

# Load climatology for each day of the year
constraint = iris.Constraint(name = 'Air Surface Temperature Climatology (Jan 1951 - Dec 1980)')
file_pattern = f"{datadir}Complete_TMAX_Daily_LatLong1_1960.nc"

example = iris.load(file_pattern)
print(example)

climatology = iris.load_cube(file_pattern, constraint)

print(climatology)

#iris.util.equalise_attributes(climatology)
#iris.util.unify_time_units(climatology)

# Calculate tailwidths over years t1-t2
for y in np.arange(t1,t2):
    
    print(f'---Processing year {y}---')
    
    if 1950 <= y < 1960:
        decade = 1950
        
    elif 1960 <= y < 1970:
        decade = 1960
        
    elif 1970 <= y < 1980:
        decade = 1970
        
    elif 1980 <= y < 1990:
        decade = 1980
        
    elif 1990 <= y < 2000:
        decade = 1990
        
    elif 2000 <= y < 2010:
        decade = 2000
        
    elif 2010 <= y < 2020:
        decade = 2010
        
    elif 2020 <= y < 2030:
        decade = 2020
        
    name_constraint = iris.Constraint(name = 'Air Surface Temperature Anomaly')
    anomaly_file = f"{datadir}Complete_TMAX_Daily_LatLong1_{decade}.nc"
    
    cubes = iris.load(anomaly_file)

    # Extract the relevant cubes
    anomaly_cube = cubes.extract_cube('surface_temperature_anomaly')
    year_cube    = cubes.extract_cube('Year')
    doy_cube     = cubes.extract_cube('Day of Year')
    
    years = np.asarray(year_cube.core_data()).astype(int)    
    doy = np.asarray(doy_cube.core_data()).astype(int) 
    
    mask_year = years == y
    
    # To handle leap year 
    mask = mask_year & (doy<=365)
    
    yearly_anomaly = anomaly_cube[mask, ...]
    yearly_doy = doy[mask]
    clim_for_year = climatology[yearly_doy - 1, ...]

    # Add climatology and anomaly for absolute Tmax 
    yearly_tmax = yearly_anomaly + clim_for_year

    print(yearly_tmax)

    time_dim = yearly_tmax.dim_coords[0]
    
    #Calculate p99-p87.5 for this year:

    data = yearly_tmax.core_data()
    
    if np.ma.isMaskedArray(data):
        data = data.filled(np.nan)
        
    data = np.array(data, copy=True)

    p99_data  = np.nanpercentile(data, 99, axis=0)
    p87p5_data = np.nanpercentile(data, 87.5, axis=0)
    tail_width_data  = p99_data - p87p5_data

    tail_width = yearly_anomaly[0, ...].copy(data=(tail_width_data).copy())
    p87p5 = yearly_anomaly[0, ...].copy(data=p87p5_data)
    
    print(tail_width)
    
    #regrid 
    tail_width_rg = tail_width.regrid(template,iris.analysis.Linear())
    tail_width_rg.rename('Yearly_tasmax_tailwidth')
    p87p5_regrid = p87p5.regrid(template,iris.analysis.Linear())
    p87p5_regrid.rename('Tamsax_p87p5')
    
    for attr in ['valid_min', 'valid_max', 'actual_range']:
        tail_width_rg.attributes.pop(attr, None)
        
    for attr in ['valid_min', 'valid_max', 'actual_range']:
        p87p5_regrid.attributes.pop(attr, None)

    # Save each year's cube separately
    out_file = os.path.join(outdir, f"BerkeleyEarth_daily_maximum_2m_temperature_{y}_tail_width.nc")
    outfile_87p5 = os.path.join(outdir, f"BerkeleyEarth_daily_maximum_2m_temperature_{y}_87p5.nc")
    iris.save(tail_width_rg, out_file)
    iris.save(p87p5_regrid, outfile_87p5)
    
    print(f'---Succesfully processed and saved year {y}---')