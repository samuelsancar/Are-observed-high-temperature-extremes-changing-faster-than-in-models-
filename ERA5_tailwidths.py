"""
17/07/25
---------
Code to merge monthly ERA5 data and calculate tail width (p99-p87.5) in daily tasmax.

"""

import iris
import numpy as np
import sys
import iris.analysis.stats
import matplotlib.pyplot as plt
import iris.plot as iplt
import iris.quickplot as qplt
import os

#Beginning & end year
t1 = 1958
t2 = 2023

datadir = '/data/users/appldata/Data/OBS-ERA5/daily/2m_temperature/daily_maximum/'
outdir = '/data/users/samuel.sancar/Heat_extremes/ERA5_tailwidths/'

# Load template grid for regridding
template = iris.load_cube('template_cube_1deg.nc')

reanalysis = 'era5'
variable = 'daily_maximum_2m_temperature'

os.makedirs(outdir, exist_ok=True)

# Calculate tailwidths over years t1-t2
for y in np.arange(t1,t2):
    
    file_pattern = f"{datadir}{reanalysis}_{variable}_{y}*.nc"
    cube_list = iris.load(file_pattern)

    iris.util.equalise_attributes(cube_list)
    iris.util.unify_time_units(cube_list)

    # Merge monthly files into yearly cube
    merged_cube = cube_list.concatenate_cube()

    #Calculate p99-p87.5 for this year:
    p99 = merged_cube.collapsed('valid_time',iris.analysis.PERCENTILE,percent=99)
    p87p5 = merged_cube.collapsed('valid_time',iris.analysis.PERCENTILE,percent=87.5)

    tail_width = p99-p87p5

    #regrid 
    tail_width_rg = tail_width.regrid(template,iris.analysis.Linear())
    tail_width_rg.rename('Yearly_tasmax_tailwidth')
    p87p5_regrid = p87p5.regrid(template,iris.analysis.Linear())
    p87p5_regrid.rename('Tamsax_p87p5')

    # Save each year's cube separately
    out_file = os.path.join(outdir, f"{reanalysis}_{variable}_{y}_tail_width.nc")
    outfile_87p5 = os.path.join(outdir, f"{reanalysis}_{variable}_{y}_87p5.nc")
    iris.save(tail_width_rg, out_file)
    iris.save(p87p5_regrid, outfile_87p5)
    

