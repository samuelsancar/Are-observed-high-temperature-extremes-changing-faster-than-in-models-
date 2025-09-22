"""

Script to calculate tailwidths of daily maximum 2m temperature from JRA-55 reanalysis data.

To run: python JRA-55_tailwidths.py

"""

import iris
import numpy as np
import sys
import iris.analysis.stats
import matplotlib.pyplot as plt
import iris.plot as iplt
import iris.quickplot as qplt
import os

datadir = '/data/users/samuel.sancar/Heat_extremes/JRA-55_tailwidths/Yearly_daily_files/'
outdir = '/data/users/samuel.sancar/Heat_extremes/JRA-55_tailwidths/Yearly_tailwidths/'
template = iris.load_cube('template_cube_1deg.nc')

# Set data, variable
reanalysis = 'jra55'
variable = '_minmax_surf.015_tmax.reg_tl319.'

os.makedirs(outdir, exist_ok=True)

# Calculate tailwidths over years t1-t2
for y in np.arange(1958,2023):

    fn = reanalysis+variable+'_'+str(y)+'_tailwidth.nc'
    year_cubes = iris.load_cube(datadir+fn)

    #Calculate p99-p87.5 for this year:
    p99 = year_cubes.collapsed('time',iris.analysis.PERCENTILE,percent=99)
    p87p5 =year_cubes.collapsed('time',iris.analysis.PERCENTILE,percent=87.5)
    tail_width = p99-p87p5

    # Save each year's cube separately
    iris.save(tail_width, outdir+fn)
    print('Saved: ', outdir+fn)


