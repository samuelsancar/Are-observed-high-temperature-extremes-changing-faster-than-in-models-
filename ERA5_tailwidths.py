"""
17/07/25
---------
Code to calculate tail width (p99-p87.5) in daily tasmax, for each era5 reanalysis 1958-2022.

to run: python ERA5tailwidth.py [1958] [2023]

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
t1 = int(sys.argv[1])
t2 = int(sys.argv[2])

datadir = '/data/users/samuel.sancar/Heat_extremes/Merged_era5_cubes/'
outdir = '/data/users/samuel.sancar/Heat_extremes/tailwidthslist/'
template = iris.load_cube('template_cube_1deg.nc')

# Set data, variable
reanalysis = 'era5'
variable = 'daily_maximum_2m_temperature'

os.makedirs(outdir, exist_ok=True)

# Calculate tailwidths over years t1-t2
for y in np.arange(t1,t2):

    fn = reanalysis+'_'+variable+'_'+str(y)+'_merged.nc'
    year_cubes = iris.load_cube(datadir+fn)

    #Calculate p99-p87.5 for this year:
    p99 = year_cubes.collapsed('valid_time',iris.analysis.PERCENTILE,percent=99)
    p87p5 =year_cubes.collapsed('valid_time',iris.analysis.PERCENTILE,percent=87.5)
    tail_width = p99-p87p5

    #regrid onto template:
    tail_width_rg = tail_width.regrid(template,iris.analysis.Linear())

    # Save each year's cube separately
    out_file = os.path.join(outdir, f"{reanalysis}_{variable}_{y}_tail_width.nc")
    iris.save(tail_width_rg, out_file)
 