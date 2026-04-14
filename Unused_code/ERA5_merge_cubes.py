"""
10/07/25  ***OLD VERSION*** 
---------
Code to merge monthly ERA5 files into yearly files (1958-2022) for comparison with models.

to run: python ERA5_merge_cubes.py [1958] [2023]
"""

import iris
import numpy as np
import sys
import os

#Beginning & end year
t1 = 1958
t2 = 2023

datadir = '/data/users/appldata/Data/OBS-ERA5/daily/2m_temperature/daily_maximum/'
outdir = '/data/users/samuel.sancar/Heat_extremes/Merged_ERA5_cubes/'

# Load template grid for regridding
template = iris.load_cube('template_cube_1deg.nc')

reanalysis = 'era5'
variable = 'daily_maximum_2m_temperature'

os.makedirs(outdir, exist_ok=True)

# Process each day
for y in range(t1, t2):
    file_pattern = f"{datadir}{reanalysis}_{variable}_{y}*.nc"
    cube_list = iris.load(file_pattern)

    iris.util.equalise_attributes(cube_list)
    iris.util.unify_time_units(cube_list)

    merged_cube = cube_list.concatenate_cube()

    # Save each year's cube separately
    out_file = os.path.join(outdir, f"{reanalysis}_{variable}_{y}_merged.nc")
    iris.save(merged_cube, out_file)
    #print(f"Saved: {out_file}")



