"""
    Code to merge yearly netCDF files for models that aren't already so.
    
    To run: python merge_yearly_tailwidths.py 

"""

import iris
import glob
import os
import sys
import numpy as np
from iris.cube import CubeList

Model_file = sys.argv[1]  # e.g., "CNRM-CM6-1-HR_hist-1950"
Model = sys.argv[2]  # e.g., "CNRM-CM6-1-HR"
ensemble_member = sys.argv[3]  # e.g., "r1i1p1f2"

data_dir = f'/data/users/samuel.sancar/Heat_extremes/Tailwidthslist/{Model_file}/'
out_dir = f'/data/users/samuel.sancar/Heat_extremes/Tail_widths/'

file_pattern = os.path.join(data_dir, f'{Model}_{ensemble_member}_tailwidth_*.nc')
print("Looking for files with pattern:", file_pattern)

yearly_files = sorted(glob.glob(file_pattern))
print("Matched files:", yearly_files)

if not yearly_files:
    print("No files matched. Check your pattern and directory.")
    sys.exit(1)

cubes = [iris.load_cube(f) for f in yearly_files]

from iris.cube import CubeList

# Load all cubes
cubes = []
for f in yearly_files:
    cube = iris.load_cube(f)
    cube.rename("tailwidth")
    cube.var_name = "tailwidth"
    cube.standard_name = None  # Remove if inconsistent
    cube.long_name = "Daily Maximum Near-Surface Air Temperature"
    cube.units = "K"  # Set explicitly if needed

    # Remove attributes that may differ
    cube.attributes = {}
    cube.cell_methods = ()
    cube.coord('time').attributes = {}

    cubes.append(cube)

# Concatenate
merged_cube = CubeList(cubes).concatenate_cube()

output_file = os.path.join(out_dir, f'{Model}_{ensemble_member}_tailwidths.nc')
iris.save(merged_cube, output_file)
print(f"Merged NetCDF file saved to {output_file}")


