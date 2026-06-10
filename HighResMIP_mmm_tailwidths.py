"""
03/04/2026
----------

Code to calculate the HighResMIP interannual multi-model mean tail widths.

To run: HighResMIP_mmm_tailwidths.py
"""

import iris
import numpy as np
import sys
from iris.cube import CubeList
import os
import sys
from Trend_calculation_functions import (create_time_coord, stack_cubes)
import pandas as pd

# Beginning and end years
t1 = 1958
t2 = 2023

years = np.arange(t1, t2)

def load_yearly_cubes(datadir, model, expt, member, years):
    
    cubes = iris.cube.CubeList()
    for year in years:
        filename = f"{model}_{expt}_{member}_tailwidth_{year}.nc"
        cube = iris.load_cube(os.path.join(datadir, filename))
        cubes.append(cube)
        
    p87p5_cubes = iris.cube.CubeList()
    for year in years:
        p87p5_filename = f"{model}_{expt}_{member}_p87p5_{year}.nc"
        p87p5_cube = iris.load_cube(os.path.join(datadir, p87p5_filename))
        p87p5_cubes.append(p87p5_cube)
        
    return cubes, p87p5_cubes

df = pd.read_excel("/data/users/samuel.sancar/Model_table.xlsx", sheet_name=0, engine="openpyxl")

# ensure consistent column names (strip spaces)
df.columns = df.columns.str.strip()

all_tw = []
all_p87p5 = []

for row in np.arange(0, 63):
    
    row_idx = row  
    print(row)
    row = df.iloc[row_idx]

    # Map columns to variables 
    Mode       = str(row["Coupled/SST"]).strip()
    Group      = str(row["Institution"]).strip()
    model      = str(row["Model"]).strip()
    experiment = str(row["Base experiment"]).strip()
    member     = str(row["Member"]).strip()

    print(Mode, Group, model, experiment, member)

    Model_dir = f'/data/users/samuel.sancar/Heat_extremes/Tailwidthslist/{Group}/'
    outtailwidth = '/data/users/samuel.sancar/Heat_extremes/Tail_widths/'

    time_coord = create_time_coord(years)
    cubes, p87p5_cubes = load_yearly_cubes(Model_dir, model, experiment, member, years)

    cubes = CubeList([iris.util.squeeze(c) for c in cubes])
    p87p5_cubes = CubeList([iris.util.squeeze(c) for c in p87p5_cubes])
    
    cube, p87p5 = stack_cubes(cubes, p87p5_cubes, time_coord)
    
    cube.long_name = f'Yearly tail widths ({t1}-{t2-1})'
    p87p5.long_name = f'Yearly p87p5 ({t1}-{t2-1})'

    cube = iris.util.new_axis(cube)
    p87p5 = iris.util.new_axis(p87p5)

    model_coord = iris.coords.DimCoord([row_idx],long_name='model',units='1')

    cube.add_dim_coord(model_coord, 0)
    p87p5.add_dim_coord(model_coord.copy(), 0)
    
    all_tw.append(cube)
    all_p87p5.append(p87p5)
    
iris.util.unify_time_units(all_tw)
iris.util.unify_time_units(all_p87p5)

iris.util.equalise_attributes(all_tw)
iris.util.equalise_attributes(all_p87p5)

tw_cube_list = iris.cube.CubeList(all_tw)
p87p5_cube_list = iris.cube.CubeList(all_p87p5)

tw_mm = tw_cube_list.concatenate_cube()
p87p5_mm = p87p5_cube_list.concatenate_cube()

tw_mmm = tw_mm.collapsed('model',iris.analysis.MEAN)
tw_std = tw_mm.collapsed('model', iris.analysis.STD_DEV)
p87p5_mmm = p87p5_mm.collapsed('model', iris.analysis.MEAN)

tw_mmm.long_name = 'Multi-model mean yearly tail width'
tw_std.long_name = 'Inter-model standard deviation of yearly tail width'
p87p5_mmm.long_name = 'Multi-model mean yearly tail width'
    
# Save merged tailwidth cube
out_file_tailwidth = os.path.join(outtailwidth, f"HighResMIP_MMM_tailwidths_{t1}-{t2-1}.nc")
iris.save(tw_mmm, out_file_tailwidth)

print(f"Saved MMM_tailwidths_{t1}-{t2-1}.nc")
    
# Save merged p87p5 cube
out_file_p87p5 = os.path.join(outtailwidth, f"HighResMIP_MMM_p87p5_{t1}-{t2-1}.nc")
iris.save(p87p5, out_file_p87p5)
    
print(f"Saved MMM_p87p5_{t1}-{t2-1}.nc")