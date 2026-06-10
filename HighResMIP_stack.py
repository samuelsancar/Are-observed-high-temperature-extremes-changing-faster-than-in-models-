"""
29/07/2026
----------

Code to stack yearly HighResMIP tailwidths into one file per model member. 
Names of Model members are imported from the Model_table.xlsx file.
Stacked cubes are saved in the Tail_widths directory, 
with one file per model member containing the tail widths for each year from 1958-2022.

To run: python HighResMIP_stack.py
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
    
    # Save merged tailwidth cube
    out_file_tailwidth = os.path.join(outtailwidth, f"{model}_{experiment}_{member}_tailwidths_{t1}-{t2-1}.nc")
    iris.save(cube, out_file_tailwidth)
    
    print(f"Saved {model}_{experiment}_{member}_tailwidths_{t1}-{t2-1}.nc")
    
    # Save merged p87p5 cube
    out_file_p87p5 = os.path.join(outtailwidth, f"{model}_{experiment}_{member}_p87p5_{t1}-{t2-1}.nc")
    iris.save(p87p5, out_file_p87p5)
    
    print(f"Saved {model}_{experiment}_{member}_p87p5_{t1}-{t2-1}.nc")
    
    