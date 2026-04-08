"""
    Stack yearly HighResMIP model tailwidths into one file per model
"""

import iris
import numpy as np
import sys
from iris.cube import CubeList
from iris.analysis import stats
import iris.quickplot as qplt
import os
from iris.coords import AuxCoord
import sys
import cf_units
from Trend_calculation_functions import (create_time_coord, create_year_cube, maskp87p5, stack_cubes, calculate_trend)
import glob
import pandas as pd

t1 = 1958
t2 = 2023
years = np.arange(int(t1), t2)

def add_year_coord(cube, first_year: int):
    n = cube.shape[0]
    sim_years = np.arange(first_year, first_year + n, dtype=int)

    # Attach as an AuxCoord along the time dimension
    time_dim = cube.coord_dims('time')[0]
    sim_year_aux = AuxCoord(sim_years, long_name='sim_year', units='year')
    cube.add_aux_coord(sim_year_aux, data_dims=(time_dim,))
    return cube

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

for row in np.arange(53,63):
    
    row_idx = row  
    print(row)
    row = df.iloc[row_idx]

    # Map columns to variables 
    Mode       = str(row["Mode"]).strip()
    Group      = str(row["Group"]).strip()
    model      = str(row["model"]).strip()
    experiment = str(row["experiment"]).strip()
    member     = str(row["member"]).strip()

    print(Mode, Group, model, experiment, member)

    Model_dir = f'/data/users/samuel.sancar/Heat_extremes/Tailwidthslist/{Group}/'
    outtailwidth = '/data/users/samuel.sancar/Heat_extremes/Tail_widths/'

    time_coord = create_time_coord(years)
    cubes, p87p5_cubes = load_yearly_cubes(Model_dir, model, experiment, member, years)

    cubes = CubeList([iris.util.squeeze(c) for c in cubes])
    p87p5_cubes = CubeList([iris.util.squeeze(c) for c in p87p5_cubes])

    data_cube, p87p5_data_cube = stack_cubes(cubes, p87p5_cubes, time_coord)
    
    cube = add_year_coord(data_cube, t1)
    p87p5 = add_year_coord(p87p5_data_cube, t1)

    cube.long_name = f'Yealy tail widths ({t1}-2022)'
    
    # Save merged tailwidth cube
    out_file_tailwidth = os.path.join(outtailwidth, f"{model}_{experiment}_{member}_tailwidths_1958-2022.nc")
    iris.save(cube, out_file_tailwidth)
    
    print(f"Saved {model}_{experiment}_{member}_tailwidths_1958-2022.nc")
    
    p87p5.long_name = f'Yearly p87p5 ({t1}-2022)'
    
    out_file_p87p5 = os.path.join(outtailwidth, f"{model}_{experiment}_{member}_p87p5_1958-2022.nc")
    iris.save(p87p5, out_file_p87p5)
    
    print(f"Saved {model}_{experiment}_{member}_p87p5_1958-2022.nc")
    
    