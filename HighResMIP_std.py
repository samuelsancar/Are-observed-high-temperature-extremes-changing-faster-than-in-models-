"""
1/04/26
--------

Code to calculate standard deviation of HighResMIP member tail width trends. 

To run: python HighResMIP_std.py
"""

import iris
import numpy as np
import os
import pandas as pd

# Beginning and end year+1
t1 = 1958
t2 = 2023

years = np.arange(t1, t2)

df = pd.read_excel("/data/users/samuel.sancar/Model_table.xlsx", sheet_name='Model Stack', engine="openpyxl")

# ensure consistent column names (strip spaces)
df.columns = df.columns.str.strip()

all_trends = []

for row in np.arange(0,63):
    
    row_idx = row  
    row = df.iloc[row_idx]

    # Map columns to variables 
    Mode       = str(row["Coupled/SST"]).strip()
    model      = str(row["Model"]).strip()
    experiment = str(row["Base experiment"]).strip()
    member     = str(row["Member"]).strip()

    print(model, experiment, member)

    Model_dir = f'/data/users/samuel.sancar/Heat_extremes/Trends/{Mode}_model_trends/'
    outdir = f'/data/users/samuel.sancar/Heat_extremes/Trends/'
    
    cube = iris.load_cube(Model_dir+f'{model}_{experiment}_{member}_tailwidth_trend_1958-{t2-1}.nc') 
    
    cube = iris.util.new_axis(cube)

    model_coord = iris.coords.DimCoord([row_idx],long_name='model',units='1')

    cube.add_dim_coord(model_coord, 0)
    
    all_trends.append(cube)
    
iris.util.unify_time_units(all_trends)
iris.util.equalise_attributes(all_trends)

all_trend_cubelist =iris.cube.CubeList(all_trends)

mm_trends = all_trend_cubelist.concatenate_cube()

mm_std = mm_trends.collapsed('model', iris.analysis.STD_DEV) 

# Save merged tailwidth cube
out_file = os.path.join(outdir, f"HighResMIP_std_{t1}-{t2-1}.nc")
iris.save(mm_std, out_file)

print(f"Saved MM_std_{t1}-{t2-1}.nc")
