"""
    Code to calculate linear trends in HighResMIP daily tasmax tailwidths. 
    Model member names are imported from the Model_table.xlsx file. 
    The code loads the stacked tail width cubes created by HighResMIP_stack.py,
"""

import iris
import numpy as np
import sys
from iris.cube import CubeList
from iris.analysis import stats
import iris.quickplot as qplt
import os
import sys
import cf_units
from Trend_calculation_functions import (create_time_coord, create_year_cube, maskp87p5, stack_cubes, calculate_trend)
import glob
import pandas as pd

# Beginning and end year
t1 = 1958
t2 = 2023

years = np.arange(int(t1), t2)

df = pd.read_excel("/data/users/samuel.sancar/Model_table.xlsx", sheet_name='Model Stack', engine="openpyxl")

# ensure consistent column names (strip spaces)
df.columns = df.columns.str.strip()

for row in np.arange(0,63):
    
    row_idx = row  
    row = df.iloc[row_idx]

    # Map columns to variables 
    Mode       = str(row["Mode"]).strip()
    model      = str(row["model"]).strip()
    experiment = str(row["experiment"]).strip()
    member     = str(row["member"]).strip()

    print(model, experiment, member)

    Model_dir = f'/data/users/samuel.sancar/Heat_extremes/Tail_widths/'
    outtrend = f'/data/users/samuel.sancar/Heat_extremes/Trends/{Mode}_model_trends/'
    
    data_cube = iris.load_cube(Model_dir+f'{model}_{experiment}_{member}_tailwidths_1958-{t2}.nc')
    p87p5_data_cube = iris.load_cube(Model_dir+f'{model}_{experiment}_{member}_p87p5_1958-{t2}.nc')
    
    year_con = iris.Constraint(sim_year=lambda y: int(t1) <= y <= 2022)
    data_cube = data_cube.extract(year_con)
    p87p5_data_cube = p87p5_data_cube.extract((year_con))
    
    year_cube, p87p5_year_cube = create_year_cube(data_cube, p87p5_data_cube, years)
    trend = calculate_trend(year_cube, data_cube)
    p87p5_trend = calculate_trend(p87p5_year_cube, p87p5_data_cube)
    masked_neg_trends = maskp87p5(p87p5_trend, trend)
        
    masked_neg_trends.rename(f'Tail width trends ({t1}-{t2})')
    p87p5_trend.rename(f'p87p5 trends ({t1}-{t2})')
    
    masked_neg_trends.remove_coord('sim_year')
    p87p5_trend.remove_coord('sim_year')
 
    # Save trend cubes
    out_file_trend = os.path.join(outtrend, f"{model}_{experiment}_{member}_tailwidth_trend_{t1}-{t2}.nc")
    iris.save(masked_neg_trends, out_file_trend)

    p87p5_out_file_trend = os.path.join(outtrend, f"{model}_{experiment}_{member}_p87p5_trend_{t1}-{t2}.nc")
    iris.save(p87p5_trend, p87p5_out_file_trend)
