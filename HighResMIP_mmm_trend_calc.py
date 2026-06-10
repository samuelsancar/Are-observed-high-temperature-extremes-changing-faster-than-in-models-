"""
06/04/2026
----------

Code to calculate linear trends in HighResMIP interannual mutli-model mean tailwidths. 

To run: python HighResMIP_mmm_trend_calc.py
"""

import iris
import numpy as np
import os
from Trend_calculation_functions import (create_year_cube, maskp87p5, calculate_trend) 
import pandas as pd

# Beginning and end year+1
t1 = 1958
t2 = 2023

years = np.arange(t1, t2)

Model_dir = f'/data/users/samuel.sancar/Heat_extremes/Tail_widths/'
outtrend = f'/data/users/samuel.sancar/Heat_extremes/Trends/'
    
data_cube = iris.load_cube(Model_dir+f'HighResMIP_MMM_tailwidths_1958-{t2-1}.nc') 
p87p5_data_cube = iris.load_cube(Model_dir+f'HighResMIP_MMM_p87p5_1958-{t2-1}.nc') 

cube_start_year = 1958 # Define years that tail width file covers
cube_end_year   = 2022

full_years = np.arange(cube_start_year, cube_end_year + 1)

mask = (full_years >= t1) & (full_years <= (t2 - 1)) # Define mask for years to calc trend over (ie, if t1 = 1980)

data_cube = data_cube[mask, :, :] # Mask years for trend calc.
p87p5_data_cube = iris.util.squeeze(p87p5_data_cube)
p87p5_data_cube = p87p5_data_cube[mask, :, :, :]
years = full_years[mask]

year_cube, p87p5_year_cube = create_year_cube(data_cube, p87p5_data_cube, years)
trend = calculate_trend(year_cube, data_cube)
p87p5_trend = calculate_trend(p87p5_year_cube, p87p5_data_cube)
masked_neg_trends = maskp87p5(p87p5_trend, trend)
        
masked_neg_trends.rename(f'Tail width trends ({t1}-{t2-1})')
p87p5_trend.rename(f'p87p5 trends ({t1}-{t2-1})')
 
# Save trend cubes
out_file_trend = os.path.join(outtrend, f"HighResMIP_MMM_tailwidth_trend_{t1}-{t2-1}.nc")
iris.save(masked_neg_trends, out_file_trend)

p87p5_out_file_trend = os.path.join(outtrend, f"HighResMIP_MMM_p87p5_trend_{t1}-{t2-1}.nc")
iris.save(p87p5_trend, p87p5_out_file_trend) # Save to plot stippling