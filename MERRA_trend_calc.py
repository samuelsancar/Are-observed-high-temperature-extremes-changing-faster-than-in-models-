"""
    Code to calculate linear trends in ERA5 daily tasmax tailwidths.
"""

import iris
import os
import numpy as np
from Trend_calculation_functions import (load_yearly_cubes, create_time_coord, create_year_cube, maskp87p5, stack_cubes, calculate_trend)

# Beginning and end year+1
t1 = 1980 
t2 = 2023 

years = np.arange(t1, t2)
    
reanalysis = 'MERRA'
variable = 'daily_maximum_2m_temperature'
    
datadir = '/data/users/samuel.sancar/Heat_extremes/MERRA_tailwidths/'
outdir = '/data/users/samuel.sancar/Heat_extremes/Trends/Reanalysis_trends/'
p87p5 = '/data/users/samuel.sancar/Heat_extremes/MERRA_tailwidths/'
    
os.makedirs(outdir, exist_ok=True)

time_coord = create_time_coord(years)
cubes, p87p5_cubes = load_yearly_cubes(datadir, p87p5 , reanalysis, variable, years)
data_cube, p87p5_data_cube = stack_cubes(cubes, p87p5_cubes, time_coord)
    
# Ensure time coordinate is attached
if not data_cube.coords('time'):
    data_cube.add_dim_coord(time_coord, 0)
if not p87p5_data_cube.coords('time'):
    p87p5_data_cube.add_dim_coord(time_coord, 0)
        
year_cube, p87p5_year_cube = create_year_cube(data_cube, p87p5_data_cube, years)
trend = calculate_trend(year_cube, data_cube)
p87p5_trend = calculate_trend(p87p5_year_cube, p87p5_data_cube)
masked_neg_trends = maskp87p5(p87p5_trend, trend)
    
masked_neg_trends.rename(f'Tail width trends ({t1}-{t2-1})')
p87p5_trend.rename(f'p87p5 trends ({t1}-{t2-1})')

# Save the trend cubes 
out_file = os.path.join(outdir, f"{reanalysis}_{variable}_tailwidth_trends_{t1}-{t2-1}.nc")
iris.save(masked_neg_trends, out_file)
    
p87p5_out_file = os.path.join(outdir, f"{reanalysis}_{variable}_p87p5_trends_{t1}-{t2-1}.nc")
iris.save(p87p5_trend, p87p5_out_file) # Save to plot stippling
