"""
17/07/25
---------
Code to calculate linear trends in ERA5 daily tasmax tailwidths.

"""

import os
import numpy as np
import iris
from iris.analysis import stats
import iris.util
from iris.cube import CubeList
from iris.coords import AuxCoord
import cf_units
from Trend_calculation_functions import (create_time_coord, create_year_cube, maskp87p5, load_yearly_cubes, stack_cubes, calculate_trend)

# Beginning & end year+1
t1 = 1979
t2 = 2023 
    
years = np.arange(t1, t2)
    
reanalysis = 'era5'
variable = 'daily_maximum_2m_temperature'
    
datadir = '/data/users/samuel.sancar/Heat_extremes/ERA5_tailwidths/'
outdir = '/data/users/samuel.sancar/Heat_extremes/Trends/Reanalysis_trends/'
p87p5 = '/data/users/samuel.sancar/Heat_extremes/ERA5_tailwidths/'
    
os.makedirs(outdir, exist_ok=True)

time_coord = create_time_coord(years)
cubes, p87p5_cubes = load_yearly_cubes(datadir, p87p5 , reanalysis, variable, years)
data_cube, p87p5_cube = stack_cubes(cubes, p87p5_cubes, time_coord)
    
# Ensure time coordinate is attached
if not data_cube.coords('time'):
    data_cube.add_dim_coord(time_coord, 0)
if not p87p5_cube.coords('time'):
        p87p5_cube.add_dim_coord(time_coord, 0)
    
data_cube.long_name = f'Yealy tail widths ({t1}-{t2-1}'
p87p5_cube.long_name = f'Yearly p87p5 ({t1}-{t2-1}'

# Save the combined yearly tail width cubes for reference        
iris.save(data_cube, os.path.join(datadir, f"{reanalysis}_{variable}_tailwidths_{t1}-{t2-1}.nc"))




    
    



