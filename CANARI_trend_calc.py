"""
03/11/2025
----------

Code to caluclate trends in tail widths and 87.5th percentiles for CANARI model data.
 
To run: python CANARI_trend_calc.py 
"""

import iris
import numpy as np
import os
from Trend_calculation_functions import (create_time_coord, create_year_cube, maskp87p5, stack_cubes, calculate_trend)
from ERA5_trend_calc import maskp87p5

# Beginning & end year+1
t1 = 1980
t2 = 2023

years = np.arange(t1, t2)

outtrend = '/data/users/samuel.sancar/Heat_extremes/Trends/CANARI_trends/'
outtailwidth = '/data/users/samuel.sancar/Heat_extremes/Tail_widths/canari/'
    
def load_yearly_cubes(datadir, year, member):
        
    cubes = iris.cube.CubeList()
    p87p5_cubes = iris.cube.CubeList()
    
    for year in years:
        
        fn = f"tailwidth_member{member}_{year}.nc"
        fn_87p5 = f"p875_member{member}_{year}.nc"
            
        year_cube = iris.load_cube(datadir + fn)
        p87p5_year_cube = iris.load_cube(datadir + fn_87p5)
        
        cubes.append(year_cube)
        p87p5_cubes.append(p87p5_year_cube)
            
    return cubes, p87p5_cubes

for member in np.arange(1, 41):

    print(f'Processing member {member}')
    datadir = f'/data/users/samuel.sancar/Heat_extremes/Tail_widths/canari/member{member}/'
        
    # Stack yearly tailwidths into single file.
    cubes = iris.cube.CubeList()
    p87p5_cubes = iris.cube.CubeList()
    
    for year in years:
        
        fn = f"tailwidth_member{member}_{year}.nc"
        fn_87p5 = f"p875_member{member}_{year}.nc"
            
        year_cube = iris.load_cube(datadir + fn)
        p87p5_year_cube = iris.load_cube(datadir + fn_87p5)
        
        cubes.append(year_cube)
        p87p5_cubes.append(p87p5_year_cube)

    time_coord = create_time_coord(years)
    data_cube, p87p5_cube = stack_cubes(cubes, p87p5_cubes, time_coord)
    
    data_cube.long_name = f'Yealy tail widths ({t1}-{t2-1})'
    p87p5_cube.long_name = f'Yearly p87p5 ({t1}-{t2-1})'

    # Save stacked tailwidth cube for later use
    out_file_tailwidth = os.path.join(outtailwidth, f"{member}_tailwidths_{t1}-{t2-1}.nc")
    iris.save(data_cube, out_file_tailwidth)

    # Calculate trend
    year_cube, p87p5_year_cube = create_year_cube(data_cube, p87p5_cube, years)
    trend = calculate_trend(year_cube, data_cube)
    p87p5_trend = calculate_trend(p87p5_year_cube, p87p5_cube)
    masked_trend = maskp87p5(p87p5_trend, trend)

    # Long name
    masked_trend.rename(f'Tail width trends ({t1}-{t2-1})')
    p87p5_trend.rename(f'p87p5 trends ({t1}-{t2-1})') 

    out_file_trend = os.path.join(outtrend, f"{member}_tailwidth_trend_{t1}-{t2-1}.nc")
    iris.save(masked_trend, out_file_trend)

    p87p5_out_file_trend = os.path.join(outtrend, f"{member}_p87p5_trend_{t1}-{t2-1}.nc")
    iris.save(p87p5_trend, p87p5_out_file_trend) # Save to plot stippling
    
    print(f'Succesfully processed member {member}')