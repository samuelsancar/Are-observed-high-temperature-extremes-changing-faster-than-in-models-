"""
03/11/2025
----------

Code to calculate trends in tail widths and 87.5th percentiles for CANARI model data.

To run: python CANARI_mean_trend_calc.py
"""

import iris
import numpy as np
import os
from Trend_calculation_functions import (create_time_coord, create_year_cube, maskp87p5, stack_cubes, calculate_trend)
from ERA5_trend_calc import maskp87p5

# Beginning & end year+1
t1 = 1958
t2 = 2023

years = np.arange(t1, t2)

outtrend = '/data/users/samuel.sancar/Heat_extremes/Trends/CANARI_trends/'
    
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

all_tw = []
all_p87p5 = []

for member in range(1, 41):

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

    data_cube = iris.util.new_axis(data_cube)
    p87p5_cube = iris.util.new_axis(p87p5_cube)

    model_coord = iris.coords.DimCoord([member],long_name='model',units='1')

    data_cube.add_dim_coord(model_coord, 0)
    p87p5_cube.add_dim_coord(model_coord.copy(), 0)
    
    all_tw.append(data_cube)
    all_p87p5.append(p87p5_cube)

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

p87p5_data_cube = iris.util.squeeze(p87p5_mmm)

# Calculate trend
year_cube, p87p5_year_cube = create_year_cube(tw_mmm, p87p5_mmm, years)
trend = calculate_trend(year_cube, data_cube)
p87p5_trend = calculate_trend(p87p5_year_cube, p87p5_cube)
masked_trend = maskp87p5(p87p5_trend, trend)

# Long name
masked_trend.rename(f'Ensemble mean tail width trends ({t1}-{t2-1})')
p87p5_trend.rename(f'Ensemble mean p87p5 trends ({t1}-{t2-1})') 

out_file_trend = os.path.join(outtrend, f"CANARI_em_tailwidth_trend_{t1}-{t2-1}.nc")
iris.save(masked_trend, out_file_trend)

p87p5_out_file_trend = os.path.join(outtrend, f"CANARI_em_p87p5_trend_{t1}-{t2-1}.nc")
iris.save(p87p5_trend, p87p5_out_file_trend) # Save to plot stippling
    
print(f'Succesfully processed')