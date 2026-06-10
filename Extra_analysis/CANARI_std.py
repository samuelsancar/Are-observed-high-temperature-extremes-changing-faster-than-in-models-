'''
15/05/26
--------

Code to calculate CANARI std in model tail width trends

To run: python CANARI_std.py
'''

import iris
import numpy as np
import os

t1 = 1958 
t2 = 2023

all_trends = []

for member in range(1, 41):

    print(f'Processing member {member}')
    filename = f"{member}_tailwidth_trend_{t1}-{t2-1}.nc"
    datadir = f'/data/users/samuel.sancar/Heat_extremes/Trends/CANARI_trends/'
    data_cube = iris.load_cube(datadir + filename)

    data_cube = iris.util.new_axis(data_cube)

    model_coord = iris.coords.DimCoord([member],long_name='model',units='1')

    data_cube.add_dim_coord(model_coord, 0)
    
    all_trends.append(data_cube)

iris.util.unify_time_units(all_trends)

iris.util.equalise_attributes(all_trends)

trend_cube_list = iris.cube.CubeList(all_trends)

mm_trend = trend_cube_list.concatenate_cube()

mm_std = mm_trend.collapsed('model',iris.analysis.STD_DEV)

mm_std.long_name = 'Ensemble std'

outtrend = '/data/users/samuel.sancar/Heat_extremes/Trends/CANARI_trends/'
out_file_trend = os.path.join(outtrend, f"CANARI_std_{t1}-{t2-1}.nc")
iris.save(mm_std, out_file_trend)

