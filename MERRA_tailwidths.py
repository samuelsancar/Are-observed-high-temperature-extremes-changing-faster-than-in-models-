"""
Code to calculate MERRA tailwidths, from 1980-2022.
"""

import iris
import iris.util
import os
from Lat_Land_mask_function import lat_land_mask
import numpy as np
import iris.analysis.stats
from iris.coords import DimCoord

years = np.arange(2020,2021)

datadir = '/data/users/samuel.sancar/Climate_models/MERRA/'
outdir = '/data/users/samuel.sancar/Heat_extremes/MERRA_tailwidths/'
template = iris.load_cube('template_cube_1deg.nc')

os.makedirs(outdir, exist_ok=True)

for y in years:
    
    print(f'Calculating tasmax tailwidth for year {y}...')
    
    constraint = iris.Constraint(cube_func =lambda c: c.var_name == 'T2MMAX' or c.long_name == '2-meter_air_temperature')
    file_pattern = f"{datadir}MERRA2_*.statD_2d_slv_Nx.{y}*.nc"
    cube_list = iris.load(file_pattern, constraint)

    iris.util.equalise_attributes(cube_list)
    iris.util.unify_time_units(cube_list)
    
    for c in cube_list:
        c.var_name == 'tasmax'
        c.rename('air_temperature')
        c.attributes = {}
        if str(c.units) == 'unknown':
            c.units = 'K'
                
        c.attributes = {}

        for coord in c.coords():
            coord.attributes = {}

        c.cell_methods = ()

    # Merge monthly files into yearly cube
    merged_cube = cube_list.concatenate_cube()

    #Calculate p99-p87.5 for this year:
    p99 = merged_cube.collapsed('time',iris.analysis.PERCENTILE,percent=99)
    p87p5 = merged_cube.collapsed('time',iris.analysis.PERCENTILE,percent=87.5)

    tail_width = p99-p87p5

    #regrid 
    tail_width_rg = tail_width.regrid(template,iris.analysis.Linear())
    tail_width_rg.rename('Yearly_tasmax_tailwidth')
    p87p5_regrid = p87p5.regrid(template,iris.analysis.Linear())
    p87p5_regrid.rename('Tamsax_p87p5')

    # Save each year's cube separately
    out_file = os.path.join(outdir, f"MERRA_daily_maximum_2m_temperature_{y}_tail_width.nc")
    outfile_87p5 = os.path.join(outdir, f"daily_maximum_2m_temperature_{y}_87p5.nc")
    iris.save(tail_width_rg, out_file)
    iris.save(p87p5_regrid, outfile_87p5)
    
    print(f'Succesfully calculated tasmax tailwidth for year {y}.')
    