"""
17/09/2025
----------

Code to aggregate JRA-55 3hourly diagnostic maximum temperature for daily maximum 2m temperature data, 
and calculate annual tail width (p99-p87.5) in tasmax.

"""

import iris
import iris.analysis
import iris.coord_categorisation
import iris.coord_systems as cs
import numpy as np
import os

reanalysis = 'jra55'
variable = 'minmax_surf.015_tmax.reg_tl319.'

template = iris.load_cube('template_cube_1deg.nc')
datadir = '/data/users/samuel.sancar/Climate_models/JRA-55/'
outdir = '/data/users/samuel.sancar/Heat_extremes/JRA-55_tailwidths/'

# Beginning & end year+1
t1 = 1958
t2 = 2023

years = np.arange(t1, t2)

# Define function to match coordinate systems 
def match_and_regrid(tailwidth, template):
    def assign_default_coord_system(cube):
        default_cs = cs.GeogCS(6371229.0)
        for coord_name in ['latitude', 'longitude']:
            coord = cube.coord(coord_name)
            if coord.coord_system is None:
                coord.coord_system = default_cs

    assign_default_coord_system(tailwidth)
    assign_default_coord_system(template)

    for coord_name in ['latitude', 'longitude']:
        merged_coord = tailwidth.coord(coord_name)
        template_coord = template.coord(coord_name)
        if merged_coord.units != template_coord.units:
            merged_coord.units = template_coord.units

    return tailwidth.regrid(template, iris.analysis.Linear())

os.makedirs(outdir, exist_ok=True)

# Calculate tailwidths over years t1-t2
for y in years:

    file_pattern = f"{datadir}{variable}{y}*"
    cube_list = iris.load(file_pattern)

    iris.util.equalise_attributes(cube_list)
    iris.util.unify_time_units(cube_list)

    merged_cube = cube_list.concatenate_cube()

    # Reshape to daily blocks (8 time steps per day)
    no_of_days = merged_cube.coord('time').shape[0] // 8
    merged_cube = merged_cube[:no_of_days * 8]

    # Aggregate to daily maximum
    iris.coord_categorisation.add_day_of_year(merged_cube, 'time', name='day')
    daily_cube = merged_cube.aggregated_by('day', iris.analysis.MAX)

    # Calculate p99-p87.5 for this year:
    p99 = daily_cube.collapsed('time',iris.analysis.PERCENTILE,percent=99)
    p87p5 =daily_cube.collapsed('time',iris.analysis.PERCENTILE,percent=87.5)
    tail_width = p99-p87p5

    #regrid onto template:
    tail_width_rg = match_and_regrid(tail_width, template)
    tail_width_rg.rename('Yearly_tasmax_tailwidth')
    p87p5_regrid = match_and_regrid(p87p5, template)
    p87p5_regrid.rename('Tamsax_p87p5')

    # Save each year's cube separately
    out_file = os.path.join(outdir, f"{reanalysis}_{variable}_{y}_tail_width.nc")
    outfile_87p5 = os.path.join(outdir, f"{reanalysis}_{variable}_{y}_87p5.nc")
    iris.save(tail_width_rg, out_file)
    iris.save(p87p5_regrid, outfile_87p5)
    
    print(f'Succesfully calculated tasmax tailwidth for year {y}.')

print("Processing complete.")





 