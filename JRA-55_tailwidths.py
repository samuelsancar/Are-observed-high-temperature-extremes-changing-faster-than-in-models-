"""

Script to calculate tailwidths of daily maximum 2m temperature from JRA-55 reanalysis data.

To run: python JRA-55_tailwidths.py

"""

import iris
import numpy as np
import iris.coord_systems as cs
import os

datadir = '/data/users/samuel.sancar/Heat_extremes/JRA-55_tailwidths/daily_maximum__2m_temperature/'
outdir = '/data/users/samuel.sancar/Heat_extremes/JRA-55_tailwidths/Yearly_tailwidths/'
template = iris.load_cube('template_cube_1deg.nc')

reanalysis = 'jra55'
variable = '_minmax_surf.015_tmax.reg_tl319.'

os.makedirs(outdir, exist_ok=True)

# Function to match coordinate systems and regrid
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

# Calculate tailwidths over years t1-t2
for y in np.arange(1958,2023):

    fn = reanalysis+variable+'_'+str(y)+'.nc'
    year_cubes = iris.load_cube(datadir+fn)

    #Calculate p99-p87.5 for this year:
    p99 = year_cubes.collapsed('time',iris.analysis.PERCENTILE,percent=99)
    p87p5 =year_cubes.collapsed('time',iris.analysis.PERCENTILE,percent=87.5)
    tail_width = p99-p87p5

    #regrid onto template:
    tail_width_rg = match_and_regrid(tail_width, template)
    p87p5_regrid = match_and_regrid(p87p5, template)

    out_file = os.path.join(outdir, f"{reanalysis}_{variable}_{y}_tail_width.nc")
    outfile_87p5 = os.path.join(outdir, f"{reanalysis}_{variable}_{y}_87p5.nc")
    iris.save(tail_width_rg, out_file)
    iris.save(p87p5_regrid, outfile_87p5)


