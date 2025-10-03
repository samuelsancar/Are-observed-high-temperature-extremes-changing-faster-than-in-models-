"""
22/09/2025
-----------

Code to merge and culculate yearly tail widths for HadGEM3-GC31-LL CANARI large ensemble.

To Run:

"""

import iris
import numpy as np
import sys
import os
import glob
import iris.analysis
from concurrent.futures import ThreadPoolExecutor, as_completed
from iris.cube import CubeList
import cf_units
import datetime
import iris.coords
import traceback
import logging

# Setup logging
logging.basicConfig(filename='tailwidths_processing.log', level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def load_template(template_path):
    try:
        return iris.load_cube(template_path)
    except Exception as e:
        logging.error(f"Failed to load template cube from {template_path}: {e}")
        logging.error(traceback.format_exc())
        return None

def file_covers_year(file_path, year):
    try:
        basename = os.path.basename(file_path)
        date_str = basename.split('_')[1]  #  ie '19501959'
        start_year = int(date_str[:4])
        end_year = int(date_str[4:])
        return start_year <= year <= end_year
    except Exception as e:
        logging.warning(f"Failed to parse year range from {file_path}: {e}")
        return False


def get_yearly_cube(y, member, datadir, template):
    search_dir = datadir
    all_files = glob.glob(os.path.join(search_dir, '*.nc'))
    matching_file = next((f for f in all_files if file_covers_year(f, y)), None)
    logging.info(f"Using file {matching_file} for year {y}")


    if not matching_file:
        logging.warning(f"No file found for year {y}, skipping.")
        return None

    try:
        cube_list = iris.load(matching_file)
        for c in cube_list:
            logging.info(f"Cube time units: {c.coord('time').units}")
        iris.util.equalise_attributes(cube_list)
        iris.util.unify_time_units(cube_list)
        cube = cube_list.concatenate_cube()
        constraint = iris.Constraint(time=lambda cell: cell.point.year == y)
        

        yearly_cube = cube.extract(constraint)

        # Match coordinate systems from model cube to match template
        yearly_cube.coord('latitude').coord_system = template.coord('latitude').coord_system
        yearly_cube.coord('longitude').coord_system = template.coord('longitude').coord_system

        if yearly_cube is None or yearly_cube.shape[0] == 0:
            logging.warning(f"No data for year {y} after extraction.")
            return None

        return yearly_cube.regrid(template, iris.analysis.Linear())
    except Exception as e:
        logging.error(f"Failed to load or process cube for year {y}: {e}")
        logging.error(traceback.format_exc())
        return None
    

def assign_year_time(cube, year):
    try:
        time_unit = cf_units.Unit('days since 1950-01-01', calendar='gregorian')
        time_point = time_unit.date2num(datetime.datetime(year, 7, 1))
        time_coord = iris.coords.DimCoord([time_point], standard_name='time', units=time_unit)
        lat = cube.coord('latitude')
        lon = cube.coord('longitude')

        return iris.cube.Cube(
            cube.data[np.newaxis, ...],
            standard_name=cube.standard_name,
            long_name=cube.long_name,
            var_name=cube.var_name,
            units=cube.units,
            attributes=cube.attributes,
            dim_coords_and_dims=[(time_coord, 0), (lat, 1), (lon, 2)],
            aux_coords_and_dims=[],
        )
    except Exception as e:
        logging.error(f"Error assigning time coordinate for year {year}: {e}")
        logging.error(traceback.format_exc())
        return None

def process_year(args):
    y, member, datadir, template = args
    cube = get_yearly_cube(y, member, datadir, template)
    if cube is None:
        return None, y
    try:
        p99 = cube.collapsed('time', iris.analysis.PERCENTILE, percent=99)
        p87p5 = cube.collapsed('time', iris.analysis.PERCENTILE, percent=87.5)
        tail_width = p99 - p87p5
        return assign_year_time(tail_width, y), y
    except Exception as e:
        logging.error(f"Error calculating tail width for year {y}: {e}")
        logging.error(traceback.format_exc())
        return None, y

def main():
    Model = sys.argv[1]
    datadir = f'/data/users/samuel.sancar/Climate_models/MASS/HadGEM3-LL_r{Model}i1p/'
    outdir = f'/data/users/samuel.sancar/Heat_extremes/Tailwidthslist/MASS/'
    os.makedirs(outdir, exist_ok=True)
    

    template = load_template('template_cube_1deg.nc')
    if template is None:
        raise RuntimeError("Template cube could not be loaded.")


    years = list(range(1958, 2023))
    args_list = [(y, Model, datadir, template) for y in years]
    member = f'r0{Model}i1p00000' 
    tail_width_list = []
    with ThreadPoolExecutor() as executor:
        future_to_year = {executor.submit(process_year, args): args[0] for args in args_list}
        for future in as_completed(future_to_year):
            result, year = future.result()
            if result is not None:
                output_path = f'{outdir}LE_{member}_tailwidth_{year}.nc'
                iris.save(result, output_path)
                tail_width_list.append((year, result))
                print(f"Saved tail width for year {year} to {output_path}")
            else:
                print(f"Skipped saving for year {year} due to missing data.")

    if not tail_width_list:
        print("[!] No tail width cubes were produced.")
        sys.exit(1)

    sorted_cubes = [cube for year, cube in sorted(tail_width_list, key=lambda x: x[0])]
    iris.util.equalise_attributes(sorted_cubes)
    merged_cube = CubeList(sorted_cubes).merge_cube()
    iris.save(merged_cube, f'{outdir}LE_{member}_tailwidths.nc')
    print(f"Saved tail width cube with shape {merged_cube.shape}")

if __name__ == "__main__":
    main()