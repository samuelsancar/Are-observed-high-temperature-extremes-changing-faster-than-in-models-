"""

Code to download, find daily maximum 2m temperature data from JRA-55 and find annual tail-width.

To run: python JRA55_tasmax.py

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
        date_str = file_path.split('_')[-1].replace('.nc', '')
        start_str, end_str = date_str.split('-')
        start_date = datetime.datetime.strptime(start_str, '%Y%m%d')
        end_date = datetime.datetime.strptime(end_str, '%Y%m%d')
        return start_date.year <= year <= end_date.year
    except Exception as e:
        logging.warning(f"Failed to parse date from {file_path}: {e}")
        return False

def get_yearly_cube(y, model, experiment, member, datadir, future_expt, template):
    expt = experiment if y < 2015 else future_expt
    subdir = f'tasmax_day_{model}_{expt}_{member}'
    search_dir = os.path.join(datadir, subdir)
    all_files = glob.glob(os.path.join(search_dir, '*.nc'))
    matching_file = next((f for f in all_files if file_covers_year(f, y)), None)

    if not matching_file:
        logging.warning(f"No file found for year {y}, skipping.")
        return None

    try:
        cube_list = iris.load(matching_file)
        iris.util.equalise_attributes(cube_list)
        iris.util.unify_time_units(cube_list)
        cube = cube_list.concatenate_cube()
        constraint = iris.Constraint(time=lambda cell: cell.point.year == y)
        yearly_cube = cube.extract(constraint)

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
    y, model, experiment, member, datadir, future_expt, template = args
    cube = get_yearly_cube(y, model, experiment, member, datadir, future_expt, template)
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
    Mode, Group, model, experiment, member = sys.argv[1:6]
    datadir = f'/data/users/samuel.sancar/Climate_models/{Mode}/{Group}/'
    outdir = f'/data/users/samuel.sancar/Heat_extremes/tailwidthslist/{Group}/'
    os.makedirs(outdir, exist_ok=True)

    template = load_template('template_cube_1deg.nc')
    if template is None:
        raise RuntimeError("Template cube could not be loaded.")

    future_expt = {
        'hist-1950': 'highres-future',
        'highresSST-present': 'highresSST-future'
    }.get(experiment)

    if future_expt is None:
        raise ValueError(f"Unknown experiment: {experiment}")

    years = list(range(1958, 2023))
    args_list = [(y, model, experiment, member, datadir, future_expt, template) for y in years]

    tail_width_list = []
    with ThreadPoolExecutor() as executor:
        future_to_year = {executor.submit(process_year, args): args[0] for args in args_list}
        for future in as_completed(future_to_year):
            result, year = future.result()
            if result is not None:
                output_path = f'{outdir}{model}_{experiment}_{member}_tailwidth_{year}.nc'
                iris.save(result, output_path)
                tail_width_list.append((year, result))
                print(f"[✓] Saved tail width for year {year} to {output_path}")
            else:
                print(f"[!] Skipped saving for year {year} due to missing data.")

    if not tail_width_list:
        print("[!] No tail width cubes were produced.")
        sys.exit(1)

    sorted_cubes = [cube for year, cube in sorted(tail_width_list, key=lambda x: x[0])]
    iris.util.equalise_attributes(sorted_cubes)
    merged_cube = CubeList(sorted_cubes).merge_cube()
    iris.save(merged_cube, f'{outdir}{model}_{experiment}_{member}_tailwidths.nc')
    print(f"[✓] Saved tail width cube with shape {merged_cube.shape}")

if __name__ == "__main__":
    main()





 