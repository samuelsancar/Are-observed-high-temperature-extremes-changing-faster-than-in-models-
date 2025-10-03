"""
17/07/25
---------

Code to concatenate SST 'daily tasmax tailwidth' into yearly files and calculate tailwidths for years 1958-2023.

To run: python SST_tailwidths.py [model] [experiment] [member]

"""

import iris
import numpy as np
import sys
import os
import glob
import iris.analysis
from concurrent.futures import ThreadPoolExecutor
from iris.cube import CubeList
import cf_units
import datetime
import iris.coords
import traceback
import logging
from iris.exceptions import MergeError

# Setup logging
logging.basicConfig(filename='tailwidths_processing.log', level=logging.INFO,
                    format='%(asctime)s - %(levelname)s - %(message)s')

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

def get_yearly_cube(y, model, experiment, member, datadir, future_expt):
    print(f"[...] Attempting to load data for year {y}")
    expt = experiment if y < 2015 else future_expt
    subdir = f'tasmax_day_{model}_{expt}_{member}'
    search_dir = os.path.join(datadir, subdir)

    all_files = glob.glob(os.path.join(search_dir, '*.nc'))
    matching_file = next((f for f in all_files if file_covers_year(f, y)), None)

    if not matching_file:
        print(f"[!] No file found for year {y}, skipping.")
        return None

    try:
        cube_list = iris.load(matching_file)
        iris.util.equalise_attributes(cube_list)
        iris.util.unify_time_units(cube_list)
        cube = cube_list.concatenate_cube()

        constraint = iris.Constraint(time=lambda cell: cell.point.year == y)
        yearly_cube = cube.extract(constraint)

        if yearly_cube is None or yearly_cube.shape[0] == 0:
            print(f"[!] No data for year {y} after extraction.")
            return None

        print(f"[✓] Loaded data for year {y}")
        return yearly_cube
    except Exception as e:
        logging.error(f"Failed to load or process cube for year {y}: {e}")
        logging.error(traceback.format_exc())
        return None

def assign_uniform_time(cube, year):
    time_unit = cf_units.Unit('days since 1950-01-01', calendar='gregorian')
    time_point = time_unit.date2num(datetime.datetime(year, 7, 1))
    time_coord = iris.coords.DimCoord([time_point], standard_name='time', units=time_unit)

    lat = cube.coord('latitude')
    lon = cube.coord('longitude')

    new_cube = iris.cube.Cube(
        cube.data[np.newaxis, ...],
        standard_name=cube.standard_name,
        long_name=cube.long_name,
        var_name=cube.var_name,
        units=cube.units,
        attributes=cube.attributes,
        dim_coords_and_dims=[(time_coord, 0), (lat, 1), (lon, 2)],
        aux_coords_and_dims=[],
    )

    return new_cube

def process_year(y, model, experiment, member, datadir, future_expt, template_path):
    print(f"[...] Starting processing for year {y}")
    cube = get_yearly_cube(y, model, experiment, member, datadir, future_expt)
    if cube is None:
        return None
    try:
        p99 = cube.collapsed('time', iris.analysis.PERCENTILE, percent=99)
        p87p5 = cube.collapsed('time', iris.analysis.PERCENTILE, percent=87.5)
        tail_width = p99 - p87p5

        template = load_template(template_path)
        if template is None:
            return None

        regridded_tail_width = tail_width.regrid(template, iris.analysis.Linear())
        regridded_p87p5 = p87p5.regrid(template, iris.analysis.Linear())
        tailwidth_cube = assign_uniform_time(regridded_tail_width, y)
        p87p5_cube = assign_uniform_time(regridded_p87p5, y)
        return tailwidth_cube, p87p5_cube
    except Exception as e:
        logging.error(f"Error calculating tail width for year {y}: {e}")
        logging.error(traceback.format_exc())
        return None

def run_process_year(args):
    return process_year(*args)

def validate_and_align_time(cube, reference_unit):
    try:
        time_coords = cube.coords('time')
        if not time_coords:
            logging.warning("Cube missing time coordinate. Assigning default time.")
            time_point = reference_unit.date2num(datetime.datetime(2000, 1, 1))
            time_coord = iris.coords.DimCoord([time_point], standard_name='time', units=reference_unit)
            cube.add_dim_coord(time_coord, 0)
        else:
            time_coord = time_coords[0]
            if time_coord.units != reference_unit:
                logging.info(f"Converting time units from {time_coord.units} to {reference_unit}")
                time_coord.convert_units(reference_unit)
            if time_coord.ndim != 1 or time_coord.shape[0] != cube.shape[0]:
                logging.warning("Malformed time coordinate. Replacing.")
                cube.remove_coord('time')
                time_point = reference_unit.date2num(datetime.datetime(2000, 1, 1))
                time_coord = iris.coords.DimCoord([time_point], standard_name='time', units=reference_unit)
                cube.add_dim_coord(time_coord, 0)
            elif time_coord in cube.aux_coords:
                iris.util.promote_aux_coord_to_dim_coord(cube, time_coord)
        return cube
    except Exception as e:
        logging.error(f"Time coordinate alignment failed: {e}")
        return None

def clean_and_merge_cubes(cube_list):
    valid_cubes = CubeList()
    time_points = []

    reference_unit = cf_units.Unit('days since 1950-01-01', calendar='gregorian')

    for cube in cube_list:
        fixed_cube = validate_and_align_time(cube, reference_unit)
        if fixed_cube:
            valid_cubes.append(fixed_cube)
            time_points.append(fixed_cube.coord('time').points[0])
        else:
            logging.warning("Skipping cube due to unresolved time issues.")

    if not valid_cubes:
        raise ValueError("No valid cubes to merge.")

    iris.util.equalise_attributes(valid_cubes)
    iris.util.unify_time_units(valid_cubes)

    sorted_indices = np.argsort(time_points)
    sorted_cubes = CubeList([valid_cubes[i] for i in sorted_indices])

    try:
        merged_cube = sorted_cubes.merge_cube()
        logging.info(f"Successfully merged {len(sorted_cubes)} cubes.")
        return merged_cube
    except MergeError as e:
        logging.error(f"MergeError: {e}")
        for i, cube in enumerate(sorted_cubes):
            logging.error(f"Cube {i} time coord: {cube.coord('time')}")
        raise

def main():
    Mode = sys.argv[1]
    Group = sys.argv[2]
    model = sys.argv[3]
    experiment = sys.argv[4]
    member = sys.argv[5]
    datadir = f'/data/users/samuel.sancar/Climate_models/{Mode}/{Group}/'
    outdir = f'/data/users/samuel.sancar/Heat_extremes/Tailwidthslist/{Group}/'
    os.makedirs(outdir, exist_ok=True)

    template_path = 'template_cube_1deg.nc'

    future_expt = {
        'hist-1950': 'highres-future',
        'highresSST-present': 'highresSST-future'
    }.get(experiment)

    if future_expt is None:
        raise ValueError(f"Unknown experiment: {experiment}")

    years = list(range(1958, 2023))
    args_list = [(y, model, experiment, member, datadir, future_expt, template_path) for y in years]

    tail_width_list = CubeList()
    lwperc_list = CubeList()

    with ThreadPoolExecutor() as executor:
        results = executor.map(run_process_year, args_list)
        for i, result in enumerate(results):
            if result is not None:
                tailwidth_cube, p87p5_cube = result
                year = years[i]

                if tailwidth_cube is not None:
                    iris.save(tailwidth_cube, f'{outdir}{model}_{experiment}_{member}_tailwidth_{year}.nc')
                    tail_width_list.append(tailwidth_cube)

                if p87p5_cube is not None:
                    iris.save(p87p5_cube, f'{outdir}{model}_{experiment}_{member}_p87p5_{year}.nc')
                    lwperc_list.append(p87p5_cube)

    for cube_list, var_name in [(tail_width_list, 'tailwidths'), (lwperc_list, '87p5')]:
        if len(cube_list) != len(years):
            print(f"[!] Warning: Only {len(cube_list)} {var_name} cubes were produced (expected {len(years)}).")
            sys.exit(1)

        try:
            merged_cube = clean_and_merge_cubes(cube_list)
            iris.save(merged_cube, f'{outdir}{model}_{experiment}_{member}_{var_name}.nc')
            print(f"[✓] Saved {var_name} cube with shape {merged_cube.shape}")
        except Exception as e:
            print(f"[!] Failed to merge {var_name} cubes: {e}")

if __name__ == "__main__":
    main()