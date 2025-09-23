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

def get_yearly_cube(y, model, experiment, member, datadir, future_expt, template_path):
    print(f"[...] Attempting to load data for year {y}")
    expt = experiment if y < 2015 else future_expt
    subdir = f'tasmax_day_{model}_{expt}_{member}'
    search_dir = os.path.join(datadir, subdir)

    template = load_template(template_path)
    if template is None:
        return None

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

        regridded_cube = yearly_cube.regrid(template, iris.analysis.Linear())
        print(f"[✓] Loaded and regridded data for year {y}")
        return regridded_cube
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
    except Exception as e:
        logging.error(f"Error assigning time coordinate for year {year}: {e}")
        logging.error(traceback.format_exc())
        return None

def process_year(y, model, experiment, member, datadir, future_expt, template_path):
    print(f"[...] Starting processing for year {y}")
    cube = get_yearly_cube(y, model, experiment, member, datadir, future_expt, template_path)
    if cube is None:
        return None
    try:
        p99 = cube.collapsed('time', iris.analysis.PERCENTILE, percent=99)
        p87p5 = cube.collapsed('time', iris.analysis.PERCENTILE, percent=87.5)
        tail_width = p99 - p87p5
        return assign_year_time(tail_width, y)
    except Exception as e:
        logging.error(f"Error calculating tail width for year {y}: {e}")
        logging.error(traceback.format_exc())
        return None

def run_process_year(args):
    return process_year(*args)

def main():
    Mode = sys.argv[1]  # 'Coupled' or 'SST'
    Group = sys.argv[2] # See directory names
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
    with ThreadPoolExecutor() as executor:
        results = executor.map(run_process_year, args_list)
        for i, result in enumerate(results):
            if result is not None:
                year = years[i]
                output_path = f'{outdir}{model}_{experiment}_{member}_tailwidth_{year}.nc'
                iris.save(result, output_path)
                print(f"[✓] Saved tail width for year {year} to {output_path}")
            else:
                print(f"[!] Skipped saving for year {years[i]} due to missing data.")


    if len(tail_width_list) != len(years):
        print(f"[!] Warning: Only {len(tail_width_list)} tail width cubes were produced (expected {len(years)}).")
        sys.exit(1)

    time_points = [cube.coord('time').points[0] for cube in tail_width_list]
    sorted_indices = np.argsort(time_points)
    tail_width_list_sorted = CubeList([tail_width_list[i] for i in sorted_indices])

    lat = tail_width_list_sorted[0].coord('latitude')
    lon = tail_width_list_sorted[0].coord('longitude')
    new_time_points = np.array(sorted(time_points), dtype=np.float64)
    standard_unit = cf_units.Unit('days since 1950-01-01', calendar='gregorian')

    for i, cube in enumerate(tail_width_list_sorted):
        # Remove all time coordinates
        for coord in cube.coords():
            if coord.standard_name == 'time':
                cube.remove_coord(coord)

        # Add clean time coordinate
        new_time_coord = iris.coords.DimCoord(
            [new_time_points[i]],
            standard_name='time',
            units=standard_unit,
            var_name='time',
            long_name='time'
        )
        cube.add_dim_coord(new_time_coord, 0)

        # Reassign latitude and longitude
        cube.remove_coord('latitude')
        cube.add_dim_coord(lat, 1)

        cube.remove_coord('longitude')
        cube.add_dim_coord(lon, 2)

    print("\n[🔍] Inspecting time coordinates before merge...\n")
    for i, cube in enumerate(tail_width_list_sorted):
        time_coord = cube.coord('time')
        print(f"Cube {i} - Year: {1958 + i}")
        print(f"  Points: {time_coord.points}")
        print(f"  Units: {time_coord.units}")
        print(f"  Calendar: {time_coord.units.calendar}")
        print(f"  Bounds: {time_coord.bounds}")
        print(f"  Attributes: {time_coord.attributes}")
        print("-" * 40)

    iris.util.equalise_attributes(tail_width_list_sorted)
    yearly_tail_widths = tail_width_list_sorted.merge_cube()
    iris.save(yearly_tail_widths, f'{outdir}{model}_{experiment}_{member}_tailwidths.nc')
    print(f"[✓] Saved tail width cube with shape {yearly_tail_widths.shape}")

if __name__ == "__main__":
    main()

