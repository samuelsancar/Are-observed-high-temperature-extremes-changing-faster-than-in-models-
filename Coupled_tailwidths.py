"""
02/07/25

Calculate tailwidths of models from 1958-2022 for Coupled models (combine historical and future runs)

To run: python Model_trend_plotter.py [model file] [model] [experiment] [member]

"""

import iris
import numpy as np
import sys
import os
import glob
import iris.analysis
from concurrent.futures import ProcessPoolExecutor
from iris.cube import CubeList

def get_yearly_cube(y, model, experiment, member, datadir, future_expt, template):
    try:
        expt = experiment if y < 2015 else future_expt
        subdir = f'tasmax_day_{model}_{expt}_{member}'

        # Define candidate file ranges
        candidate_ranges = [
            (y // 10 * 10, y // 10 * 10 + 9),  # Decadal
            (y // 4 * 4, y // 4 * 4 + 3),      # 4-year
        ]

        # Try each candidate range
        for start, end in candidate_ranges:
            pattern = f'tasmax_day_{model}_{expt}_{member}_*_{start}0101-{end}1231.nc'
            files = glob.glob(os.path.join(datadir, subdir, pattern))
            if files:
                print(f"[✓] Found file for {start}-{end}")
                cube = iris.load_cube(files[0])
                constraint = iris.Constraint(time=lambda cell: cell.point.year == y)
                cube = cube.extract(constraint)
                return cube

        # Check for yearly file
        yearly_pattern = f'tasmax_day_{model}_{expt}_{member}_*_{y}0101-{y}*.nc'
        yearly_files = glob.glob(os.path.join(datadir, subdir, yearly_pattern))
        if yearly_files:
            print(f"[✓] Found existing yearly file for {y}")
            return iris.load_cube(yearly_files[0])

        # Merge monthly files
        monthly_pattern = f'tasmax_day_{model}_{expt}_{member}_*_{y}*.nc'
        monthly_files = glob.glob(os.path.join(datadir, subdir, monthly_pattern))
        if not monthly_files:
            print(f"[!] No monthly files found for year {y}, skipping.")
            return None

        cube_list = iris.load(monthly_files)
        iris.util.equalise_attributes(cube_list)
        iris.util.unify_time_units(cube_list)
        merged_cube = cube_list.concatenate_cube()
        regridded_cube = merged_cube.regrid(template, iris.analysis.Linear())
        print(f"[+] Merged monthly files for {y}")
        return regridded_cube

    except Exception as e:
        print(f"[✗] Error loading or merging year {y}: {e}")
        return None

# Tail width calculation
def process_year(cube, template):
    try:
        if cube is None:
            return None
        p99 = cube.collapsed('time', iris.analysis.PERCENTILE, percent=99)
        p87p5 = cube.collapsed('time', iris.analysis.PERCENTILE, percent=87.5)
        tail_width = p99 - p87p5
        tail_width_rg = tail_width.regrid(template, iris.analysis.Linear())
        return tail_width_rg
    except Exception as e:
        print(f"[✗] Error calculating tail width: {e}")
        return None

def unpack_and_process(args):
    return process_year(*args)

# Main script
if __name__ == "__main__":
    Group = sys.argv[1]
    model = sys.argv[2]
    experiment = sys.argv[3]
    member = sys.argv[4]

    datadir = f'/data/users/samuel.sancar/Climate_models/Coupled/{Group}/'
    outdir = '/data/users/samuel.sancar/Heat_extremes/tailwidthslist/'

    os.makedirs(outdir, exist_ok=True)

    template = iris.load_cube('template_cube_1deg.nc')

    # Determine future experiment
    if experiment == 'hist-1950':
        future_expt = 'highres-future'
    elif experiment == 'highresSST-present':
        future_expt = 'highresSST-future'
    else:
        raise ValueError(f"Unknown experiment: {experiment}")

    # Step 1: Load or merge yearly cubes
    cubes = []
    for y in range(1958, 2023):
        cube = get_yearly_cube(y, model, experiment, member, datadir, future_expt, template)
        if cube:
            cubes.append((cube, template))

    # Step 2: Parallel tail width calculation
    tail_width_list = CubeList()
    with ProcessPoolExecutor() as executor:
        results = executor.map(unpack_and_process, cubes)
        for result in results:
            if result is not None:
                tail_width_list.append(result)

    if len(tail_width_list) == 0:
        print("No tail width cubes were produced.")
        sys.exit(1)

    iris.util.equalise_attributes(tail_width_list)
    yearly_tail_widths = tail_width_list.merge()[0]
    iris.save(yearly_tail_widths, f'{outdir}{model}_{experiment}_{member}_tailwidths.nc')
