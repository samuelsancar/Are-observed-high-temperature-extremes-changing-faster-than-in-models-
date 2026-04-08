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
import pandas as pd

def file_covers_year(file_path, year):
    try:
        date_str = file_path.split('_')[-1].replace('.nc', '')
        start_str, end_str = date_str.split('-')
        start_date = datetime.datetime.strptime(start_str, '%Y%m%d')
        end_date = datetime.datetime.strptime(end_str, '%Y%m%d')
        return start_date.year <= year <= end_date.year
    except Exception as e:
        print(f"Failed to parse date from {file_path}: {e}")
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
        return yearly_cube
    
    except Exception as e:
        print(f"Failed to load or process cube for year {y}: {e}")
        return None

def assign_year_time(cube, year):
    try:
        time_unit = cf_units.Unit('days since 1950-01-01', calendar='gregorian') #change to 'days since 1950-01-01'
        time_point = time_unit.date2num(datetime.datetime(year, 7, 1))
        time_coord = iris.coords.DimCoord([time_point], standard_name='time', units=time_unit)

        lat = cube.coord('latitude')
        lon = cube.coord('longitude')

        new_cube = iris.cube.Cube(          # Can maybe get rid and and put inplace cube.rename('Yearly tasmax tailwidth')
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
        print(f"Error assigning time coordinate for year {year}: {e}")
        return None

def process_year(y, model, experiment, member, datadir, future_expt, template_path):
    
    print(f"[...] Starting processing tailwidth for year {y}")
    
    cube = get_yearly_cube(y, model, experiment, member, datadir, future_expt)
    template = iris.load_cube(template_path)
    
    if cube is None:
        return None
    try:
        p99 = cube.collapsed('time', iris.analysis.PERCENTILE, percent=99)
        p87p5 = cube.collapsed('time', iris.analysis.PERCENTILE, percent=87.5)
        tail_width = p99 - p87p5
        
        regridded_tailwidth = tail_width.regrid(template, iris.analysis.Linear())
        regridded_p87p5 = p87p5.regrid(template, iris.analysis.Linear())
        print(f"[✓] Tailwidth calculated and regridded data for year {y}")
        
        return assign_year_time(regridded_tailwidth, y)
    except Exception as e:
        print(f"Error calculating tail width for year {y}: {e}")
        return None
    
def p87p5(y, model, experiment, member, datadir, future_expt, template_path):
    print(f"[...] Starting processing p87p5 for year {y}")
    
    cube = get_yearly_cube(y, model, experiment, member, datadir, future_expt)
    template = iris.load_cube(template_path)
    
    if cube is None:
        return None
    try:
        p87p5 = cube.collapsed('time', iris.analysis.PERCENTILE, percent=87.5)
        regridded_p87p5 = p87p5.regrid(template, iris.analysis.Linear())
        print(f"[✓] p87p5 calculated and regridded data for year {y}")
        
        return assign_year_time(regridded_p87p5, y)
    except Exception as e:
        print(f"Error calculating p87p5 for year {y}: {e}")
        return None

def run_process_year(args):
    return process_year(*args)

def run_process_p87p5(args):
    return p87p5(*args)

def main():
    
    df = pd.read_excel("/data/users/samuel.sancar/Model_table.xlsx", sheet_name='CNRM', engine="openpyxl")

    # ensure consistent column names (strip spaces)
    df.columns = df.columns.str.strip()

    for row in np.arange(11,12): ## Need to also process 0-3
        
        row_idx = row  
        row = df.iloc[row_idx]

        # Map columns to variables 
        Mode       = str(row["Mode"]).strip()
        Group      = str(row["Group"]).strip()
        model      = str(row["model"]).strip()
        experiment = str(row["experiment"]).strip()
        member     = str(row["member"]).strip()
        
        print(f'Calculating yearly tail widths for {model}_{experiment}_{member}')

        datadir = f'/data/users/samuel.sancar/Climate_models/{Mode}/{Group}/'
        outdir = f'/data/users/samuel.sancar/Heat_extremes/Tailwidthslist/{Group}/'
        os.makedirs(outdir, exist_ok=True)

        template_path = 'template_cube_1deg.nc'

        future_expt = {'hist-1950': 'highres-future',
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

        with ThreadPoolExecutor() as executor:
            results_p87p5 = executor.map(run_process_p87p5, args_list)
            for i, results_p87p5 in enumerate(results_p87p5):
                if results_p87p5 is not None:
                    year = years[i]
                    output_path = f'{outdir}{model}_{experiment}_{member}_p87p5_{year}.nc'
                    iris.save(results_p87p5, output_path)
                    print(f"[✓] Saved p87p5 for year {year} to {output_path}")
                else:
                    print(f"[!] Skipped saving for year {years[i]} due to missing data.")

if __name__ == "__main__":
    main()