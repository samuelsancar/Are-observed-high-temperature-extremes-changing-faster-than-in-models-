"""
17/07/25
---------
Code to calculate tail width (p99-p87.5) from HighResMIP model members.

Note_1: The code calculates tail widths for for model members listed in the Model_table.xlsx file, 
which includes both Coupled and SST-based HighResMIP Experiments.These are automatically imported from the Excel file, 
and the code iterates through each row to process the specified model, experiment, and member.

Note_2: The code concatenates/merges tasmax data: hist-1950 (1958-2014) and highres-future (2015-2022) for Coupled experiments, 
and highresSST-present (1958-2014) and highresSST-future (2015-2022) for SST-based experiments.

Note_3: The code can calculate from files covering varying time periods, as long as they include the year being processed. 
It will automatically select the correct file for each year.

Note_4: The code saves yearly tail width and p87.5 cubes separately for each year, model, experiment, and member. 
Directories are created for each institution under the Tailwidthslist directory.

To run: HighResMIP_tailwidths.py
"""

import iris
import numpy as np
import os
import glob
import iris.analysis
from concurrent.futures import ThreadPoolExecutor
from iris.cube import CubeList
import cf_units
import datetime
import iris.coords
import pandas as pd

def file_covers_year(file_path, year):
    try:
        date_str = file_path.split('_')[-1].replace('.nc', '')
        start_str, end_str = date_str.split('-')
        start_date = datetime.datetime.strptime(start_str, '%Y%m%d')
        end_date = datetime.datetime.strptime(end_str, '%Y%m%d')
        return start_date.year <= year <= end_date.year
    except Exception as e:
        print(f"Failed to find date from {file_path}: {e}")
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
        time_unit = cf_units.Unit('days since 1970-01-01', calendar='gregorian') #Check if works, otherwise use 'days since 1950-01-01'
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
        print(f"Error assigning time coordinate for year {year}: {e}")
        return None
    
def process_year(y, model, experiment, member, datadir, future_expt, template_path):
    print(f"[...] Starting processing for year {y}")

    cube = get_yearly_cube(y, model, experiment, member, datadir, future_expt)
    if cube is None:
        return None, None

    template = iris.load_cube(template_path)

    try:
        # Percentiles
        p99 = cube.collapsed('time', iris.analysis.PERCENTILE, percent=99)
        p87p5 = cube.collapsed('time', iris.analysis.PERCENTILE, percent=87.5)

        # Tail width
        tail_width = p99 - p87p5

        # Regrid (once each)
        p87p5 = p87p5.regrid(template, iris.analysis.Linear())
        tail_width = tail_width.regrid(template, iris.analysis.Linear())

        print(f"Percentiles + tail width calculated and regridded for {y}")

        return (
            assign_year_time(tail_width, y),
            assign_year_time(p87p5, y),
        )

    except Exception as e:
        print(f"Error processing year {y}: {e}")
        return None, None

def run_process_year(args):
    return process_year(*args)

def main():
    
    # Beginning and end year+1
    t1 = 1958
    t2 = 2023
    
    df = pd.read_excel("/data/users/samuel.sancar/Model_table.xlsx", sheet_name='Model Stack', engine="openpyxl")

    # ensure consistent column names (strip spaces)
    df.columns = df.columns.str.strip()

    for row in np.arange(0,63): # Range of rows to process (0-62 for all 63 rows)
        
        row_idx = row  
        row = df.iloc[row_idx]

        # Map columns to variables 
        Mode       = str(row["Coupled/SST"]).strip()
        Group      = str(row["Institution"]).strip()
        model      = str(row["Model"]).strip()
        experiment = str(row["Base experiment"]).strip()
        member     = str(row["Member"]).strip()
        
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

        years = list(range(t1, t2))
        args_list = [(y, model, experiment, member, datadir, future_expt, template_path)
                    for y in years]

        with ThreadPoolExecutor() as executor:
            results = executor.map(run_process_year, args_list)

            for i, (tw, p87) in enumerate(results):
                year = years[i]

                if tw is not None:
                    tw_path = f'{outdir}{model}_{experiment}_{member}_tailwidth_{year}.nc'
                    iris.save(tw, tw_path)
                    print(f"Saved tailwidth {year}")

                else:
                    print(f"[!] Tailwidth missing for {year}")

                if p87 is not None:
                    p87_path = f'{outdir}{model}_{experiment}_{member}_p87p5_{year}.nc'
                    iris.save(p87, p87_path)
                    print(f"Saved p87p5 {year}")

                else:
                    print(f"[!] p87p5 missing for {year}")

if __name__ == "__main__":
    main()