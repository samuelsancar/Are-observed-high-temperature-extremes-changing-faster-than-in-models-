"""
31/10/2025
----------

Code to calculate tailwidths from CANARI data.

"""

import iris
import glob
import os
import re
from tqdm import tqdm
from concurrent.futures import ThreadPoolExecutor, as_completed

# Paths
base_data_path = "/data/users/samuel.sancar/Climate_models/CANARI"
output_base_dir = "/data/users/samuel.sancar/Heat_extremes/Tail_widths/canari"
template_path = "template_cube_1deg.nc"

os.makedirs(output_base_dir, exist_ok=True)
template_cube = iris.load_cube(template_path)

def process_file(file_path):
    try:
        cube = iris.load_cube(file_path)
        time_coord = 'valid_time' if cube.coords('valid_time') else 'time'

        # Compute percentiles
        p99_cube = cube.collapsed(time_coord, iris.analysis.PERCENTILE, percent=99)
        p875_cube = cube.collapsed(time_coord, iris.analysis.PERCENTILE, percent=87.5)

        # Tailwidth
        tailwidth_data = p99_cube.data - p875_cube.data
        tailwidth_cube = p99_cube.copy(data=tailwidth_data)
        tailwidth_cube.rename("tailwidth_of_daily_max_temperature")

        # Extract year from filename
        year_match = re.search(r'(\d{4})', os.path.basename(file_path))
        year = int(year_match.group(1)) if year_match else None
        if year:
            aux_coord = iris.coords.AuxCoord(year, long_name='year', units='year')
            tailwidth_cube.add_aux_coord(aux_coord)
            p875_cube.add_aux_coord(aux_coord)

        # Regrid
        tailwidth_cube = tailwidth_cube.regrid(template_cube, iris.analysis.Linear())
        p875_cube = p875_cube.regrid(template_cube, iris.analysis.Linear())

        return year, tailwidth_cube, p875_cube
    except Exception as e:
        print(f"[!] Failed to process {file_path}: {e}")
        return None

# Loop over members
for member in range(36, 41):
    data_path = f"{base_data_path}/member{member}"
    member_output_dir = os.path.join(output_base_dir, f"member{member}")
    os.makedirs(member_output_dir, exist_ok=True)

    file_list = sorted(glob.glob(f"{data_path}/*.nc"))
    if not file_list:
        print(f"[!] No files found for member {member}. Skipping.")
        continue

    results = {}
    max_workers = min(8, os.cpu_count())
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = {executor.submit(process_file, f): f for f in file_list}

        for future in tqdm(as_completed(futures), total=len(futures),
                           desc=f"Member {member}: Processing files", unit="file"):
            file_path = futures[future]
            result = future.result()
            if result:
                year, tailwidth_cube, p875_cube = result
                results[file_path] = (year, tailwidth_cube, p875_cube)

    # Sort by file name to preserve chronological order
    tailwidth_cubes, p875_cubes = [], []
    for file_path in sorted(results.keys()):
        year, tailwidth_cube, p875_cube = results[file_path]
        tailwidth_cubes.append(tailwidth_cube)
        p875_cubes.append(p875_cube)

        # Save individual yearly cubes
        tailwidth_year_path = os.path.join(member_output_dir, f"tailwidth_member{member}_{year}.nc")
        p875_year_path = os.path.join(member_output_dir, f"p875_member{member}_{year}.nc")
        iris.save(tailwidth_cube, tailwidth_year_path)
        iris.save(p875_cube, p875_year_path)

    # Merge cubes for convenience
    try:
        merged_tailwidth = iris.cube.CubeList(tailwidth_cubes).concatenate_cube()
    except iris.exceptions.ConcatenateError:
        print("[!] Concatenate failed for tailwidth cubes. Saving as CubeList.")
        merged_tailwidth = iris.cube.CubeList(tailwidth_cubes)

    try:
        merged_p875 = iris.cube.CubeList(p875_cubes).concatenate_cube()
    except iris.exceptions.ConcatenateError:
        print("[!] Concatenate failed for p875 cubes. Saving as CubeList.")
        merged_p875 = iris.cube.CubeList(p875_cubes)

    # Save merged files
    tailwidth_path = os.path.join(member_output_dir, f"tailwidth_member{member}_merged.nc")
    p875_path = os.path.join(member_output_dir, f"p875_member{member}_merged.nc")
    iris.save(merged_tailwidth, tailwidth_path)
    iris.save(merged_p875, p875_path)

    print(f"Completed for member {member}.\nSaved in {member_output_dir}")
