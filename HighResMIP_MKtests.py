"""
10/08/2025
----------

Code that applies MannKendall test on HighResMIP tail widths.
Output is KendallTau and p-values (for all trends and positive only trends).

To run: python HighResMIP_MKtests.py
"""

import iris
import numpy as np
import pandas as pd
import os
import xarray as xr
from iris.coords import AuxCoord
from MK_functions import compute_kendall_tau

# Beginning & end year+1
t1 = 1958
t2 = 2023

years = np.arange(t1, t2)

lat_minimum = -60
lat_maximum = 90

Model_dir = '/data/users/samuel.sancar/Heat_extremes/Tail_widths/'
outdir = '/data/users/samuel.sancar/Heat_extremes/Trends/Trend_significance/P-vals/'

template = iris.load_cube('template_cube_1deg.nc')

# Read in excel sheet containing model names
df = pd.read_excel("/data/users/samuel.sancar/Model_table.xlsx", sheet_name='Model Stack', engine="openpyxl")

# ensure consistent column names (strip spaces)
df.columns = df.columns.str.strip()

processed = 0

for row in np.arange(0,63):
    
    row_idx = row  
    row = df.iloc[row_idx]

    # Map columns to variables 
    Mode       = str(row["Coupled/SST"]).strip()
    model      = str(row["Model"]).strip()
    experiment = str(row["Base experiment"]).strip()
    member     = str(row["Member"]).strip()

    print(model, experiment, member)
    
    cube = iris.load_cube(f"{Model_dir}{model}_{experiment}_{member}_tailwidths_1958-{t2-1}.nc")

    p87p5_dir = f"/data/users/samuel.sancar/Heat_extremes/Trends/{Mode}_model_trends/"
    p87p5_cube = iris.load_cube(f"{p87p5_dir}{model}_{experiment}_{member}_p87p5_trend_{t1}-{t2-1}.nc")

    cube_start_year = 1958 # Define years that tail width file covers (Keep the same even if wanting p-vals over years 1980-2022)
    cube_end_year   = 2022

    full_years = np.arange(cube_start_year, cube_end_year + 1)

    mask = (full_years >= t1) & (full_years <= (t2 - 1)) # Define mask for years to calc trend over (ie, if t1 = 1980)

    cube = cube[mask, :, :] # Mask years for trend calc.
    years = full_years[mask]
        
    # Compute p-values
    tau_matrix, p_matrix, lat_points, lon_points = compute_kendall_tau(years, cube)
        
    # Mask regions with negative trends in p87p5 percentile.
    p_matrix_masked = np.ma.masked_where(p87p5_cube.data < 0, p_matrix)
    tau_matrix_masked = np.ma.masked_where(p87p5_cube.data < 0, tau_matrix)
    
    p_out  = np.ma.filled(p_matrix_masked, np.nan)
    tau_out = np.ma.filled(tau_matrix_masked, np.nan)
           
    tau_save = xr.DataArray(tau_out, coords={"latitude": lat_points, "longitude": lon_points}, dims=["latitude", "longitude"], name = "Tailwidth_mktau")
    p_save = xr.DataArray(p_out, coords={"latitude": lat_points, "longitude": lon_points}, dims=["latitude", "longitude"], name = "Tailwidth_trend_p-values")
    tau_save.to_netcdf(os.path.join(outdir, f"{model}_{experiment}_{member}_mktau_{t1}-{t2-1}.nc"))
    p_save.to_netcdf(os.path.join(outdir, f"{model}_{experiment}_{member}_p_{t1}-{t2-1}.nc"))
    
    pos_p_matrix = np.where(tau_out > 0, p_out, np.nan)
    pos_p_save = xr.DataArray(pos_p_matrix, coords={"latitude": lat_points, "longitude": lon_points}, dims=["latitude", "longitude"], name="p_value_given_positive_tau")
    pos_p_save.to_netcdf(os.path.join(outdir, f"{model}_{experiment}_{member}_pos_p_{t1}-{t2-1}.nc"))
            
    processed += 1

    print(f"Finished processing. Total files successfully processed: {processed}")









