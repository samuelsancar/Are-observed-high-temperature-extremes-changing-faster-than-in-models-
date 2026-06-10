"""
1/11/25
---------

Code that applies MannKendall test on CANARI tail widths.
Output is KendallTau and p-values (for all trends and positive only trends).

To run: python CANARI_MKtest.py
"""

import iris
import numpy as np
import os
import xarray as xr
from MK_functions import compute_kendall_tau

# Beginning & end year+1
t1 = 1958
t2 = 2023

lat_minimum = -60
lat_maximum = 90

years = np.arange(t1, t2)

Model_dir = '/data/users/samuel.sancar/Heat_extremes/Tail_widths/canari/'
outdir = '/data/users/samuel.sancar/Heat_extremes/Trends/Trend_significance/CANARI/P-vals/'

template = iris.load_cube('template_cube_1deg.nc')

processed = 1

for i in np.arange(1,41): # 40 members
    
    print(f"Processing member: {i}")

    if t1 == 1958:
        cubes = iris.load_cube(f'{Model_dir}{i}_tailwidths_{t1}-{t2-1}.nc')
        
    else:
        cubes = iris.load_cube(f'{Model_dir}{i}_tailwidths_{t1}-{t2-1}.nc')
        
    p87p5_dir = f'/data/users/samuel.sancar/Heat_extremes/Trends/CANARI_trends/'
    p87p5_cube = iris.load_cube(p87p5_dir + f'{i}_p87p5_trend_{t1}-{t2-1}.nc')
            
    tau_matrix, p_matrix, lat_points, lon_points = compute_kendall_tau(years, cubes)   
    
    # Mask regions with negative trends in p87p5 percentile.
    p_matrix_masked = np.ma.masked_where(p87p5_cube.data < 0, p_matrix)
    tau_matrix_masked = np.ma.masked_where(p87p5_cube.data < 0, tau_matrix)
    
    p_out  = np.ma.filled(p_matrix_masked, np.nan)
    tau_out = np.ma.filled(tau_matrix_masked, np.nan)    
      
    tau_save = xr.DataArray(tau_out, coords={"latitude": lat_points, "longitude": lon_points}, dims=["latitude", "longitude"], name = "Tailwidth_mktau")
    p_save = xr.DataArray(p_out, coords={"latitude": lat_points, "longitude": lon_points}, dims=["latitude", "longitude"], name = "Tailwidth_trend_p-values")
    tau_save.to_netcdf(os.path.join(outdir, f"{i}_mktau_{t1}-{t2-1}.nc"))
    p_save.to_netcdf(os.path.join(outdir, f"{i}_p_{t1}-{t2-1}.nc"))
    
    pos_p_matrix = np.where(tau_out > 0, p_out, np.nan)
    pos_p_save = xr.DataArray(pos_p_matrix, coords={"latitude": lat_points, "longitude": lon_points}, dims=["latitude", "longitude"], name="p_value_given_positive_tau")
    pos_p_save.to_netcdf(os.path.join(outdir, f"{i}_pos_p_{t1}-{t2-1}.nc"))

    print(f"Finished processing. Total files successfully processed: {processed}/40")
    
    processed += 1


