"""
24/07/25
---------

Code that applies MannKendall test on yearly ERA5 tailwidths.
Output is KendallTau and p-values (for all trends and positive only trends).

To run: python ERA5_MKtest.py 
"""

import numpy as np
import iris
from MK_functions import load_yearly_data, create_time_coord,stack_cubes,compute_kendall_tau, save_results, save_positive_p_values

# Beginning & end year+1
t1 = 1958
t2 = 2023

lat_min = -60
lat_max = 90

years = np.arange(t1, t2)    
    
reanalysis = 'era5'
variable = 'daily_maximum_2m_temperature'
datadir = '/data/users/samuel.sancar/Heat_extremes/ERA5_tailwidths/'
outdir = '/data/users/samuel.sancar/Heat_extremes/Trends/Reanalysis_trends/'
p87p5 = '/data/users/samuel.sancar/Heat_extremes/Trends/Reanalysis_trends/'

template = iris.load_cube('template_cube_1deg.nc')
p87p5_cube = iris.load_cube(f"{p87p5}era5_daily_maximum_2m_temperature_p87p5_trends_{t1}-{t2-1}.nc")

# Prepare cubes
time_coord = create_time_coord(years)
cubes = load_yearly_data(datadir, reanalysis, variable, years)
data_cube = stack_cubes(cubes, time_coord)
    
# Compute p-values
tau_matrix, p_matrix, lat_points, lon_points = compute_kendall_tau(years, data_cube)

# Mask regions with negative trends in p87p5 percentile.
p_matrix_masked = np.ma.masked_where(p87p5_cube.data < 0, p_matrix)
tau_matrix_masked = np.ma.masked_where(p87p5_cube.data < 0, tau_matrix)

p_out  = np.ma.filled(p_matrix_masked, np.nan)
tau_out = np.ma.filled(tau_matrix_masked, np.nan)

#Save p-values
save_results(tau_out, p_out, lat_points, lon_points, outdir, reanalysis, variable, t1, t2)
save_positive_p_values(tau_out, p_out, lat_points, lon_points, outdir, reanalysis, variable, t1, t2)



