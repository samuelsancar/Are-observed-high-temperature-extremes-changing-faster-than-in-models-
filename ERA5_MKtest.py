"""

24/07/25
---------
Code that applies MannKendall test on ERA5 tailwidth trends. Output is KendallTau and p-values.

To run: python ERA5_MK_tests.py [1958] [2023]

"""

from MK_functions import load_and_mask_yearly_data, create_time_coord, stack_cubes, create_year_cube, compute_kendall_tau
import sys
import numpy as np
import xarray as xr
import os

def save_results(tau_matrix, p_matrix, lat_points, lon_points, outdir, reanalysis, variable, t1, t2):
    tau_da = xr.DataArray(tau_matrix, coords={"latitude": lat_points, "longitude": lon_points}, dims=["latitude", "longitude"])
    p_da = xr.DataArray(p_matrix, coords={"latitude": lat_points, "longitude": lon_points}, dims=["latitude", "longitude"])
    os.makedirs(outdir, exist_ok=True)
    tau_da.to_netcdf(os.path.join(outdir, f"{reanalysis}_{variable}_tailwidth_mktau_{t1}-{t2}.nc"))
    p_da.to_netcdf(os.path.join(outdir, f"{reanalysis}_{variable}_tailwidth_p_{t1}-{t2}.nc"))

def main():
    datadir = '/data/users/samuel.sancar/Heat_extremes/ERA5_tailwidths/'
    outdir = '/data/users/samuel.sancar/Heat_extremes/Trends/Reanalysis_trends'
    reanalysis = 'era5'
    variable = 'daily_maximum_2m_temperature'
    t1 = int(sys.argv[1])
    t2 = int(sys.argv[2])
    years = np.arange(t1, t2)

    time_coord = create_time_coord(years)
    cubes = load_and_mask_yearly_data(datadir, reanalysis, variable, years)
    data_cube = stack_cubes(cubes, time_coord)
    year_cube = create_year_cube(data_cube, years)
    tau_matrix, p_matrix, lat_points, lon_points = compute_kendall_tau(year_cube, data_cube)
    save_results(tau_matrix, p_matrix, lat_points, lon_points, outdir, reanalysis, variable, t1, t2)

if __name__ == "__main__":
    main()