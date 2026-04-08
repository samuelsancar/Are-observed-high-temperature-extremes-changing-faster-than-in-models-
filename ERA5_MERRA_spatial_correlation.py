"""
Code to find spatial correlation between MERRA and ERA5.
"""

import os
import numpy as np
import iris.util
from iris.analysis import stats
from Lat_Land_mask_function import lat_land_mask

era5 = iris.load_cube('/data/users/samuel.sancar/Heat_extremes/Trends/Reanalysis_trends/era5_daily_maximum_2m_temperature_tailwidth_trends_1980-2023.nc')
merra = iris.load_cube('/data/users/samuel.sancar/Heat_extremes/Trends/Reanalysis_trends/MERRA_daily_maximum_2m_temperature_tailwidth_trends_1980-2023.nc')

era5_masked, _ = lat_land_mask(era5,-60, 90)
merra_masked, _ = lat_land_mask(merra, -60, 90)

era5_data = np.ma.masked_invalid(era5_masked.data)
jra_data = np.ma.masked_invalid(merra_masked.data)

common_mask = np.ma.getmaskarray(era5_data) | np.ma.getmaskarray(jra_data)

era5_pair = era5_masked.copy()
merra_pair = merra_masked.copy()

era5_pair.data = np.ma.array(era5_data, mask=common_mask)
merra_pair.data  = np.ma.array(jra_data,  mask=common_mask)

r = iris.analysis.stats.pearsonr(era5_pair, merra_pair, corr_coords = ('latitude', 'longitude'))

print(r.data)