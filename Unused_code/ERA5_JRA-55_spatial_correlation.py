"""
Code to find spatial correlation between JRA-55 and ERA5.
"""

import os
import numpy as np
import iris.util
from iris.analysis import stats
from Lat_Land_mask_function import lat_land_mask

era5 = iris.load_cube('/data/users/samuel.sancar/Heat_extremes/Trends/Reanalysis_trends/era5_daily_maximum_2m_temperature_tailwidth_trends_1958-2023.nc')
jra = iris.load_cube('/data/users/samuel.sancar/Heat_extremes/Trends/Reanalysis_trends/jra55__minmax_surf.015_tmax.reg_tl319._tailwidth_trends_1958-2023.nc')

era5_masked, _ = lat_land_mask(era5,-60, 90)
jra_masked, _ = lat_land_mask(jra, -60, 90)

era5_data = np.ma.masked_invalid(era5_masked.data)
jra_data = np.ma.masked_invalid(jra_masked.data)

common_mask = np.ma.getmaskarray(era5_data) | np.ma.getmaskarray(jra_data)

era5_pair = era5_masked.copy()
jra_pair = jra_masked.copy()

era5_pair.data = np.ma.array(era5_data, mask=common_mask)
jra_pair.data  = np.ma.array(jra_data,  mask=common_mask)

r = iris.analysis.stats.pearsonr(era5_pair, jra_pair, corr_coords = ('latitude', 'longitude'))

print(r.data)