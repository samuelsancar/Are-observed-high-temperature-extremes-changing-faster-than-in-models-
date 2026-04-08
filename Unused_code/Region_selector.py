"""
14/11/2025
----------
Region selector based on conditions 1, 2 & 6.

"""

import iris
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import iris.plot as iplt

era5 = iris.load_cube('/data/users/samuel.sancar/Heat_extremes/Trends/Reanalysis_trends/era5_daily_maximum_2m_temperature_tailwidth_trends_1958-2023.nc')
jra55 = iris.load_cube('/data/users/samuel.sancar/Heat_extremes/Trends/Reanalysis_trends/jra55__minmax_surf.015_tmax.reg_tl319._tailwidth_trends_1958-2023.nc')
era5_array = np.array(era5.data)
jra55_array = np.array(jra55.data)

era5_sig = iris.load_cube('/data/users/samuel.sancar/Heat_extremes/Trends/Reanalysis_trends/era5_daily_maximum_2m_temperature_tailwidth_p_1958-2023.nc')

era5_array = era5_array.astype(float)
jra55_array = jra55_array.astype(float)
#condition_1 = np.where(era5_sig.data < 0.05, True, False)
condition_2 = np.where((era5_array > 0) & (jra55_array > 0), True, False)

# Combine conditions
combined_conditions = condition_2

# Mask ERA5 data using combined conditions
masked_data = np.ma.masked_where(~combined_conditions, era5.data)

lats = era5.coord('latitude').points
lons = era5.coord('longitude').points

fig = plt.figure(figsize=(10, 6))
ax = plt.axes(projection=ccrs.Robinson(central_longitude=0))

ax.add_feature(cfeature.COASTLINE)
ax.add_feature(cfeature.BORDERS, linestyle=':')
ax.set_global()

mesh = ax.pcolormesh(lons, lats, masked_data, transform=ccrs.PlateCarree(),
                     cmap='Reds', shading='auto')


cb = plt.colorbar(mesh, orientation='horizontal', pad=0.05)
cb.set_label('ERA5 Tailwidth Trend (°C per decade)')

plt.title('Significant Positive Tailwidth Trends (ERA5 & JRA55)')
plt.show()





























