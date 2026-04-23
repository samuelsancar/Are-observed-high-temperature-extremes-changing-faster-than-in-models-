""" 
Code to plot contour plots of obs and model data.
"""

import iris
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import numpy as np
import os
import iris.plot as iplt
from Lat_Land_mask_function import lat_land_mask
import glob

era5 = iris.load_cube('/data/users/samuel.sancar/Heat_extremes/Trends/Reanalysis_trends/era5_daily_maximum_2m_temperature_tailwidth_trends_1958-2023.nc')
jra = iris.load_cube('/data/users/samuel.sancar/Heat_extremes/Trends/Reanalysis_trends/jra55__minmax_surf.015_tmax.reg_tl319._tailwidth_trends_1958-2023.nc')

era5_masked, _ = lat_land_mask(era5, -60, 90)
jra_masked, _ = lat_land_mask(jra, -60, 90)

member = np.arange(1,41)
canari_datadir = '/data/users/samuel.sancar/Heat_extremes/Trends/CANARI_trends/'
coupled_datadir = sorted(glob.glob('/data/users/samuel.sancar/Heat_extremes/Trends/Coupled_model_trends/*_tailwidth_trend.nc'))
sst_datadir = sorted(glob.glob('/data/users/samuel.sancar/Heat_extremes/Trends/SST_model_trends/*_tailwidth_trend.nc'))

fig = plt.figure(figsize=(10,8))
ax = plt.axes(projection=ccrs.Robinson(central_longitude=0))
iplt.contour(era5_masked,levels=[-0.2,0.2],colors='black',linewidths = [3])
iplt.contour(jra_masked,levels=[-0.2,0.2],colors='black',linewidths=[3])

for i in member:
    print(i)
    canari_cube = iris.load_cube(os.path.join(canari_datadir, f'{i}_tailwidth_trend.nc'))
    cube_masked, _ = lat_land_mask(canari_cube, -60, 90)
    iplt.contour(cube_masked,levels=[-0.2,0.2],colors=['blue','red'])

ax = plt.gca()
ax.coastlines()
plt.title('Tailwidth trends (1958-2022)')
plt.tight_layout()
plt.show()
plt.close()


fig = plt.figure(figsize=(10,8))
ax = plt.axes(projection=ccrs.Robinson(central_longitude=0))
iplt.contour(era5_masked,levels=[-0.2,0.2],colors='black',linewidths = [3])
iplt.contour(jra_masked,levels=[-0.2,0.2],colors='black',linewidths=[3])

for i in coupled_datadir:
    print(i)
    coupled_cube = iris.load_cube(i)
    cube_masked, _ = lat_land_mask(coupled_cube, -60,90)
    iplt.contour(cube_masked,levels=[-0.2,0.2],colors=['blue','red'])
    
for i in sst_datadir:
    print(i)
    sst_cube = iris.load_cube(i)
    cube_masked, _ = lat_land_mask(sst_cube, -60,90)
    iplt.contour(cube_masked,levels=[-0.2,0.2],colors=['blue','red'])

ax = plt.gca()
ax.coastlines()
plt.title('Tailwidth trends (1958-2022)')
plt.tight_layout()
plt.show()