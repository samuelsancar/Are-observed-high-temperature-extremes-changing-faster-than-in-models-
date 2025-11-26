import iris
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import iris.plot as iplt
import sys

variable = sys.argv[1]

if variable == 'random':
    members = np.random.choice(np.arange(1, 41, 1), size=2, replace=False)
    member1, _ = members[0], members[1]
else:
    member1 = variable

# Import Reanalysis trends
member1_trend = iris.load_cube(f'/data/users/samuel.sancar/Heat_extremes/Trends/CANARI_trends/{member1}_tailwidth_trend.nc')
jra55_trend = iris.load_cube(f'/data/users/samuel.sancar/Heat_extremes/Trends/Reanalysis_trends/jra55__minmax_surf.015_tmax.reg_tl319._tailwidth_trends_1958-2023.nc')
member1_agreement = np.load(f'/data/users/samuel.sancar/Heat_extremes/Trends/CANARI_trends/{member1}_trend_agreement.npy')
member1_array = np.array(member1_trend.data)
member2_array = np.array(jra55_trend.data)

# Import ERA5 significant trends
member2_sig = iris.load_cube(f'/data/users/samuel.sancar/Heat_extremes/Trends/Trend_significance/CANARI/P-vals/{member1}_tailwidths_p_1958.nc')

member1_array = member1_array.astype(float)
member2_array = member2_array.astype(float)
condition_1 = np.where(member2_sig.data < 0.05, True, False)
condition_2 = np.where((member1_array > 0) & (member2_array > 0), True, False)
condition_3 = np.where(member1_agreement>=0.9, True, False)

# Combine conditions
combined_conditions = condition_1 & condition_2 & condition_3

# Mask ERA5 data using combined conditions
masked_data = np.ma.masked_where(~combined_conditions, member1_trend.data)

lats = member1_trend.coord('latitude').points
lons = member1_trend.coord('longitude').points

fig = plt.figure(figsize=(10, 6))
ax = plt.axes(projection=ccrs.Robinson(central_longitude=0))

# Add map features
ax.add_feature(cfeature.COASTLINE)
ax.add_feature(cfeature.BORDERS, linestyle=':')
ax.set_global()

# Plot masked data
mesh = ax.pcolormesh(lons, lats, masked_data, transform=ccrs.PlateCarree(),
                     cmap='Reds', shading='auto')

# Add colourbar
cb = plt.colorbar(mesh, orientation='horizontal', pad=0.05)
cb.set_label('ERA5 Tailwidth Trend (°C per decade)')

plt.title(f'Regions passing conditions for CANARI member {member1}, and positive trends for member JRA-55')
plt.show()