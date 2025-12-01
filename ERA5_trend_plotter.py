"""
17/07/25
---------
Plot trend in tail width calculated by ERA5_tailwidth.py
To run: python plot_tailwidth_trend.py [1958] [2023] [fdr|p]
"""

"""
17/07/25
---------
Plot trend in tail width calculated by ERA5_tailwidth.py
To run: python plot_tailwidth_trend.py [1958] [2023] [fdr|p]
"""

import iris
import sys
import os
import matplotlib.pyplot as plt
import iris.plot as iplt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.cm as mpl_cm
import numpy as np
import xarray as xr
from Lat_Land_mask_function import lat_land_mask
from ERA5_MKtest_Bootstrap import significance_plot

# Define colour map
cmap_rb = mpl_cm.get_cmap('coolwarm', 15)

reanalysis = 'era5'
variable = 'daily_maximum_2m_temperature'
t1 = 1958
t2 = 2023
sig = sys.argv[1] # either 'fdr' or 'pvals'
alpha = float(sys.argv[2])
lat_minimum = -60
lat_maximum = int(sys.argv[3])


# Load trend data
datadir = '/data/users/samuel.sancar/Heat_extremes/Trends/Reanalysis_trends/'
trend = iris.load_cube(f"{datadir}{reanalysis}_{variable}_tailwidth_trends_{t1}-{t2}.nc")
p87p5_trend = iris.load_cube(f"{datadir}{reanalysis}_{variable}_p87p5_trends_{t1}-{t2}.nc")
template = iris.load_cube('template_cube_1deg.nc')

# Apply latitude and land mask using function
model_mask, _ = lat_land_mask(trend, lat_minimum, lat_maximum)
p87p5_masked, _ = lat_land_mask(p87p5_trend, lat_minimum, lat_maximum)

# Extract grid point locations for stippling
lon = template.coord('longitude').points
lat = template.coord('latitude').points
lon2d, lat2d = np.meshgrid(lon, lat)

# Load significance mask (already masked for land and latitude)
if sig == 'fdr':
    statsig = xr.DataArray(
        data=np.load(f"{datadir}{reanalysis}_{variable}_tailwidth_fdr_{t1}-{t2}.npy"),
        dims=["latitude", "longitude"],
        coords={"latitude": lat, "longitude": lon},
        name="significance_mask"
    )
    stipple_lats, stipple_lons = lat2d[statsig], lon2d[statsig]
else:
    statsig_ds, _, _ = significance_plot(lat_minimum, lat_maximum)
    statsig = statsig_ds[0]
    sig_mask_np = statsig_ds < alpha
    stipple_lats, stipple_lons = lat2d[sig_mask_np], lon2d[sig_mask_np]

neg_indices = np.where(p87p5_masked.data < 0)
neg_lats = lat2d[neg_indices]
neg_lons = lon2d[neg_indices]

fig = plt.figure(figsize=(12, 6))
ax = plt.axes(projection=ccrs.Robinson(central_longitude=0))
ax.set_extent([-180, 180, lat_minimum, lat_maximum], crs=ccrs.PlateCarree(central_longitude=0))

iplt.pcolormesh(model_mask, cmap=cmap_rb, vmin=-0.3, vmax=0.3)

ax.coastlines(resolution='50m', linewidth=0.5)
ax.add_feature(cfeature.NaturalEarthFeature(category='cultural', name='admin_0_boundary_lines_land', facecolor='none', scale='50m', edgecolor='0', lw=.05))
ax.add_feature(cfeature.LAKES.with_scale('110m'), edgecolor='0', facecolor='none', lw=.05)
ax.add_feature(cfeature.STATES.with_scale('110m'), edgecolor='0', facecolor='none', lw=.2)
ax.add_feature(cfeature.RIVERS.with_scale('110m'), edgecolor='tab:blue', facecolor='none', lw=.4)
ax.add_feature(cfeature.LAND, facecolor='lightgray')
ax.add_feature(cfeature.OCEAN, facecolor='white')

# Plot stippling
ax.scatter(stipple_lons, stipple_lats, color='black', s=0.4, alpha=0.3, transform=ccrs.PlateCarree(central_longitude=0), label='Stippling')
# Overlay white dots where p87p5_trend is negative
ax.scatter(neg_lons, neg_lats, color='white', s=0.9, alpha=0.6, transform=ccrs.PlateCarree(), label='p87p5 trend < 0')

gl = ax.gridlines(draw_labels=True, linewidth=0.5, color='gray', alpha=0.5, linestyle='--')
gl.top_labels = False
gl.right_labels = False

cbar = plt.colorbar(ax.collections[0], orientation='horizontal', pad=0.05, aspect=50)
cbar.set_label('Trend in Tailwidth (°C per decade)', fontsize=12)

plt.title(f"Trend in Tailwidth for ERA5 reanalysis (1958-2023)", fontsize=14, weight='bold')

plt.tight_layout()
plt.show()


