
"""
17/07/25
---------
Plots regional trends in tail width calculated by ERA5_tailwidth.py for: 'Northwest Europe', 'Central China', 'Southern South America', 'Arabian Peninsula',
'Eastern Australia', 'Japan/Korea', 'Arctic Canada/Greenland', 'Northwest Canada','North Africa', 'Central Siberia', 'Eastern Siberia', 'Greenland Ice Sheet'
To run: python plot_tailwidth_trend.py [1958] [2023]
"""

import iris
import sys
import matplotlib.pyplot as plt
import iris.plot as iplt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.cm as mpl_cm
import xarray as xr
import numpy as np

reanalysis ='era5'
t1 = int(sys.argv[1])
t2 = int(sys.argv[2])
variable = 'daily_maximum_2m_temperature'

# Define color map
cmap_rb = mpl_cm.get_cmap('coolwarm', 15)

# Define regions
regions = ['europe', 'china', 'argentina', 'oman', 'australia', 'japan', 'arctic',
           'nwcanada', 'nafrica', 'central siberia', 'eastern siberia', 'gris']
regtitles = ['Northwest Europe', 'Central China', 'Southern South America', 'Arabian Peninsula',
             'Eastern Australia', 'Japan/Korea', 'Arctic Canada/Greenland', 'Northwest Canada',
             'North Africa', 'Central Siberia', 'Eastern Siberia', 'Greenland Ice Sheet']
reglats = [[45, 58], [27, 37], [-53, -35], [15, 25], [-40, -30], [30, 39],
           [78, 85], [62, 72], [20, 35], [55, 70], [65, 75], [64, 75]]
reglons = [[-5, 12], [95, 110], [-77, -66], [45, 60], [135, 150], [125, 142],
           [-95, -10], [-135, -115], [10, 35], [75, 100], [120, 140], [-50, -40]]

# Load trend data
datadir = '/data/users/samuel.sancar/Heat_extremes/trends/Reanalysis_trends/'
trend = xr.open_dataset(f"{datadir}{reanalysis}_{variable}_tailwidth_trends_{t1}-{t2}.nc")


# Download sea mask
lsm = xr.open_dataset('lsm_1279l4_0.1x0.1.grb_v4_unpack.nc')['lsm'].squeeze()
lsm_interp = lsm.interp(latitude=trend.latitude, longitude=trend.longitude, method="nearest")

# Download Significance stippling
lon = trend['longitude'].values
lat = trend['latitude'].values
lon2d, lat2d = np.meshgrid(lon, lat)

statsig = xr.DataArray(data=np.load(f"{datadir}{reanalysis}_{variable}_tailwidth_fdr_{t1}-{t2}.npy"),dims=["latitude", "longitude"],coords={"latitude": lat, "longitude": lon},name="significance_mask")
stipple_lats, stipple_lons = lat2d[statsig], lon2d[statsig]

# Apply mask: keep land (1), mask sea (0)
masked_trends = trend.where(lsm_interp > 0.5)

data_var = list(masked_trends.data_vars)[0]

# Plot each region
for i, region in enumerate(regions):
    fig, ax = plt.subplots(figsize=(8, 6), subplot_kw={'projection': ccrs.PlateCarree()})
    ax.set_extent([reglons[i][0], reglons[i][1], reglats[i][0], reglats[i][1]], crs=ccrs.PlateCarree())
    mesh = masked_trends[data_var].plot.pcolormesh(ax=ax,cmap=cmap_rb,vmin=-0.4,vmax=0.4,add_colorbar=False,transform=ccrs.PlateCarree())
    cbar = plt.colorbar(mesh, orientation='horizontal', pad=0.05, aspect=50)
    cbar.set_label('Trend in Tail Width (°C per decade)', fontsize=12)
    ax.scatter(stipple_lons, stipple_lats, s=1, color='black', transform=ccrs.PlateCarree())
    ax.coastlines(resolution='50m', linewidth=0.5); ax.add_feature(cfeature.NaturalEarthFeature(category='cultural',name='admin_0_boundary_lines_land',facecolor='none',scale='50m',edgecolor='0',lw=.05)); ax.add_feature(cfeature.LAKES.with_scale('110m'),edgecolor='0',facecolor='none',lw=.05); ax.add_feature(cfeature.STATES.with_scale('110m'),edgecolor='0',facecolor='none',lw=.2); ax.add_feature(cfeature.RIVERS.with_scale('110m'),edgecolor='tab:blue',facecolor='none',lw=.4)
    ax.add_feature(cfeature.BORDERS, linewidth=0.5)
    ax.add_feature(cfeature.LAND, facecolor='lightgray')
    ax.add_feature(cfeature.OCEAN, facecolor='white')

    plt.title(f"{regtitles[i]} ({region})")
    plt.tight_layout()
    plt.show()

