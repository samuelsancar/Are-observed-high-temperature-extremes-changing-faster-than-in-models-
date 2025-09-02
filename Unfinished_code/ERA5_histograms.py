"""
17/07/25
---------
Plot histogram in tail widths trends of ERA5 data.
To run: python plot_tailwidth_trend.py [1958] [2023]
"""

import iris
import sys
import matplotlib.pyplot as plt
import matplotlib.cm as mpl_cm
import numpy as np
import xarray as xr

# Read command-line arguments
reanalysis = 'era5'
t1 = int(sys.argv[1])
t2 = int(sys.argv[2])
sig = sys.argv[3]
variable = 'daily_maximum_2m_temperature'

# Load trend data
datadir = '/data/users/samuel.sancar/Heat_extremes/trends/Reanalysis_trends/'
trend = xr.open_dataset(f"{datadir}{reanalysis}_{variable}_tailwidth_trends_{t1}-{t2}.nc")
template = xr.open_dataset('template_cube_1deg.nc')

# Regrid ERA5 trend
trend_interp = trend.interp(latitude=template.latitude, longitude=template.longitude, method="nearest")

# Download sea mask
lsm = xr.open_dataset('lsm_1279l4_0.1x0.1.grb_v4_unpack.nc')['lsm'].squeeze()
lsm_interp = lsm.interp(latitude=template.latitude, longitude=template.longitude, method="nearest")

# Create a 2D mask for latitudes >= -60 (exclude Antarctica)
latitudes = template.latitude.values
lat_mask = latitudes >= -60
lat_mask_2d = np.repeat(lat_mask[:, np.newaxis], lsm_interp.shape[1], axis=1)

# Apply land and latitude mask
combined_mask = (lsm_interp > 0.5) & lat_mask_2d
masked_trend = trend_interp.where(combined_mask)

# Extract trend data
trend_var = list(masked_trend.data_vars)[0]
trend_data = masked_trend[trend_var].values.flatten()
trend_data = trend_data[~np.isnan(trend_data)]  # Remove NaNs

# Plot histogram
plt.figure(figsize=(10, 6))
plt.hist(trend_data, bins=100, range=(trend_data.min(), trend_data.max()))
#plt.yscale('log')
plt.xlabel('Trend magnitude (°C per decade)')
plt.ylabel('Number of grid points')
plt.title('ERA5 Tail Width Trend Histogram (1958–2023)')
plt.grid(True)
plt.tight_layout()
plt.show()

lon = trend['longitude'].values
lat = trend['latitude'].values
lon2d, lat2d = np.meshgrid(lon, lat)

# Load significance mask based on user input
if sig == 'fdr':
    statsig = xr.DataArray(data=np.load(f"{datadir}{reanalysis}_{variable}_tailwidth_fdr_{t1}-{t2}.npy"),dims=["latitude", "longitude"],coords={"latitude": lat, "longitude": lon},name="significance_mask")
    stipple_lats, stipple_lons = lat2d[statsig], lon2d[statsig]
else:
    statsig_ds = xr.open_dataset(f"{datadir}{reanalysis}_{variable}_tailwidth_p_{t1}-{t2}.nc")
    statsig_list = list(statsig_ds.data_vars.values())[0]
    sig_mask_np = statsig_list.values < 0.05 
    statsig = (lsm_interp > 0.5) & sig_mask_np
    stipple_lats, stipple_lons = lat2d[statsig], lon2d[statsig]

stat_sig_points = masked_trend.where(statsig)

# Extract trend data
sig_trend_var = list(stat_sig_points.data_vars)[0]
sig_trend_data = stat_sig_points[sig_trend_var].values.flatten()
sig_trend_data = sig_trend_data[~np.isnan(sig_trend_data)]  # Remove NaNs

# Plot histogram
plt.figure(figsize=(10, 6))
plt.hist(sig_trend_data, bins=50, range=(trend_data.min(), trend_data.max()))
#plt.yscale('log')
plt.xlabel('Trend magnitude (°C per decade)')
plt.ylabel('Number of grid points')
plt.title('ERA5 statistically significant tailwidth trends histogram (1958–2023)')
plt.grid(True)
plt.tight_layout()
plt.show()
