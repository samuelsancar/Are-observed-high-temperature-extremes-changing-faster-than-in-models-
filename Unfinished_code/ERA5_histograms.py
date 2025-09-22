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
from Lat_Land_mask_function import lat_land_mask

# Read command-line arguments
reanalysis = 'era5'
t1 = int(sys.argv[1])
t2 = int(sys.argv[2])
#sig = sys.argv[3]
variable = 'daily_maximum_2m_temperature'

# Load trend data
datadir = '/data/users/samuel.sancar/Heat_extremes/trends/Reanalysis_trends/'
trend = xr.open_dataset(f"{datadir}{reanalysis}_{variable}_tailwidth_trends_{t1}-{t2}.nc")

# Apply latitude and land mask using function
masked_trend, _ = lat_land_mask(trend, -60)

# Extract masked_trend data
trend_var = list(masked_trend.data_vars)[0]
trend_data = masked_trend[trend_var].values.flatten()
trend_data = trend_data[~np.isnan(trend_data)]  # Remove NaNs

# Load and regrid model cubes
coupled_filenames = sorted(glob.glob('/data/users/samuel.sancar/Heat_extremes/Trends/Coupled_model_trends/*'))
sst_filenames = sorted(glob.glob('/data/users/samuel.sancar/Heat_extremes/Trends/SST_model_trends/*'))

coupled_data = [iris.load_cube(f) for f in coupled_filenames]
sst_data = [iris.load_cube(f) for f in sst_filenames]



plt.figure(figsize=(10, 6))
for f in coupled_data:
    masked_model, _ = lat_land_mask(f, -60)
    model_data = masked_model.data.flatten()
    model_data = model_data[~np.isnan(model_data)]
    plt.hist(model_data, bins=100, histtype= 'line', alpha=0.3)
    
plt.hist(trend_data_data, bins=100, histtype= 'step', color='red', alpha=0.3)
plt.xlabel('Trend magnitude (°C per decade)')
plt.ylabel('Number of grid points')
plt.title('ERA5 Tail Width Trend Histogram (1958–2022)')
plt.grid(True)
plt.tight_layout()
plt.show() 
    
    
plt.figure(figsize=(10, 6))
for f in sst_data:
    masked_model, _ = lat_land_mask(f, -60)
    model_data = masked_model.data.flatten()
    model_data = model_data[~np.isnan(model_data)]
    plt.hist(model_data, bins=100, histtype= 'line', alpha=0.3)
    
plt.hist(trend_data_data, bins=100, histtype= 'step', color='red', alpha=0.3)
plt.xlabel('Trend magnitude (°C per decade)')
plt.ylabel('Number of grid points')
plt.title('ERA5 Tail Width Trend Histogram (1958–2022)')
plt.grid(True)
plt.tight_layout()
plt.show()


lon = masked_trend['longitude'].values
lat = masked_trend['latitude'].values
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

# Extract masked_trend data
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
