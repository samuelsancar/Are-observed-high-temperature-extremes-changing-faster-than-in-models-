"""
10/11/2025
----------
Code to generate analysis plots of models.

"""

import glob
import iris
import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from Lat_Land_mask_function import lat_land_mask

# Load template and model data
template = iris.load_cube('template_cube_1deg.nc')

coupled_files = sorted(glob.glob('/data/users/samuel.sancar/Heat_extremes/Trends/Coupled_model_trends/*_tailwidth_trend.nc'))
sst_files = sorted(glob.glob('/data/users/samuel.sancar/Heat_extremes/Trends/SST_model_trends/*_tailwidth_trend.nc'))

coupled_data = [iris.load_cube(f) for f in coupled_files]
sst_data = [iris.load_cube(f) for f in sst_files]

# Combine data into array
trend_array = np.array([cube.data for cube in coupled_data + sst_data])

# Mask invalid and extreme values (both positive and negative)
trend_array = np.ma.masked_invalid(trend_array)
trend_array = np.ma.masked_outside(trend_array, -1e5, 1e5)  

# Compute mean, std, and CLT error
trend_mean = trend_array.mean(axis=0)
trend_std = trend_array.std(axis=0)
clt_error = trend_std / np.sqrt(trend_array.shape[0])

# Convert masked arrays to normal arrays for np.where
trend_mean = trend_mean.filled(np.nan)
trend_std = trend_std.filled(np.nan)
clt_error = clt_error.filled(np.nan)

print("=== DEBUG INFO AFTER MASKING ===")
print(f"trend_array shape: {trend_array.shape}")
print(f"trend_mean stats: min={np.nanmin(trend_mean):.4f}, max={np.nanmax(trend_mean):.4f}, mean={np.nanmean(trend_mean):.4f}")
print(f"trend_std stats: min={np.nanmin(trend_std):.4f}, max={np.nanmax(trend_std):.4f}")
print(f"clt_error stats: min={np.nanmin(clt_error):.4f}, max={np.nanmax(clt_error):.4f}")
print("===============================")

# Robust trend logic
robust_trends_final = np.where(
    (trend_mean > 0) & ((trend_mean - clt_error) > 0), 1,
    np.where((trend_mean < 0) & ((trend_mean + clt_error) < 0), -1, 0)
)

# Check distribution of robust trend signs
unique, counts = np.unique(robust_trends_final, return_counts=True)
print("Robust trend sign counts:", dict(zip(unique, counts)))

# Histogram of trend mean
plt.figure(figsize=(6, 4))
plt.hist(trend_mean.flatten(), bins=50, color='steelblue', edgecolor='black')
plt.title("Distribution of trend_mean (after masking)")
plt.xlabel("Trend Mean")
plt.ylabel("Frequency")
plt.show()

masked_result = lat_land_mask(robust_trends_final, -60, 90)

if isinstance(masked_result, tuple):
    data = masked_result[0]
else:
    data = masked_result

# If it's a masked array, extract data
if hasattr(data, 'data'):
    data = data.data

# Ensure numeric type
data = np.array(data, dtype=float)

print("Final data shape:", data.shape)

def plot_map(data, title, cmap='RdBu_r', vmin=None, vmax=None, label='', save_path=None):
    # Mask invalid values
    data = ma.masked_invalid(data)

    if data.ndim != 2:
        raise ValueError(f"Data must be 2D, got shape {data.shape}")
    
    lats = template.coord('latitude').points
    lons = template.coord('longitude').points
    
    # Adjust latitude range 
    if data.shape[0] != len(lats):
        lat_mask = (lats >= -60) & (lats <= 90)
        lats = lats[lat_mask]
    
    lon2d, lat2d = np.meshgrid(lons, lats)

    print(f"lon2d shape: {lon2d.shape}, lat2d shape: {lat2d.shape}, data shape: {data.shape}")
    print(f"Valid cells: {np.sum(~data.mask)}")
    if np.sum(~data.mask) == 0:
        print("Warning: All cells are masked. Map will be blank.")

    fig = plt.figure(figsize=(12, 6))
    ax = plt.axes(projection=ccrs.PlateCarree())
    c = ax.pcolormesh(lon2d, lat2d, data, cmap=cmap, vmin=vmin, vmax=vmax, shading='auto')

    ax.coastlines(resolution='50m', linewidth=1)
    ax.add_feature(cfeature.LAND, facecolor='none')
    ax.add_feature(cfeature.OCEAN, facecolor='white')
    
    cb = plt.colorbar(c, ax=ax, orientation='vertical', shrink=0.7, pad=0.05)
    cb.set_label(label)
    
    plt.title(title, fontsize=14)
    plt.tight_layout()
    plt.show()

# Plot robust trends 
plot_map(
    data,
    title='Robust Trends (Sign Only)',
    cmap='RdBu_r',
    vmin=-1,
    vmax=1,
    label='Trend Sign',
    save_path='robust_trends_map.png'
)

