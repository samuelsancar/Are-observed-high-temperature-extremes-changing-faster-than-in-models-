"""
24/07/25
---------
Plot and save Stat_sig (p-values) in MannKendall results for ERA5 tail width trends

To run: python ERA5_MKstatsig_plot.py [1958] [2023]

"""

import os
import sys
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib.cm as mpl_cm
import cartopy.crs as ccrs
import cartopy.feature as cfeature

def load_p_values(datadir, reanalysis, variable, t1, t2):
    filename = f"{reanalysis}_{variable}_tailwidth_p_{t1}-{t2}.nc"
    ds = xr.open_dataset(os.path.join(datadir, filename))
    # Extract the first variable in the dataset
    data_var = list(ds.data_vars.values())[0]
    return data_var

def create_stippling_mask(p_values, threshold=0.05):
    return p_values.values < threshold

def plot_stippling(lat, lon, mask, reanalysis, variable, t1, t2):
    lon2d, lat2d = np.meshgrid(lon, lat)
    stipple_lats = lat2d[mask]
    stipple_lons = lon2d[mask]

    fig = plt.figure(figsize=(12, 6))
    ax = plt.axes(projection=ccrs.PlateCarree())

    ax.scatter(stipple_lons, stipple_lats, color='black', s=0.5, alpha=0.5,
               transform=ccrs.PlateCarree(), label='Stippling')

    ax.coastlines(resolution='50m', linewidth=1)
    ax.add_feature(cfeature.LAND, facecolor='white')
    ax.add_feature(cfeature.OCEAN, facecolor='white')

    gl = ax.gridlines(draw_labels=True, linewidth=0.5, color='gray',
                      alpha=0.5, linestyle='--')
    gl.top_labels = False
    gl.right_labels = False

    title = f"p < 0.05 in tailwidth trends for\n{reanalysis} {variable.replace('_', ' ')} ({t1}-{t2})"
    plt.title(title, fontsize=14, weight='bold')
    plt.tight_layout()
    plt.show()

def main():
    reanalysis = 'era5'
    t1 = int(sys.argv[1])
    t2 = int(sys.argv[2])
    variable = 'daily_maximum_2m_temperature'

    datadir = '/data/users/samuel.sancar/Heat_extremes/trends/Reanalysis_trends/'
    out_path = datadir

    p_values = load_p_values(datadir, reanalysis, variable, t1, t2)
    sig_mask_np = create_stippling_mask(p_values)

    np.save(os.path.join(out_path, f"{reanalysis}_{variable}_tailwidth_p_mask_{t1}-{t2}.npy"), sig_mask_np)

    lat = p_values['latitude'].values
    lon = p_values['longitude'].values
    plot_stippling(lat, lon, sig_mask_np, reanalysis, variable, t1, t2)

if __name__ == "__main__":
    main()
