"""
In progress

"""

import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
import glob
import iris
import iris.analysis as ia
import iris.analysis.cartography as iac  # Correct import
from iris.coord_categorisation import add_year
import warnings

def get_area_weights(cube):
    for coord in ['latitude','longitude']:
        if cube.coord(coord).has_bounds() is False:
            cube.coord(coord).guess_bounds()
    weights = iac.area_weights(cube)
    return weights

def areaMean(cube):
    # Only collapse over latitude and longitude if they exist
    collapse_dims = [dim for dim in ['latitude', 'longitude'] if dim in [coord.name() for coord in cube.coords()]]
    if collapse_dims:
        weights = get_area_weights(cube)
        cube = cube.collapsed(collapse_dims, ia.MEAN, weights=weights)
    return cube

# Configuration
variable = 'daily_maximum_2m_temperature'
climate_models = '/data/users/samuel.sancar/Heat_extremes/trends/Coupled_model_trends/*'
era5_file = '/data/users/samuel.sancar/Heat_extremes/trends/Reanalysis_trends/era5_daily_maximum_2m_temperature_tailwidth_trends_1958-2023.nc'

models = glob.glob(climate_models)
era5 = iris.load_cube(era5_file)

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

# Prepare subplots
fig, axes = plt.subplots(3, 4, figsize=(20, 12), sharex=True, sharey=True)
axes = axes.flatten()

for i, region_name in enumerate(regions):
    ax = axes[i]
    for model in models:
        trend = iris.load_cube(model)
        lsm = xr.open_dataset('lsm_1279l4_0.1x0.1.grb_v4_unpack.nc')['lsm'].squeeze()
        lat_points = trend.coord('latitude').points
        lon_points = trend.coord('longitude').points
        lsm_interp = lsm.interp(latitude=lat_points, longitude=lon_points, method="nearest")
        trend.data = np.ma.masked_where(lsm_interp.values <= 0.5, trend.data)
        masked_trend = trend
        region = masked_trend.intersection(latitude=(reglats[i][0], reglats[i][1]),
                                           longitude=(reglons[i][0], reglons[i][1]),
                                           ignore_bounds=True)
        try:
            if region is None or region.data.size == 0:
                print(f"Empty region for model {model} in region {region_name}")
                continue
            avg_trend = areaMean(region)
            if avg_trend is None or avg_trend.data.size == 0:
                print(f"Empty avg_trend for model {model} in region {region_name}")
                continue
            # Check if avg_trend has a time dimension
            if not avg_trend.coords('time') or avg_trend.coord('time').ndim == 0 or avg_trend.data.ndim == 0:
                print(f"avg_trend is scalar or missing time for model {model} in region {region_name}")
                continue
            time_coord = avg_trend.coord('time')
            years = np.array([dt.year for dt in time_coord.units.num2date(time_coord.points)])
            data = avg_trend.data
            mask = (years >= 1958) & (years <= 2022)
            years = years[mask]
            data = data[mask]
            valid = (~np.ma.getmaskarray(data)) & (~np.isnan(data))
            if np.sum(valid) > 1:
                fit_coeffs = np.polyfit(years[valid], data[valid], 1)
                fit_line = np.polyval(fit_coeffs, years[valid])
                ax.plot(years[valid], fit_line, label=model.split('/')[-1], linewidth=1.0)
            else:
                print(f"Not enough valid data for model {model} in region {region_name}")
        except Exception as e:
            print(f"Error for model {model} in region {region_name}: {e}")
    ax.set_title(regtitles[i])
    ax.grid(True)

# Set common labels and legend
fig.text(0.5, 0.04, 'Year', ha='center', fontsize=14)
fig.text(0.04, 0.5, 'Trend', va='center', rotation='vertical', fontsize=14)
handles, labels = axes[0].get_legend_handles_labels()
fig.legend(handles, labels, loc='upper center', ncol=4, fontsize=10)
plt.tight_layout(rect=[0, 0.03, 1, 0.95])
plt.show()








