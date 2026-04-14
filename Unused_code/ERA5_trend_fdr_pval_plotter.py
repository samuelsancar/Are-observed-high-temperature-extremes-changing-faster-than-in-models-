import iris
import matplotlib.pyplot as plt
import iris.plot as iplt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.cm as mpl_cm
import numpy as np
import xarray as xr
from Lat_Land_mask_function import lat_land_mask
from Save_figure import save_figure

# Parameters
reanalysis = 'era5'
t1 = 1958
t2 = 2023
variable = 'daily_maximum_2m_temperature'
datadir = '/data/users/samuel.sancar/Heat_extremes/trends/Reanalysis_trends/'

# Colour map
cmap_rb = mpl_cm.get_cmap('coolwarm', 15)

# Load data
trend = iris.load_cube(f"{datadir}{reanalysis}_{variable}_tailwidth_trends_{t1}-{t2}.nc")
template = iris.load_cube('template_cube_1deg.nc')
model_mask, _ = lat_land_mask(trend, -60, 90)

lon = template.coord('longitude').points
lat = template.coord('latitude').points
lon2d, lat2d = np.meshgrid(lon, lat)

# Load significance masks
fdr_mask = xr.DataArray(
    data=np.load(f"{datadir}{reanalysis}_{variable}_tailwidth_fdr_{t1}-{t2}.npy"),
    dims=["latitude", "longitude"],
    coords={"latitude": lat, "longitude": lon}
)
fdr_lats, fdr_lons = lat2d[fdr_mask], lon2d[fdr_mask]

pval_ds = xr.open_dataset(f"{datadir}{reanalysis}_{variable}_tailwidth_p_{t1}-{t2}.nc")
pval_mask = list(pval_ds.data_vars.values())[0].values < 0.05
pval_lats, pval_lons = lat2d[pval_mask], lon2d[pval_mask]

# Create figure with GridSpec
fig = plt.figure(figsize=(14, 12))
gs = fig.add_gridspec(2, 2, width_ratios=[20, 1], wspace=0.02, hspace=0.1)

axs = [fig.add_subplot(gs[i, 0], projection=ccrs.PlateCarree()) for i in range(2)]
cax = fig.add_subplot(gs[:, 1])  # Colourbar axis

stipple_coords = [(fdr_lons, fdr_lats), (pval_lons, pval_lats)]
labels = ['(a)', '(b)']

for ax, label, (st_lon, st_lat) in zip(axs, labels, stipple_coords):
    ax.set_extent([-180, 180, -60, 90], crs=ccrs.PlateCarree())
    mesh = iplt.pcolormesh(model_mask, cmap=cmap_rb, vmin=-0.3, vmax=0.3, axes=ax)
    ax.coastlines(resolution='50m', linewidth=0.5)
    ax.add_feature(cfeature.LAND, facecolor='lightgray')
    ax.add_feature(cfeature.OCEAN, facecolor='white')
    ax.scatter(st_lon, st_lat, color='black', s=0.5, alpha=0.5, transform=ccrs.PlateCarree())
    gl = ax.gridlines(draw_labels=True, linewidth=0.5, color='gray', alpha=0.5, linestyle='--')
    gl.top_labels = False
    gl.right_labels = False
    ax.text(-175, 85, label, fontsize=14, weight='bold', transform=ccrs.PlateCarree())
    
cbar = fig.colorbar(mesh, cax=cax, orientation='vertical', shrink = 0.3, pad=0.03)
cbar.set_label('Trend in Tail Width (°C per decade)', fontsize=12)

fig.suptitle(f"Trend in Tail Width\n{reanalysis} {variable.replace('_', ' ')} ({t1}-{t2})", fontsize=16)
save_figure(plt, 'ERA5_p&fdr_trends', 'Maps')
plt.show()


