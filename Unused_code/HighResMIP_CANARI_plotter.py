
import sys
import numpy as np
import iris
import iris.plot as iplt
import matplotlib.pyplot as plt
import matplotlib.cm as mpl_cm
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import xarray as xr

sys.path.append('/home/users/samuel.sancar/Heat_Extremes/Code/Functions/')
from Lat_Land_mask_function import lat_land_mask
from Save_figure_function import save_figure

LAT_MIN, LAT_MAX = -60, 90
CMAP = mpl_cm.get_cmap('coolwarm', 15)
VMIN, VMAX = -0.3, 0.3
CENTRAL_LON = 0

def _add_common_features(ax):
    ax.coastlines(resolution='50m', linewidth=0.5)
    ax.add_feature(cfeature.NaturalEarthFeature(category='cultural',
                                                name='admin_0_boundary_lines_land',
                                                facecolor='none', scale='50m',
                                                edgecolor='0', lw=.05))
    ax.add_feature(cfeature.LAKES.with_scale('110m'), edgecolor='0', facecolor='none', lw=.05)
    ax.add_feature(cfeature.STATES.with_scale('110m'), edgecolor='0', facecolor='none', lw=.2)
    ax.add_feature(cfeature.RIVERS.with_scale('110m'), edgecolor='tab:blue', facecolor='none', lw=.4)
    ax.add_feature(cfeature.LAND, facecolor='lightgray')
    ax.add_feature(cfeature.OCEAN, facecolor='white')
    gl = ax.gridlines(draw_labels=True, linewidth=0.5, color='gray', alpha=0.5, linestyle='--')
    gl.xlabel_style = {'size': 14, 'color': 'black'}
    gl.ylabel_style = {'size': 14, 'color': 'black'}
    gl.top_labels = False
    gl.right_labels = False

# HighResMIP panel 

def plot_highresmip_panel(ax, mode, model, experiment, member):

    datadir = f'/data/users/samuel.sancar/Heat_extremes/Trends/{mode}_model_trends/'
    fname_trend = f"{model}_{experiment}_{member}_tailwidth_trend.nc"
    cube = iris.load_cube(datadir + fname_trend)
    p87p5_trend = iris.load_cube(datadir + fname_trend.replace('tailwidth', 'p87p5'))

    masked_cube, _ = lat_land_mask(cube, LAT_MIN, LAT_MAX)
    p87p5_masked, _ = lat_land_mask(p87p5_trend, LAT_MIN, LAT_MAX)

    lon = cube.coord('longitude').points
    lat = cube.coord('latitude').points
    lon2d, lat2d = np.meshgrid(lon, lat)

    # p87p5 negative dots (white)
    neg_mask = p87p5_masked.data < 0
    neg_lats = lat2d[neg_mask]
    neg_lons = lon2d[neg_mask]

    # p-values < 0.05 (black stippling)
    sig_datadir = f'/data/users/samuel.sancar/Heat_extremes/Trends/Trend_significance/P-vals/'
    p_path = f"{sig_datadir}{model}_{experiment}_{member}_tailwidths_p_1958.nc"
    statsig_ds = xr.open_dataset(p_path)
    first_var = list(statsig_ds.data_vars.values())[0]
    sig_mask = (first_var.values < 0.05)

    ax.set_title(f"{model} {experiment} {member}", fontsize=16, weight='bold')
    ax.set_global()
    ax.set_extent([-180, 180, LAT_MIN, LAT_MAX], crs=ccrs.PlateCarree())
    plt.sca(ax)
    mappable = iplt.pcolormesh(masked_cube, cmap=CMAP, vmin=VMIN, vmax=VMAX)

    ax.scatter(neg_lons, neg_lats, color='white', s=0.9, alpha=0.6,
               transform=ccrs.PlateCarree(), label='p87p5 trend < 0')
    stipple_lats = lat2d[sig_mask]
    stipple_lons = lon2d[sig_mask]
    ax.scatter(stipple_lons, stipple_lats, color='black', s=0.1, alpha=0.8,
               transform=ccrs.PlateCarree(), label='p < 0.05')

    _add_common_features(ax)
    return mappable

# CANARI panel 

def plot_canari_panel(ax, member):
    datadir = f'/data/users/samuel.sancar/Heat_extremes/Trends/CANARI_trends/'
    fname_trend = f"{member}_tailwidth_trend.nc"
    cube = iris.load_cube(datadir + fname_trend)
    p87p5_trend = iris.load_cube(datadir + fname_trend.replace('tailwidth', 'p87p5'))

    # Mask
    masked_cube, _ = lat_land_mask(cube, LAT_MIN, LAT_MAX)
    p87p5_masked, _ = lat_land_mask(p87p5_trend, LAT_MIN, LAT_MAX)

    lon = cube.coord('longitude').points
    lat = cube.coord('latitude').points
    lon2d, lat2d = np.meshgrid(lon, lat)

    # p87p5 negative dots (white)
    neg_mask = p87p5_masked.data < 0
    neg_lats = lat2d[neg_mask]
    neg_lons = lon2d[neg_mask]

    # p-values < 0.05 (black stippling)
    sig_datadir = f'/data/users/samuel.sancar/Heat_extremes/Trends/Trend_significance/CANARI/P-vals/'
    p_path = f"{sig_datadir}{member}_tailwidths_p_1958.nc"
    statsig_ds = xr.open_dataset(p_path)
    first_var = list(statsig_ds.data_vars.values())[0]
    sig_mask = (first_var.values < 0.05)

    ax.set_title(f"CANARI member {member}", fontsize=18, weight='bold')
    ax.set_global()
    ax.set_extent([-180, 180, LAT_MIN, LAT_MAX], crs=ccrs.PlateCarree())
    plt.sca(ax)
    mappable = iplt.pcolormesh(masked_cube, cmap=CMAP, vmin=VMIN, vmax=VMAX)
    
    ax.scatter(neg_lons, neg_lats, color='white', s=0.5, alpha=0.6,
               transform=ccrs.PlateCarree(), label='p87p5 trend < 0')
    stipple_lats = lat2d[sig_mask]
    stipple_lons = lon2d[sig_mask]
    ax.scatter(stipple_lons, stipple_lats, color='black', s=0.1, alpha=0.8,
               transform=ccrs.PlateCarree(), label='p < 0.05')

    _add_common_features(ax)
    return mappable

def main():
    """
    Usage:
        python plot_tailwidth_4panel_pvals.py <Mode> <model> <experiment> <hr_member1> <hr_member2> <canari_member1> <canari_member2>

    Example:
        python plot_tailwidth_4panel_pvals.py Coupled HadGEM3-LL highres-future r1i1p1f1 r3i1p1f1 CANARI-01 CANARI-05
    """

    fig, axs = plt.subplots(
        2, 2,
        figsize=(22, 11),
        subplot_kw={'projection': ccrs.Robinson(central_longitude=CENTRAL_LON)}
    )

    fig.subplots_adjust(
        left=0.05, right=0.95, top=0.92, bottom=0.12,
        wspace=0.08, hspace=0.12
    )

    # Panels
    mappables = []
    mappables.append(plot_highresmip_panel(axs[0, 0], 'SST', 'MRI-AGCM3-2-H', 'highresSST-present', 'r1i1p1f1'))
    mappables.append(plot_highresmip_panel(axs[0, 1], 'Coupled', 'CNRM-CM6-1','hist-1950','r3i1p1f2'))
    mappables.append(plot_canari_panel(axs[1, 0], 12))
    mappables.append(plot_canari_panel(axs[1, 1], 28))

    # Colourbar strip
    cbar_ax = fig.add_axes([0.15, 0.05, 0.70, 0.025])
    cbar = fig.colorbar(mappables[0], cax=cbar_ax, orientation='horizontal')
    cbar.ax.tick_params(labelsize=16)
    cbar.set_label('Trend in Tailwidth (°C per decade)', fontsize=16)

    # Title and save
    #fig.suptitle("Trend in Tailwidth (1958–2022)\nHighResMIP (top) and CANARI (bottom)",
                 #fontsize=14, weight='bold', y=0.98)
    #plt.tight_layout()
    save_figure(plt, f'HighResMIP_CANARI_tailwidths_trends', 'Maps')
    plt.show()

if __name__ == "__main__":
    main()

