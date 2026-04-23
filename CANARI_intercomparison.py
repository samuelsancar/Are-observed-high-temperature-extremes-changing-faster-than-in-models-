"""
13/11/2025
----------
CANARI intercomparison: Loops through each member, and calculates at each grid point, 
what percentile in the ensemble spread said grid point lies.
"""

import iris
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import glob
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.colors as mcolors
from Lat_Land_mask_function import lat_land_mask
from Save_figure_function import save_figure

member = sys.argv[1]

# Beginning & end year+1 
t1 = 1958
t2 = 2023

lat_minimum = -60
lat_maximum = 90

template = iris.load_cube("template_cube_1deg.nc")

# Colour bar 
cmap = mcolors.ListedColormap([
    "#5AB4AC",
    "#99D594",
    "#E6F5B0",
    "#F6B661",
    "#F87234",
])
cmap.set_under("#544B94")  # below first bound
cmap.set_over("#8D1B15DB")   # above last bound

bounds = [0.1, 5, 25, 75, 95, 100]
norm = mcolors.BoundaryNorm(bounds, cmap.N, clip=False)

datadir = "/data/users/samuel.sancar/Heat_extremes/Trends/Reanalysis_trends/"
era5 = iris.load_cube(f"{datadir}era5_daily_maximum_2m_temperature_tailwidth_trends_{t1}-{t2-1}.nc")
era5_masked_cube, _ = lat_land_mask(era5, template, lat_minimum, lat_maximum)
era5_masked = era5_masked_cube.data 

CANARI_filenames = sorted(glob.glob(f"/data/users/samuel.sancar/Heat_extremes/Trends/CANARI_trends/*_tailwidth_trend_{t1}-{t2-1}.nc"))

CANARI_files = [(iris.load_cube(fname), fname) for fname in CANARI_filenames]

CANARI_masked_cubes = []
CANARI_member_names = []

for cube, fname in CANARI_files:
    cube_m, _ = lat_land_mask(cube, template, lat_minimum, lat_maximum)
    CANARI_masked_cubes.append(cube_m)

    filename = os.path.splitext(os.path.basename(fname))[0]
    CANARI_member_names.append(filename.split("_")[0])  

# Lat/lon for plotting
lons = template.coord("longitude").points
lats = template.coord("latitude").points
lon2d, lat2d = np.meshgrid(lons, lats)

# ERA5 branch: percentile of ERA5 relative to CANARI
if member.lower() == "era5":

    # Load p87p5 trend (mask it for stippling)
    p87p5_trend = iris.load_cube(f"{datadir}era5_daily_maximum_2m_temperature_p87p5_trends_{t1}-{t2-1}.nc")
    
    p87p5_masked_cube, _ = lat_land_mask(p87p5_trend, template, lat_minimum, lat_maximum)
    p87p5_masked = p87p5_masked_cube.data

    neg_indices = np.where(np.ma.filled(p87p5_masked < 0, False))
    neg_lats = lat2d[neg_indices]
    neg_lons = lon2d[neg_indices]

    # Model stack
    # shape: (n_members + 1, lat, lon)
    model_stack = np.ma.stack([c.data for c in CANARI_masked_cubes], axis=0)

    # Percent of stack <= ERA5 
    agreement = np.mean(np.ma.filled(model_stack <= era5_masked, False), axis=0)
    agreement = np.ma.array(agreement, mask=np.ma.getmaskarray(era5_masked))  # keep domain mask

    # Plot map
    fig = plt.figure(figsize=(12, 6))
    ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=0))
    ax.set_extent([-180, 180, lat_minimum, lat_maximum], crs=ccrs.PlateCarree())

    for spine in ax.spines.values():
        spine.set_linewidth(2)
        spine.set_edgecolor("black")

    c = ax.pcolormesh(
        lon2d, lat2d, agreement * 100,
        cmap=cmap, norm=norm,
        shading="auto",
        transform=ccrs.PlateCarree()
    )

    ax.scatter(
        neg_lons, neg_lats,
        color="white", s=0.9, alpha=0.9,
        transform=ccrs.PlateCarree(),
        label="p87p5 trend < 0"
    )

    ax.coastlines(resolution="50m", linewidth=0.8)
    ax.add_feature(cfeature.LAND, facecolor="lightgray")
    ax.add_feature(cfeature.OCEAN, facecolor="aliceblue")

    #cbar = plt.colorbar(
    #    c,
    #    ax=ax,
    #    orientation='horizontal',
    #    pad=0.05,
    #    aspect=50,
    #    extend='both')
    #cbar.set_ticklabels(['0', '5', '25', '75', '95', '100'])
    #cbar.set_label(f'Percent of CANARI-LE ensemble member trends less than ERA5 [%]', fontsize=11)
    
    plt.title(f'ERA5', fontsize=22, weight='bold')
    plt.tight_layout()
    save_figure(plt, "ERA5_percentile_map", "Maps")
    plt.show()

# CANARI member branch: percentile of a chosen member relative to (other CANARI + ERA5)
else:

    # loop through each member, but only plot the one asked for.
    for i, (canari_cube_masked, name) in enumerate(zip(CANARI_masked_cubes, CANARI_member_names)):

        model_name = name
        if model_name != member:
            continue

        # Load p87p5 for this member and mask it for stippling
        canari_p87p5 = iris.load_cube(
            f"/data/users/samuel.sancar/Heat_extremes/Trends/CANARI_trends/{model_name}_p87p5_trend_{t1}-{t2-1}.nc"
        )
        p87p5_masked_cube, _ = lat_land_mask(canari_p87p5, template, lat_minimum, lat_maximum)
        p87p5_masked = p87p5_masked_cube.data

        neg_indices = np.where(np.ma.filled(p87p5_masked < 0, False))
        neg_lats = lat2d[neg_indices]
        neg_lons = lon2d[neg_indices]

        # Model stack (keep ERA5 in it, as requested) 
        # stack = all other CANARI members + ERA5
        others = [CANARI_masked_cubes[j].data for j in range(len(CANARI_masked_cubes)) if j != i]
        model_stack = np.ma.stack(others + [era5_masked], axis=0)

        target = canari_cube_masked.data

        # Percent of stack <= target
        agreement = np.mean(np.ma.filled(model_stack <= target, False), axis=0)
        agreement = np.ma.array(agreement, mask=np.ma.getmaskarray(target))  # preserve domain

        # Plot
        fig = plt.figure(figsize=(12, 6))
        ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=0))
        ax.set_extent([-180, 180, lat_minimum, lat_maximum], crs=ccrs.PlateCarree())

        for spine in ax.spines.values():
            spine.set_linewidth(2)
            spine.set_edgecolor("black")

        c = ax.pcolormesh(
            lon2d, lat2d, agreement * 100,
            cmap=cmap, norm=norm,
            shading="auto",
            transform=ccrs.PlateCarree()
        )

        ax.scatter(
            neg_lons, neg_lats,
            color="white", s=0.9, alpha=0.6,
            transform=ccrs.PlateCarree(),
            label="p87p5 trend < 0"
        )

        ax.coastlines(resolution="50m", linewidth=0.8)
        ax.add_feature(cfeature.LAND, facecolor="lightgray")
        ax.add_feature(cfeature.OCEAN, facecolor="aliceblue")

                
        cbar = plt.colorbar(c,ax=ax,orientation='horizontal',pad=0.05,aspect=50,extend='both')

        cbar.set_ticklabels(['0', '5', '25', '75', '95', '100'])
        cbar.ax.tick_params(labelsize=16)
        cbar.set_label(f'Percent of CANARI-LE trends less than member [%]', fontsize=18)
        
        plt.title(f'CANARI-LE member {model_name}', fontsize=22, weight='bold')
        plt.tight_layout()
        save_figure(plt, f'{model_name}_percentile_map', 'Maps')
        plt.show()



        

