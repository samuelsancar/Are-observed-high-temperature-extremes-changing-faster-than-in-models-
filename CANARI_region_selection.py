import iris
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import iris.plot as iplt
import sys
from Lat_Land_mask_function import lat_land_mask

def pass_conditions(member1):
    
    member1_trend = iris.load_cube(f'/data/users/samuel.sancar/Heat_extremes/Trends/CANARI_trends/{member1}_tailwidth_trend.nc')
    member1_masked, _ = lat_land_mask(member1_trend, -60, 90)
    member1_agreement = np.load(f'/data/users/samuel.sancar/Heat_extremes/Trends/CANARI_trends/{member1}_trend_agreement.npy')
    member1_array = np.array(member1_masked.data)
    jra55_trend = iris.load_cube(f'/data/users/samuel.sancar/Heat_extremes/Trends/Reanalysis_trends/jra55__minmax_surf.015_tmax.reg_tl319._tailwidth_trends_1958-2023.nc')

    member1_sig = iris.load_cube(f'/data/users/samuel.sancar/Heat_extremes/Trends/Trend_significance/CANARI/P-vals/{member1}_tailwidths_p_1958.nc')

    member1_array = member1_array.astype(float)
    condition_1 = member1_sig.data < 0.05
    condition_2 = member1_array > 0
    condition_3 = member1_agreement >= 0.9

    # Combine conditions
    combined_conditions = condition_1 & condition_2 & condition_3
    masked_data = np.ma.masked_where(~combined_conditions, member1_masked.data)
    number_of_gridpoints = np.sum(masked_data)
    condition_count = condition_1.astype(int) + condition_2.astype(int) + condition_3.astype(int)
    
    return combined_conditions, number_of_gridpoints, member1_masked, masked_data, condition_count

variable = sys.argv[1]

if variable == 'iterate':
    member_passes = []

    for i in range(1, 41): 
        combined_conditions, number_of_gridpoints, member1_masked, masked_data = pass_conditions(i)
        print(f'Member {i} has {number_of_gridpoints} gridpoints that pass all 3 conditions.')
        member_passes.append(number_of_gridpoints)

    member_passes_array = np.array(member_passes)

    plt.figure(figsize=(12, 7))
    plt.hist(member_passes_array, bins=10, histtype='bar', edgecolor='black')
    plt.xlabel('Number of Gridpoints')
    plt.ylabel('Frequency')
    plt.title('Distribution of Gridpoints Passing All Conditions')
    plt.tight_layout()
    plt.show()
    
else:
    member1 = variable
    
    combined_conditions, number_of_gridpoints, member1_masked, masked_data, condition_count = pass_conditions(member1)
     
    print(f'Member {member1} has {number_of_gridpoints} gridpoints that pass all 3 conditions.')

    lats = member1_masked.coord('latitude').points
    lons = member1_masked.coord('longitude').points

    fig = plt.figure(figsize=(10, 6))
    ax = plt.axes(projection=ccrs.Robinson(central_longitude=0))
    ax.set_extent([-180, 180, -60, 90], crs=ccrs.PlateCarree(central_longitude=0))
    ax.add_feature(cfeature.COASTLINE)
    ax.add_feature(cfeature.BORDERS, linestyle=':')
    ax.set_global()

    mesh = ax.pcolormesh(lons, lats, masked_data, transform=ccrs.PlateCarree(central_longitude=0),
                        cmap='Reds', shading='auto')

    cb = plt.colorbar(mesh, orientation='horizontal', pad=0.05)
    cb.set_label('Tailwidth Trend (°C per decade)')

    plt.title(f'Regions passing conditions for CANARI member {member1}')
    plt.show()
    
    plt.figure(figsize=(10, 6))
    ax = plt.axes(projection=ccrs.Robinson(central_longitude= 0))
    ax.set_extent([-180, 180, -60, 90], crs=ccrs.PlateCarree(central_longitude=0))
    ax.coastlines(resolution='50m', linewidth=0.5)
    ax.add_feature(cfeature.NaturalEarthFeature(category='cultural', name='admin_0_boundary_lines_land', facecolor='none', scale='50m', edgecolor='0', lw=.05))
    ax.add_feature(cfeature.LAKES.with_scale('110m'), edgecolor='0', facecolor='none', lw=.05)
    ax.add_feature(cfeature.STATES.with_scale('110m'), edgecolor='0', facecolor='none', lw=.2)
    ax.add_feature(cfeature.RIVERS.with_scale('110m'), edgecolor='tab:blue', facecolor='none', lw=.4)
    ax.add_feature(cfeature.LAND, facecolor='lightgray')
    ax.add_feature(cfeature.OCEAN, facecolor='white')

    lon_grid, lat_grid = np.meshgrid(lons, lats)

    gl = ax.gridlines(draw_labels=True, linewidth=0.5, color='gray', alpha=0.5, linestyle='--')
    gl.top_labels = False
    gl.right_labels = False

    mesh = ax.pcolormesh(lon_grid, lat_grid, condition_count, cmap='viridis', transform=ccrs.PlateCarree(central_longitude=0))
    plt.colorbar(mesh, orientation='horizontal', label='Number of Conditions Met')
    plt.title('Conditions Met per Grid Point')
    plt.show()