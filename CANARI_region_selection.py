import iris
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import iris.plot as iplt
from matplotlib import cm
import matplotlib.colors as mcolors
from matplotlib.colors import ListedColormap, BoundaryNorm
import cf_units
import sys
from Lat_Land_mask_function import lat_land_mask

def pass_conditions(member1):
    
    member1_trend = iris.load_cube(f'/data/users/samuel.sancar/Heat_extremes/Trends/CANARI_trends/{member1}_tailwidth_trend.nc')
    member1_masked, _ = lat_land_mask(member1_trend, -60, 90)
    member1_agreement = np.load(f'/data/users/samuel.sancar/Heat_extremes/Trends/CANARI_trends/{member1}_trend_agreement.npy')
    member1_array = np.array(member1_masked.data)

    member1_sig = iris.load_cube(f'/data/users/samuel.sancar/Heat_extremes/Trends/Trend_significance/CANARI/P-vals/{member1}_tailwidths_p_1958.nc')

    member1_array = member1_trend.data.astype(float)
    member1_sig = iris.load_cube(f'/data/users/samuel.sancar/Heat_extremes/Trends/Trend_significance/CANARI/P-vals/{member1}_tailwidths_p_1958.nc')

    condition_1 = member1_sig.data < 0.05
    condition_2 = member1_array > 0
    condition_3 = member1_agreement >= 0.9
    
    if not member1_trend.coord('latitude').has_bounds():
        member1_trend.coord('latitude').guess_bounds()
    if not member1_trend.coord('longitude').has_bounds():
        member1_trend.coord('longitude').guess_bounds()
    
    # Fix units if missing or unknown
    if member1_trend.coord('latitude').units.is_unknown():
        member1_trend.coord('latitude').units = cf_units.Unit('degrees')
    if member1_trend.coord('longitude').units.is_unknown():
        member1_trend.coord('longitude').units = cf_units.Unit('degrees')

    area_weights = iris.analysis.cartography.area_weights(member1_trend)
    
    combined_conditions = condition_1 & condition_2 & condition_3
    trend_copy = member1_trend.copy()
    trend_copy.data = np.ma.masked_where(~combined_conditions, member1_trend.data)
    plot = np.ma.masked_where(~combined_conditions, member1_trend.data)
    masked_area = trend_copy.collapsed(['latitude', 'longitude'], iris.analysis.SUM, weights=area_weights).data
    number_of_gridpoints = np.sum(combined_conditions)
    condition_count = condition_1.astype(int) + condition_2.astype(int) + condition_3.astype(int)
    
    _0p9 = np.ma.masked_where(~condition_3, member1_trend.data)
    no_0p9 = np.sum(_0p9)

    return member1_masked, masked_area, combined_conditions, number_of_gridpoints, member1_trend, trend_copy, condition_count, no_0p9, plot

variable = sys.argv[1]

if variable == 'iterate':
    member_passes = []
    _90th_percentile = []
    for i in range(1, 41): 
        member1_masked, passed_area, combined_conditions, number_of_gridpoints, member1_trend, trend_copy, condition_count, no_0p9, plot = pass_conditions(i)
        print(f'Member {i} has {number_of_gridpoints} gridpoints that pass all 3 conditions.')
        #member_passes.append(number_of_gridpoints)
        if i == 32:
            chosen_area = passed_area
        member_passes.append(passed_area)
        _90th_percentile.append(no_0p9)

    member_passes_array = np.array(member_passes)
    era5 = iris.load_cube('/data/users/samuel.sancar/Heat_extremes/Trends/Reanalysis_trends/era5_daily_maximum_2m_temperature_tailwidth_trends_1958-2023.nc')
    era5_sig = iris.load_cube('/data/users/samuel.sancar/Heat_extremes/Trends/Reanalysis_trends/era5_daily_maximum_2m_temperature_tailwidth_p_1958-2023.nc')
    era5_agreement = np.load('/home/users/samuel.sancar/Documents/Heat_Extremes/Model_agreement/CANARI_ERA5_gridpoint_fractions.npy')
    era5_array = era5.data.astype(float)

    era5_condition_1 =  era5_sig.data < 0.05
    era5_condition_2 =  era5_array > 0 
    era5_condition_3 =  era5_agreement >= 0.9
    
    if not era5.coord('latitude').has_bounds():
        era5.coord('latitude').guess_bounds()
    if not era5.coord('longitude').has_bounds():
        era5.coord('longitude').guess_bounds()
    
    # Fix units if missing or unknown
    if era5.coord('latitude').units.is_unknown():
        era5.coord('latitude').units = cf_units.Unit('degrees')
    if era5.coord('longitude').units.is_unknown():
        era5.coord('longitude').units = cf_units.Unit('degrees')

    era5_area_weights = iris.analysis.cartography.area_weights(era5)
    
    era5_combined_conditions = era5_condition_1 & era5_condition_2 & era5_condition_3
    era5_copy = era5.copy()
    era5_copy.data = np.ma.masked_where(~era5_combined_conditions, era5.data)
    era5_masked_area = era5_copy.collapsed(['latitude', 'longitude'], iris.analysis.SUM, weights=era5_area_weights).data
      
    era5_masked_data = np.ma.masked_where(~era5_combined_conditions, era5.data)
    era5_number_of_gridpoints = np.sum(era5_combined_conditions)
    
    era5_0p9 = np.ma.masked_where(~era5_condition_3, era5.data)
    era5_no_0p9 = np.sum(era5_0p9)
    
    members_array = np.array(member_passes)
    era5_percentile = (np.sum(members_array < era5_masked_area))/ (len(members_array)) * 100
    
    print(f'ERA5 has {era5_masked_area} metres squared that passed all three conditions. It is in the {np.round(era5_percentile, 1)} % of the model members spread for number of grid points that passed all three conditions.')

    #_90th_array = np.array(_90th_percentile)
    #era5_90th_percentile = (np.sum(_90th_array < era5_no_0p9)) / (len(_90th_array)) * 100
    #print(f'ERA5 has {era5_90th_percentile} gridpoints where it has trends greater than the 90% of models. It is in the {np.round(era5_90th_percentile)} % of the model members spread for number of grid points that passed all three conditions.')
     
    plt.style.use('ggplot')
    plt.figure(figsize=(12, 7))
    plt.hist(member_passes, bins='auto', edgecolor='black')
    plt.axvline(x=era5_masked_area, color="#35D726", linestyle='--', linewidth=2.5, label='ERA5 Value')
    plt.axvline(x=chosen_area, color="#2644D7", linestyle='--', linewidth=2.5, label='Member 32')
    plt.xlabel('Total area to pass all 3 conditions (m2)')
    plt.ylabel('Frequency')
    plt.title('Distribution of Gridpoints Passing All Conditions')
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.show()
    
else:
    member1 = variable
    bounds = [0,1,2,3,4]
    cmap = cm.get_cmap('viridis', len(bounds)-1)
    norm = mcolors.BoundaryNorm(bounds, cmap.N)
    
    member1_masked, passed_area, combined_conditions, number_of_gridpoints, member1_trend, trend_copy, condition_count, no_0p9, plot = pass_conditions(member1)
     
    print(f'Member {member1} has {number_of_gridpoints} gridpoints that pass all 3 conditions.')

    lats = member1_masked.coord('latitude').points
    lons = member1_masked.coord('longitude').points

    #fig = plt.figure(figsize=(10, 6))
    #ax = plt.axes(projection=ccrs.Robinson(central_longitude=0))
    #ax.set_extent([-180, 180, -60, 90], crs=ccrs.PlateCarree(central_longitude=0))
    #ax.add_feature(cfeature.COASTLINE)
    #ax.add_feature(cfeature.BORDERS, linestyle=':')
    #ax.set_global()

    #mesh = ax.pcolormesh(lons, lats, masked_data, transform=ccrs.PlateCarree(central_longitude=0),
                        #cmap='Reds', shading='auto')

    #cb = plt.colorbar(mesh, orientation='horizontal', pad=0.05)
    #cb.set_label('Tailwidth Trend (°C per decade)')

    #plt.title(f'Regions passing conditions for CANARI member {member1}')
    #plt.show()
    
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

    mesh = ax.pcolormesh(lon_grid, lat_grid, condition_count,
                        cmap=cmap, norm=norm,
                        transform=ccrs.PlateCarree(central_longitude=0))

    cb = plt.colorbar(mesh, orientation='horizontal', pad=0.05, aspect=50)
    cb.set_label('Number of conditions true')

    # Compute midpoints for tick positions
    bounds = [0, 1, 2, 3, 4] 
    midpoints = [(bounds[i] + bounds[i+1]) / 2 for i in range(len(bounds)-1)]

    cb.set_ticks(midpoints)
    cb.set_ticklabels(['0', '1', '2', '3'])
    plt.title('Large Ensemble member 32')
    plt.show()