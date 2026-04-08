import iris
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from matplotlib import cm
import matplotlib.colors as mcolors
from matplotlib.colors import ListedColormap, BoundaryNorm
from Lat_Land_mask_function import lat_land_mask
import iris.analysis
import sys
import glob
import os
import cf_units

coupled = sorted(glob.glob('/data/users/samuel.sancar/Heat_extremes/Trends/Coupled_model_trends/*tailwidth*'))
sst = sorted(glob.glob('/data/users/samuel.sancar/Heat_extremes/Trends/SST_model_trends/*tailwidth*'))

coupled_files = [(iris.load_cube(fname), fname) for fname in coupled]
sst_files = [(iris.load_cube(fname), fname) for fname in sst]

all_files = [("Coupled", cube, fname) for cube, fname in coupled_files] + \
            [("SST", cube, fname) for cube, fname in sst_files]

def pass_conditions(member1):
    sst_path = f'/data/users/samuel.sancar/Heat_extremes/Trends/SST_model_trends/{member1}_tailwidth_trend.nc'
    coupled_path = f'/data/users/samuel.sancar/Heat_extremes/Trends/Coupled_model_trends/{member1}_tailwidth_trend.nc'

    if os.path.exists(sst_path):
        member1_trend = iris.load_cube(sst_path)
    elif os.path.exists(coupled_path):
        member1_trend = iris.load_cube(coupled_path)
    else:
        raise FileNotFoundError(f"No trend file found for {member1}")

    basename1 = member1.split('_')[0]
    basename2 = member1.split('_')[2]
    basename = f'{basename1}_{basename2}'
    member1_agreement = np.load(f'/data/users/samuel.sancar/Heat_extremes/Trends/Trend_agreement/{basename}_agreement.npy')
    member1_trend_masked, _ = lat_land_mask(member1_trend, -60, 90)
    member1_array = member1_trend.data.astype(float)
    member1_sig = iris.load_cube(f'/data/users/samuel.sancar/Heat_extremes/Trends/Trend_significance/P-vals/{member1}_tailwidths_p_1958.nc')

    condition_1 = member1_sig.data < 0.05
    condition_2 = member1_array > 0
    condition_3 = member1_agreement >= 0.9
    
    if not member1_trend_masked.coord('latitude').has_bounds():
        member1_trend_masked.coord('latitude').guess_bounds()
    if not member1_trend_masked.coord('longitude').has_bounds():
        member1_trend_masked.coord('longitude').guess_bounds()
    
    # Fix units if missing or unknown
    if member1_trend_masked.coord('latitude').units.is_unknown():
        member1_trend_masked.coord('latitude').units = cf_units.Unit('degrees')
    if member1_trend_masked.coord('longitude').units.is_unknown():
        member1_trend_masked.coord('longitude').units = cf_units.Unit('degrees')

    area_weights = iris.analysis.cartography.area_weights(member1_trend_masked)
    
    combined_conditions = condition_1 & condition_2 & condition_3 
    masked_area = np.sum(area_weights[combined_conditions])
    number_of_gridpoints = np.sum(combined_conditions)
    condition_count = condition_1.astype(int) + condition_2.astype(int) + condition_3.astype(int)
    
    _0p9 = np.ma.masked_where(~condition_3, member1_trend_masked.data)
    no_0p9 = np.sum(_0p9)

    return masked_area, combined_conditions, number_of_gridpoints, member1_trend, condition_count, no_0p9

variable = sys.argv[1]

if variable == 'iterate':
    member_passes = []
    _90th_percentile = []
    for tag, cube, fname in all_files:
        member_name = os.path.basename(fname).split('_tailwidth')[0]
        passed_area, combined_conditions, number_of_gridpoints, member1_trend, condition_count, no_0p9 = pass_conditions(member_name)
        print(f'Member {member_name} ({tag}) has {passed_area}metres squared passing all conditions.')
        if member_name == "EC-Earth3P-HR_highresSST-present_r3i1p1f1":
            chosen_area = passed_area
        member_passes.append(passed_area)
        _90th_percentile.append(no_0p9)
        
    era5 = iris.load_cube('/data/users/samuel.sancar/Heat_extremes/Trends/Reanalysis_trends/era5_daily_maximum_2m_temperature_tailwidth_trends_1958-2023.nc')
    era5_masked, _ = lat_land_mask(era5, -60, 90)
    era5_sig = iris.load_cube('/data/users/samuel.sancar/Heat_extremes/Trends/Reanalysis_trends/era5_daily_maximum_2m_temperature_tailwidth_p_1958-2023.nc')
    era5_agreement = np.load('/home/users/samuel.sancar/Documents/Heat_Extremes/Model_agreement/ERA5_HighResMIP_agreement.npy')
    era5_array = era5_masked.data.astype(float)

    era5_condition_1 =  era5_sig.data < 0.05
    era5_condition_2 =  era5_array > 0 
    era5_condition_3 =  era5_agreement >= 0.9
    
    if not era5_masked.coord('latitude').has_bounds():
        era5_masked.coord('latitude').guess_bounds()
    if not era5_masked.coord('longitude').has_bounds():
        era5_masked.coord('longitude').guess_bounds()
    
    # Fix units if missing or unknown
    if era5_masked.coord('latitude').units.is_unknown():
        era5_masked.coord('latitude').units = cf_units.Unit('degrees')
    if era5_masked.coord('longitude').units.is_unknown():
        era5_masked.coord('longitude').units = cf_units.Unit('degrees')

    era5_area_weights = iris.analysis.cartography.area_weights(era5_masked)
    
    era5_combined_conditions = era5_condition_1 & era5_condition_2 & era5_condition_3
    era5_masked_area = np.sum(era5_area_weights[era5_combined_conditions])

    era5_number_of_gridpoints = np.sum(era5_combined_conditions)
    
    era5_0p9 = np.ma.masked_where(~era5_condition_3, era5_masked.data)
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
    plt.axvline(x=chosen_area, color="#2644D7", linestyle='--', linewidth=2.5, label='EC-Earth3P-HR_\nr1i1p1f1 (AMIP) Value')
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
    
    passed_area, combined_conditions, number_of_gridpoints, member1_trend, condition_count, _ = pass_conditions(member1)
     
    print(f'Member {member1} has {number_of_gridpoints} gridpoints that pass all 3 conditions.')

    lats = member1_trend.coord('latitude').points
    lons = member1_trend.coord('longitude').points

    fig = plt.figure(figsize=(10, 6))
    ax = plt.axes(projection=ccrs.Robinson(central_longitude=0))
    ax.set_extent([-180, 180, -60, 90], crs=ccrs.PlateCarree(central_longitude=0))
    ax.add_feature(cfeature.COASTLINE)
    ax.add_feature(cfeature.BORDERS, linestyle=':')
    ax.set_global()

    mesh = ax.pcolormesh(lons, lats, transform=ccrs.PlateCarree(central_longitude=0),cmap='RdYlBu', shading='auto')

    cb = plt.colorbar(mesh, orientation='horizontal', pad=0.05)
    cb.set_label('Tailwidth Trend (°C per decade)')

    plt.title(f'Regions passing conditions for HighResMIP member\n {member1}')
    plt.show()
    

    plt.figure(figsize=(10, 6))
    ax = plt.axes(projection=ccrs.Robinson(central_longitude=0))
    ax.set_extent([-180, 180, -60, 90], crs=ccrs.PlateCarree(central_longitude=0))
    ax.coastlines()
    ax.add_feature(cfeature.BORDERS, linestyle=':')

    lon_grid, lat_grid = np.meshgrid(lons, lats)

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

    plt.title(f'HighResMIP member\n {member1}')
    plt.show()
