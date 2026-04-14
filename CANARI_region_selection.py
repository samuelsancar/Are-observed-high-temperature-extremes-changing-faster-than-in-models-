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
from Save_figure_function import save_figure
from Lat_Land_mask_function import lat_land_mask

def guess_bounds(cube):
    
    if not cube.coord('latitude').has_bounds():
        cube.coord('latitude').guess_bounds()
    if not cube.coord('longitude').has_bounds():
        cube.coord('longitude').guess_bounds()
    
    # Fix units if missing or unknown
    if cube.coord('latitude').units.is_unknown():
        cube.coord('latitude').units = cf_units.Unit('degrees')
    if cube.coord('longitude').units.is_unknown():
        cube.coord('longitude').units = cf_units.Unit('degrees')
        
    return cube
    
def pass_conditions(member1):
    
    member1_trend = iris.load_cube(f'/data/users/samuel.sancar/Heat_extremes/Trends/CANARI_trends/{member1}_tailwidth_trend_1958-2023.nc')
    member1_masked, _ = lat_land_mask(member1_trend, -60, 80)
    member1_agreement = np.load(f'/data/users/samuel.sancar/Heat_extremes/Trends/CANARI_trends/{member1}_trend_agreement.npy')
    member1_trend_1980 = iris.load_cube(f'/data/users/samuel.sancar/Heat_extremes/Trends/CANARI_trends/{member1}_tailwidth_trend_1980-2023.nc')
    member1_1980_masked, _ = lat_land_mask(member1_trend_1980, -60, 80)
    
    #member1_sig = iris.load_cube(f'/data/users/samuel.sancar/Heat_extremes/Trends/Trend_significance/CANARI/FDR_corrections_All/{member1}_tailwidths_p_1958_fdr_significant_all.nc')
    member1_sig = iris.load_cube(f'/data/users/samuel.sancar/Heat_extremes/Trends/Trend_significance/CANARI/P-vals/{member1}_tailwidths_p_1958.nc')
    
    member1_masked = guess_bounds(member1_masked)
    
    # Ensure a "valid data" mask from the primary cube
    valid = ~np.ma.getmaskarray(member1_masked.data)

    # 4 Objective Conditions
    #condition_1 = member1_sig.data.astype(bool)
    condition_1 = np.ma.filled(member1_sig.data < 0.05, False)
    condition_2 = np.ma.filled(member1_masked.data > 0, False)
    condition_3 = np.ma.filled(member1_agreement >= 0.9, False)
    condition_4 = np.ma.filled(member1_1980_masked.data > 0, False)

    combined_conditions = valid & condition_1 & condition_2 & condition_3 & condition_4

    # Count grid points
    number_of_gridpoints = int(np.count_nonzero(combined_conditions))

    # Area weighting
    aw = np.ma.array(
        iris.analysis.cartography.area_weights(member1_masked),
        mask=~combined_conditions
    )
    masked_area = float(aw.sum()) 
    
    number_of_gridpoints = np.sum(combined_conditions)
    condition_count = condition_1.astype(int) + condition_2.astype(int) + condition_3.astype(int) + condition_4.astype(int)
    
    _0p9 = np.ma.masked_where(~condition_3, member1_masked.data)
    no_0p9 = np.sum(_0p9)

    return member1_masked, masked_area, combined_conditions, number_of_gridpoints, condition_count, no_0p9

variable = sys.argv[1]

if variable == 'iterate':
    member_passes = []
    _80th_percentile = []
    for i in range(1, 41): 
        
        member1_masked, passed_area, combined_conditions, number_of_gridpoints, condition_count, no_0p9 = pass_conditions(i)
        print(f'Member {i} has {number_of_gridpoints} gridpoints that pass all 3 conditions.')
        #member_passes.append(number_of_gridpoints)
        
        if i == 29:
            chosen_area = passed_area
        member_passes.append(passed_area)
        _80th_percentile.append(no_0p9)

    member_passes_array = np.array(member_passes)
    era5 = iris.load_cube('/data/users/samuel.sancar/Heat_extremes/Trends/Reanalysis_trends/era5_daily_maximum_2m_temperature_tailwidth_trends_1958-2023.nc')
    era5_masked, _ = lat_land_mask(era5, -60, 80)       
    #era5_sig = iris.load_cube('/data/users/samuel.sancar/Heat_extremes/Trends/Reanalysis_trends/era5_daily_maximum_2m_temperature_tailwidth_fdr_all_1958-2023.nc')
    era5_sig = iris.load_cube('/data/users/samuel.sancar/Heat_extremes/Trends/Reanalysis_trends/era5_daily_maximum_2m_temperature_tailwidth_p_1958-2023.nc')

    era5_agreement = np.load('/home/users/samuel.sancar/Documents/Temperature_extremes/Model_agreement/CANARI_ERA5_gridpoint_fractions.npy')
    era5_1980 = iris.load_cube('/data/users/samuel.sancar/Heat_extremes/Trends/Reanalysis_trends/era5_daily_maximum_2m_temperature_tailwidth_trends_1980-2023.nc')
    era5_1980_masked, _ = lat_land_mask(era5_1980, -60, 80)
    era5_1980_array = era5_1980_masked.data.astype(float)
    era5_array = era5_masked.data.astype(float)
    
    merra = iris.load_cube('/data/users/samuel.sancar/Heat_extremes/Trends/Reanalysis_trends/MERRA_daily_maximum_2m_temperature_tailwidth_trends_1980-2023.nc')
    jra = iris.load_cube('/data/users/samuel.sancar/Heat_extremes/Trends/Reanalysis_trends/jra55__minmax_surf.015_tmax.reg_tl319._tailwidth_trends_1958-2023.nc')
    jra_1980 = iris.load_cube('/data/users/samuel.sancar/Heat_extremes/Trends/Reanalysis_trends/jra55__minmax_surf.015_tmax.reg_tl319._tailwidth_trends_1980-2023.nc')
    merra, _= lat_land_mask(merra, -60, 80)
    jra, _ = lat_land_mask(jra, -60, 80)
    jra_1980, _ = lat_land_mask(jra_1980, -60, 80)
    
    
    era5_masked = guess_bounds(era5_masked)
    merra_masked = guess_bounds(merra)
    jra_masked = guess_bounds(jra)
    jra_1980_masked = guess_bounds(jra_1980)
    
    # Ensure a "valid data" mask from the primary cube
    era5_valid = ~np.ma.getmaskarray(era5_masked.data)

    # 4 Objective Conditions
    #era5_condition_1 = era5_sig.data.astype(bool)
    era5_condition_1 = np.ma.filled(era5_sig.data < 0.05, False)
    era5_condition_2 = np.ma.filled(era5_masked.data > 0, False) #& np.ma.filled(merra.data > 0, False) & np.ma.filled(jra.data > 0, False) & np.ma.filled(jra_1980.data > 0, False)
    era5_condition_3 = np.ma.filled(era5_agreement >= 0.9, False)
    era5_condition_4 = np.ma.filled(era5_1980_masked.data > 0, False)

    era5_combined_conditions = era5_valid & era5_condition_1 & era5_condition_2 & era5_condition_3 & era5_condition_4

    # Count grid points
    number_of_gridpoints = int(np.count_nonzero(era5_combined_conditions))

    # Area weighting
    era5_aw = np.ma.array(
        iris.analysis.cartography.area_weights(era5_masked),
        mask=~era5_combined_conditions
    )
    
    era5_masked_area = float(era5_aw.sum()) 
    
    era5_number_of_gridpoints = np.sum(era5_combined_conditions)
    
    era5_0p9 = np.ma.masked_where(~era5_condition_3, era5_masked.data)
    era5_no_0p9 = np.sum(era5_0p9)
    
    members_array = np.array(member_passes)
    era5_percentile = (np.sum(members_array < era5_masked_area))/ (len(members_array)) * 100
    member_percentile = (np.sum(members_array < chosen_area))/ (len(members_array)) * 100
    
    print(f'ERA5 has {era5_masked_area} metres squared, {era5_number_of_gridpoints} gridpoints, that passed all three conditions. It is in the {np.round(era5_percentile, 1)} % of the model members spread for number of grid points that passed all three conditions.')
    print(f'Member 29 is in the {np.round(member_percentile, 1)} % of the model members spread for number of grid points that passed all three conditions.')
    #_80th_array = np.array(_80th_percentile)
    #era5_80th_percentile = (np.sum(_80th_array < era5_no_0p9)) / (len(_80th_array)) * 100
    #print(f'ERA5 has {era5_80th_percentile} gridpoints where it has trends greater than the 80% of models. It is in the {np.round(era5_80th_percentile)} % of the model members spread for number of grid points that passed all three conditions.')
     
    #plt.style.use('ggplot')
    plt.figure(figsize=(12, 6))
    plt.hist(member_passes, bins = 'auto', histtype='bar',
    color="#4682B4",  
    alpha=0.85,
    edgecolor='white',
    linewidth=2
    , label='CANARI LE members')
    
    ax = plt.gca()
    for spine in ax.spines.values():
        spine.set_linewidth(2)

    plt.axvline(x=era5_masked_area, color='black', linestyle='-', linewidth=2.5, label='ERA5')
    plt.axvline(x=chosen_area, color='#d55e00', linestyle='-', linewidth=2.5, label='Member 29')
    plt.xlabel('Total area to pass all 4 conditions (\u00D710\u00B9\u00B2 m\u00B2)', fontsize=21)
    plt.ylabel('Frequency', fontsize=21)
    plt.xticks(fontsize=19)
    plt.yticks(fontsize=19)
    plt.title('Distribution of area passing all conditions', fontsize=22, fontweight='bold')
    plt.grid(True, linestyle='--', alpha=0.5, linewidth=1.2)
    plt.legend(fontsize=19)
    plt.tight_layout()
    save_figure(plt, f'CANARI_Histogram_region_selection', 'Graphs')
    plt.show()
    
# Plot individual member maps
else:
    member1 = variable
    bounds = [0,1,2,3,4,5]
    cmap = cm.get_cmap('viridis', len(bounds)-1)
    norm = mcolors.BoundaryNorm(bounds, cmap.N)
    
    if member1 =='era5':
        
        era5 = iris.load_cube('/data/users/samuel.sancar/Heat_extremes/Trends/Reanalysis_trends/era5_daily_maximum_2m_temperature_tailwidth_trends_1958-2023.nc')
        era5_masked, _ = lat_land_mask(era5, -60, 80)       
        #era5_sig = iris.load_cube('/data/users/samuel.sancar/Heat_extremes/Trends/Reanalysis_trends/era5_daily_maximum_2m_temperature_tailwidth_fdr_all_1958-2023.nc')
        era5_sig = iris.load_cube('/data/users/samuel.sancar/Heat_extremes/Trends/Reanalysis_trends/era5_daily_maximum_2m_temperature_tailwidth_p_1958-2023.nc')

        era5_agreement = np.load('/home/users/samuel.sancar/Documents/Temperature_extremes/Model_agreement/CANARI_ERA5_gridpoint_fractions.npy')
        era5_1980 = iris.load_cube('/data/users/samuel.sancar/Heat_extremes/Trends/Reanalysis_trends/era5_daily_maximum_2m_temperature_tailwidth_trends_1980-2023.nc')
        era5_1980_masked, _ = lat_land_mask(era5_1980, -60, 80)
        era5_masked = guess_bounds(era5_masked)
        
        # Ensure a "valid data" mask from the primary cube
        era5_valid = ~np.ma.getmaskarray(era5_masked.data)

        # 4 Objective Conditions
        #era5_condition_1 = era5_sig.data.astype(bool)
        era5_condition_1 = np.ma.filled(era5_sig.data < 0.05, False)
        era5_condition_2 = np.ma.filled(era5_masked.data > 0, False) #& np.ma.filled(merra.data > 0, False) & np.ma.filled(jra.data > 0, False) & np.ma.filled(jra_1980.data > 0, False)
        era5_condition_3 = np.ma.filled(era5_agreement >= 0.9, False)
        era5_condition_4 = np.ma.filled(era5_1980_masked.data > 0, False)

        era5_combined_conditions = era5_valid & era5_condition_1 & era5_condition_2 & era5_condition_3 & era5_condition_4

        condition_count = era5_condition_1.astype(int) + era5_condition_2.astype(int) + era5_condition_3.astype(int) + era5_condition_4.astype(int)
        
        lats = era5_masked.coord('latitude').points
        lons = era5_masked.coord('longitude').points
        
        plt.figure(figsize=(12,6))
        #ax = plt.axes(projection=ccrs.Robinson(central_longitude= 0))
        ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=0))
        ax.set_extent([-180, 180, -60, 80], crs=ccrs.PlateCarree(central_longitude=0))
        
        for spine in ax.spines.values():
            spine.set_linewidth(3)   
            spine.set_edgecolor('black')
            
        ax.coastlines()
        ax.add_feature(cfeature.LAND, facecolor='lightgray')
        ax.add_feature(cfeature.OCEAN, facecolor='white', zorder = 3)

        lon_grid, lat_grid = np.meshgrid(lons, lats)

        mesh = ax.pcolormesh(lon_grid, lat_grid, condition_count,
                            cmap=cmap, norm=norm,
                            transform=ccrs.PlateCarree(central_longitude=0))

        #cb = plt.colorbar(mesh, orientation='horizontal', pad=0.05, aspect=50)
        #cb.set_label('Number of conditions true', fontsize=13)

        # Compute midpoints for tick positions
        #midpoints = [(bounds[i] + bounds[i+1]) / 2 for i in range(len(bounds)-1)]

        #cb.set_ticks(midpoints)
        #cb.ax.tick_params(labelsize=12)
        #cb.set_ticklabels(['0', '1', '2', '3', '4'])
        plt.title(f'ERA5', fontsize=18, fontweight='bold')
        save_figure(plt, f'ERA5_region_selection', 'Maps')
        plt.show()
        
    else:
        member1_masked, passed_area, combined_conditions, number_of_gridpoints, condition_count, no_0p9 = pass_conditions(member1)
        
        print(f'Member {member1} has {number_of_gridpoints} gridpoints that pass all 3 conditions.')

        lats = member1_masked.coord('latitude').points
        lons = member1_masked.coord('longitude').points
        
        plt.figure(figsize=(12,6))
        #ax = plt.axes(projection=ccrs.Robinson(central_longitude= 0))
        ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=0))
        
        for spine in ax.spines.values():
            spine.set_linewidth(3)  
            spine.set_edgecolor('black')
    
        ax.set_extent([-180, 180, -60, 80], crs=ccrs.PlateCarree(central_longitude=0))
        ax.coastlines()
        ax.add_feature(cfeature.LAND, facecolor='lightgray')
        ax.add_feature(cfeature.OCEAN, facecolor='white', zorder = 3)

        lon_grid, lat_grid = np.meshgrid(lons, lats)

        mesh = ax.pcolormesh(lon_grid, lat_grid, condition_count,
                            cmap=cmap, norm=norm,
                            transform=ccrs.PlateCarree(central_longitude=0))

        cb = plt.colorbar(mesh, orientation='horizontal', pad=0.05, aspect=50)
        cb.set_label('Number of conditions true', fontsize=16)

        # Compute midpoints for tick positions
        midpoints = [(bounds[i] + bounds[i+1]) / 2 for i in range(len(bounds)-1)]

        cb.set_ticks(midpoints)
        cb.ax.tick_params(labelsize=15)
        cb.set_ticklabels(['0', '1', '2', '3', '4'])
        plt.title(f'CANARI-LE member {member1}', fontsize=18, fontweight='bold')
        save_figure(plt, f'CANARI_member_{member1}_region_selection', 'Maps')
        plt.show()