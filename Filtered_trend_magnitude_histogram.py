import iris
import iris.analysis
import cf_units
import sys
import matplotlib.pyplot as plt
import numpy as np
import glob
from Lat_Land_mask_function import lat_land_mask
from Save_figure_function import save_figure

# Command-line arguments
t1 = int(sys.argv[1])  # Start year
t2 = int(sys.argv[2])  # End year
variable = 'daily_maximum_2m_temperature'
reanalysis = 'era5'
jra = 'jra55'

# Directories
trend_dir = '/data/users/samuel.sancar/Heat_extremes/Trends/Reanalysis_trends/'
sig_dir = '/data/users/samuel.sancar/Heat_extremes/Trends/Trend_significance/'

# --- Load ERA5 trend and significance ---
era5_trend = iris.load_cube(f"{trend_dir}{reanalysis}_{variable}_tailwidth_trends_{t1}-{t2}.nc")
era5_sig = iris.load_cube(f"{trend_dir}{reanalysis}_{variable}_tailwidth_p_{t1}-{t2}.nc")

masked_trend_era5, _ = lat_land_mask(era5_trend, -60, 90)
masked_sig_era5, _ = lat_land_mask(era5_sig, -60, 90)
era5_data = masked_trend_era5.data[masked_sig_era5.data < 0.05].flatten()

# --- Load JRA-55 trend and significance ---
jra_trend = iris.load_cube(f"{trend_dir}{jra}__minmax_surf.015_tmax.reg_tl319._tailwidth_trends_{t1}-{t2}.nc")
jra_sig = iris.load_cube(f"{trend_dir}{jra}__minmax_surf.015_tmax.reg_tl319._tailwidth_p_{t1}-{t2}.nc")

masked_trend_jra, _ = lat_land_mask(jra_trend, -60, 90)
masked_sig_jra, _ = lat_land_mask(jra_sig, -60, 90)
jra_data = masked_trend_jra.data[masked_sig_jra.data < 0.05].flatten()

# --- Load model trends and significance ---
coupled_files = sorted(glob.glob('/data/users/samuel.sancar/Heat_extremes/Trends/Coupled_model_trends/*_tailwidth_trend.nc'))
sst_files = sorted(glob.glob('/data/users/samuel.sancar/Heat_extremes/Trends/SST_model_trends/*_tailwidth_trend.nc'))

coupled_sig_files = sorted(glob.glob(f"{sig_dir}P-vals/*_hist-1950_*_p_{t1}.nc"))
sst_sig_files = sorted(glob.glob(f"{sig_dir}P-vals/*_highresSST-present_*_p_{t1}.nc"))

bins = np.arange(-0.5, 1.5, 0.1)

# Lists for area and gridpoint calculations
significant_area = []
significant_gridpoints = []

plt.figure(figsize=(14, 6))

# --- Coupled models ---
for trend_file, sig_file in zip(coupled_files, coupled_sig_files):
    trend_cube = iris.load_cube(trend_file)
    sig_cube = iris.load_cube(sig_file)

    # Apply land mask
    masked_trend, _ = lat_land_mask(trend_cube, -60, 90)
    masked_sig, _ = lat_land_mask(sig_cube, -60, 90)

    # Filter significant trends
    model_data = masked_trend.data[masked_sig.data < 0.05].flatten()
    model_data = model_data[~np.isnan(model_data)]
    plt.hist(model_data, bins=bins, histtype='stepfilled', alpha=0.2, color='skyblue', edgecolor='black')

    # Ensure bounds for area weights
    if not sig_cube.coord('latitude').has_bounds():
        sig_cube.coord('latitude').guess_bounds()
    if not sig_cube.coord('longitude').has_bounds():
        sig_cube.coord('longitude').guess_bounds()
    if sig_cube.coord('latitude').units.is_unknown():
        sig_cube.coord('latitude').units = cf_units.Unit('degrees')
    if sig_cube.coord('longitude').units.is_unknown():
        sig_cube.coord('longitude').units = cf_units.Unit('degrees')

    area_weights = iris.analysis.cartography.area_weights(sig_cube)
    sig_mask = sig_cube.data < 0.05
    sig_cube_masked = sig_cube.copy()
    sig_cube_masked.data = np.ma.masked_where(~sig_mask, sig_cube.data)

    sig_area = sig_cube_masked.collapsed(['latitude', 'longitude'], iris.analysis.SUM, weights=area_weights).data
    significant_area.append(sig_area)
    significant_gridpoints.append(np.count_nonzero(sig_mask))

# --- SST models ---
for trend_file, sig_file in zip(sst_files, sst_sig_files):
    trend_cube = iris.load_cube(trend_file)
    sig_cube = iris.load_cube(sig_file)

    masked_trend, _ = lat_land_mask(trend_cube, -60, 90)
    masked_sig, _ = lat_land_mask(sig_cube, -60, 90)

    model_data = masked_trend.data[masked_sig.data < 0.05].flatten()
    model_data = model_data[~np.isnan(model_data)]
    plt.hist(model_data, bins=bins, histtype='stepfilled', alpha=0.2, color='lightgreen', edgecolor='black')

    if not sig_cube.coord('latitude').has_bounds():
        sig_cube.coord('latitude').guess_bounds()
    if not sig_cube.coord('longitude').has_bounds():
        sig_cube.coord('longitude').guess_bounds()
    if sig_cube.coord('latitude').units.is_unknown():
        sig_cube.coord('latitude').units = cf_units.Unit('degrees')
    if sig_cube.coord('longitude').units.is_unknown():
        sig_cube.coord('longitude').units = cf_units.Unit('degrees')

    area_weights = iris.analysis.cartography.area_weights(sig_cube)
    sig_mask = sig_cube.data < 0.05
    sig_cube_masked = sig_cube.copy()
    sig_cube_masked.data = np.ma.masked_where(~sig_mask, sig_cube.data)

    sig_area = sig_cube_masked.collapsed(['latitude', 'longitude'], iris.analysis.SUM, weights=area_weights).data
    significant_area.append(sig_area)
    significant_gridpoints.append(np.count_nonzero(sig_mask))

# --- ERA5 and JRA-55 histograms ---
plt.hist(era5_data, bins=bins, histtype='step', color='red', linewidth=2.5, linestyle='--', label='ERA5')
plt.hist(jra_data, bins=bins, histtype='step', color='blue', linewidth=2.5, linestyle='-.', label='JRA-55')

plt.xlabel('Trend Magnitude (°C per Decade)', fontsize=15, fontweight='bold')
plt.ylabel('Number of Grid Points', fontsize=15, fontweight='bold')
plt.title(f'Tail Width Trend Histogram (Significant Trends Only, {t1}–{t2})', fontsize=18, fontweight='bold')
plt.legend(loc='upper right', fontsize=13)
plt.grid(True, linestyle='--', alpha=0.6)
plt.tight_layout()
save_figure(plt, 'Significant_Tailwidth_Trend_Histogram', 'Graphs')
plt.show()

# --- Compute ERA5 and JRA-55 area and gridpoints ---
def compute_area_and_gridpoints(cube):
    if not cube.coord('latitude').has_bounds():
        cube.coord('latitude').guess_bounds()
    if not cube.coord('longitude').has_bounds():
        cube.coord('longitude').guess_bounds()
    if cube.coord('latitude').units.is_unknown():
        cube.coord('latitude').units = cf_units.Unit('degrees')
    if cube.coord('longitude').units.is_unknown():
        cube.coord('longitude').units = cf_units.Unit('degrees')
    area_weights = iris.analysis.cartography.area_weights(cube)
    sig_mask = cube.data < 0.05
    sig_cube = cube.copy()
    sig_cube.data = np.ma.masked_where(~sig_mask, cube.data)
    sig_area = sig_cube.collapsed(['latitude', 'longitude'], iris.analysis.SUM, weights=area_weights).data
    return sig_area, np.count_nonzero(sig_mask)

era5_area, era5_gridpoints = compute_area_and_gridpoints(era5_sig)
jra_area, jra_gridpoints = compute_area_and_gridpoints(jra_sig)

# --- Percentile calculations ---
significant_area_array = np.array(significant_area)
significant_gridpoints_array = np.array(significant_gridpoints)

era5_area_percentile = (np.sum(significant_area_array < era5_area) / len(significant_area_array)) * 100
jra_area_percentile = (np.sum(significant_area_array < jra_area) / len(significant_area_array)) * 100
era5_gridpoint_percentile = (np.sum(significant_gridpoints_array < era5_gridpoints) / len(significant_gridpoints_array)) * 100
jra_gridpoint_percentile = (np.sum(significant_gridpoints_array < jra_gridpoints) / len(significant_gridpoints_array)) * 100

print(f"ERA5 significant grid points: {era5_gridpoints}, Percentile: {era5_gridpoint_percentile:.2f}%")
print(f"JRA-55 significant grid points: {jra_gridpoints}, Percentile: {jra_gridpoint_percentile:.2f}%")
print(f"ERA5 significant area: {era5_area:.2e}, Percentile: {era5_area_percentile:.2f}%")
print(f"JRA-55 significant area: {jra_area:.2e}, Percentile: {jra_area_percentile:.2f}%")

# --- Second Panel: Histogram of significant area ---
plt.figure(figsize=(12, 7))
counts, bins, patches = plt.hist(
    significant_area_array,
    bins=10,
    range=(significant_area_array.min(), 3.2e11),
    histtype='bar',
    color='#4682B4',
    alpha=0.85,
    edgecolor='white',
    linewidth=2
)

plt.axvline(x=era5_area, color='#D7263D', linestyle='--', linewidth=2.5, label='ERA5')
plt.axvline(x=jra_area, color='#1B998B', linestyle='--', linewidth=2.5, label='JRA-55')

plt.xlabel('Area of Significant Trends (m²)', fontsize=15, fontweight='bold')
plt.ylabel('Number of Models', fontsize=15, fontweight='bold')
plt.title(f'Histogram of Total Area with Significant Tail-Width Trends ({t1}–{t2})', fontsize=18, fontweight='bold')
plt.legend(fontsize=13, loc='upper right', frameon=False)
plt.grid(True, linestyle='--', alpha=0.5, linewidth=1.2)
plt.tight_layout()
save_figure(plt, 'Significant_Area_Histogram', 'Graphs')
plt.show()