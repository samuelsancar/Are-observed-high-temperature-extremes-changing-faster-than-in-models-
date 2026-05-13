""" 
Threshold curve: to go at end of CANARI_intercomparison.py
"""

np.save(outdir + f'{model_name}_agreement.npy', agreement_land)

    # Calculate threshold curve
    thresholds = np.arange(0, 1.01, 0.01)
    fraction_gridpoints = [
        100 * np.count_nonzero((agreement_land >= thresh) & ~np.isnan(agreement_land)) / n_gridpoints
        for thresh in thresholds
    ]
    
    OHpcnt_gridpoints = 100 * np.count_nonzero((agreement_land >= 1) & ~np.isnan(agreement_land)) / n_gridpoints
    
    np.save(outdir + f'{model_name}_agreement.npy', agreement_land)
    print(f"Model {model_name} has {OHpcnt_gridpoints}% of the globe covered by grid points where it's trend is greater than 100% of models.")
    
    all_fraction_gridpoints.append((fname, fraction_gridpoints))
    OHpcnt_list.append((fname, OHpcnt_gridpoints))
    
    if OHpcnt_gridpoints >= ERA5_OHpcnt:
        labels_list.append(fname)
    
    full = 100 * np.count_nonzero((agreement_land >= 1) & ~np.isnan(agreement_land)) / n_gridpoints
    number = np.count_nonzero((agreement_land >= 1) & ~np.isnan(agreement_land))
    print(f"Model {model_name} has {full}% of the globe ({number} grid points) covered by grid points where it's trend is greater than 100% of models.")
    
    fraction_area = []
    
    for thresh in thresholds:
        
        grid_points = (agreement_land >= thresh) & (~np.isnan(agreement_land))
        cube_masked = agreement_land_cube.copy()
        cube_masked.data = grid_points.astype(int)  

        if not cube_masked.coord('latitude').has_bounds():
            cube_masked.coord('latitude').guess_bounds()
        if not cube_masked.coord('longitude').has_bounds():
            cube_masked.coord('longitude').guess_bounds()
        
        # Fix units if missing or unknown
        if cube_masked.coord('latitude').units.is_unknown():
            cube_masked.coord('latitude').units = cf_units.Unit('degrees')
        if cube_masked.coord('longitude').units.is_unknown():
            cube_masked.coord('longitude').units = cf_units.Unit('degrees')

        area_weights = iris.analysis.cartography.area_weights(cube_masked)

        area = cube_masked.collapsed(['latitude', 'longitude'], iris.analysis.SUM, weights=area_weights).data
        
        # Convert to percentage of total area
        total_area = np.sum(area_weights)
        fraction_area.append(100 * area / total_area)
        
    fraction_area_list.append((fname, fraction_area))

# Plot all threshold curves together
plt.figure(figsize=(8, 8))
y = -np.arange(0, 101, 1) + 100

for fname, fraction_gridpoints in all_fraction_gridpoints:
    plt.plot(np.arange(0, 101, 1), fraction_gridpoints, marker='.', color = 'grey') 

plt.plot(np.arange(0, 101, 1), y, linestyle='--', color='black', label='1:1 Line')
plt.plot(np.arange(0, 101, 1), ERA5, marker='.', color='red', label='ERA5 agreement')
#plt.plot(np.arange(0, 101, 1), JRA55, marker='.', color='blue', label='JRA-55 agreement')

plt.xlabel('Percentage of models where trends are less X (where X is an ensemble member, ERA5 or JRA-55) per grid point')
plt.ylabel('Percentage of Land Grid Points (%)')
plt.title('Model Agreement Curves for ERA5 and SST Models')
plt.grid(True)
plt.legend(loc='best', fontsize='small')
plt.tight_layout()
plt.show()

# Plot all threshold curves together
plt.figure(figsize=(8, 8))
y = -np.arange(0, 101, 1) + 100

for fname, fraction_area in fraction_area_list:
    plt.plot(np.arange(0, 101, 1), fraction_area, marker='.', color = 'grey') 

plt.plot(np.arange(0, 101, 1), y, linestyle='--', color='black', label='1:1 Line')
plt.plot(np.arange(0, 101, 1), ERA5_area, marker='.', color='red', label='ERA5 agreement')
#plt.plot(np.arange(0, 101, 1), JRA55, marker='.', color='blue', label='JRA-55 agreement')

plt.xlabel('Percentage of models where trends are less X (where X is an ensemble member, ERA5 or JRA-55) per grid point')
plt.ylabel('Percentage of Land Grid Points (%)')
plt.title('Model Agreement Curves for ERA5 and SST Models')
plt.grid(True)
plt.legend(loc='best', fontsize='small')
plt.tight_layout()
plt.show()

OHpcnt_values = [val for _, val in OHpcnt_list]  

OHpcnt_array = np.array(OHpcnt_values, dtype=float)
ERA5_OHpcnt = float(ERA5_OHpcnt)  

# Calculate ERA5 percentile
era5_percentile = (np.sum(OHpcnt_array < ERA5_OHpcnt)) / len(OHpcnt_array) * 100
print(f"ERA5 is in the {era5_percentile:.2f} percentile of model members for fraction of grid points where it has the greatest trend")

plt.style.use('ggplot')

fig, ax = plt.subplots(figsize=(7, 7))

# Boxplot with custom colours
bp = ax.boxplot(
    OHpcnt_array,
    patch_artist=True,
    boxprops=dict(facecolor='#4682B4', alpha=0.85, edgecolor='white', linewidth=2),
    whiskerprops=dict(color='black', linewidth=2),
    capprops=dict(color='black', linewidth=2),
    medianprops=dict(color='black', linewidth=2)
)
ax.scatter([1] * len(OHpcnt_array), OHpcnt_array, color = 'blue', alpha = 0.5, marker = 'o')

# Scatter point for ERA5
ax.scatter(1, ERA5_OHpcnt, color='#D7263D', edgecolor='#D7263D', alpha=1, s=100, label='ERA5', marker='x')

# Labels and title
ax.set_xlabel('', fontsize=15, fontweight='bold')
ax.set_ylabel('Percentage of Grid Points with Greatest Trend (%)', fontsize=15, fontweight='bold')
ax.set_title('Boxplot of % of Grid Points with Greatest Trend', fontsize=18, fontweight='bold')

# Ticks and grid
ax.tick_params(axis='both', labelsize=13)
ax.grid(True, linestyle='--', alpha=0.5, linewidth=1.2)

# Legend
ax.legend(fontsize=13, loc='upper right', frameon=True)

plt.tight_layout()
plt.show()

