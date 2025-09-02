"""
17/07/25
---------
Code to demonstrate how FDR and raw p-values changes with alpha level.

To run: python FDRvsP-Value.py

"""

import iris
import numpy as np
import matplotlib.pyplot as plt
from statsmodels.stats.multitest import fdrcorrection
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from Lat_Land_mask_function import lat_land_mask

# Define alpha levels
alpha_levels = [0.2, 0.15, 0.10, 0.05, 0.04, 0.03, 0.02, 0.01, 0.005, 0.0025, 0.001, 0.00025, 0]

# Load p-value data
datadir = '/data/users/samuel.sancar/Heat_extremes/trends/Reanalysis_trends/'
reanalysis = 'era5'
variable = 'daily_maximum_2m_temperature'
t1, t2 = 1958, 2023
pval_cube = iris.load_cube(f"{datadir}{reanalysis}_{variable}_tailwidth_p_{t1}-{t2}.nc")

# Apply lat_land_mask to p-values
masked_pval_cube, num_land_points = lat_land_mask(pval_cube, -60)
masked_pvals = masked_pval_cube.data
pvals_flat = masked_pvals.flatten()
pvals_flat = pvals_flat[~np.isnan(pvals_flat)]

# Initialise lists to store FDR and raw p-value percentages
fdr_percentages = []
raw_pval_percentages = []

# Compute FDR and raw p-value coverage
for alpha in alpha_levels:
    # FDR correction
    rejected, _ = fdrcorrection(pvals_flat, alpha=alpha)
    fdr_percentage = 100 * np.sum(rejected) / num_land_points
    fdr_percentages.append(fdr_percentage)

    # Raw p-value thresholding
    raw_rejected = pvals_flat < alpha
    raw_percentage = 100 * np.sum(raw_rejected) / num_land_points
    raw_pval_percentages.append(raw_percentage)

# Plot FDR vs raw p-value coverage
plt.figure(figsize=(12, 7))
plt.plot(alpha_levels, fdr_percentages, marker='o', linestyle='-', linewidth=2, markersize=8,
         label='FDR Corrected', color='royalblue')
plt.plot(alpha_levels, raw_pval_percentages, marker='s', linestyle='--', linewidth=2, markersize=8,
         label='Raw p-value', color='darkorange')
plt.xlabel('Alpha Level', fontsize=14)
plt.ylabel('Percentage of Globe Covered (%)', fontsize=14)
plt.title('Trend Significance Coverage by Alpha Level', fontsize=16, weight='bold')
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.legend(fontsize=12)
plt.grid(True, linestyle='--', alpha=0.6)
plt.tight_layout()
plt.savefig("fdr_vs_raw_pvalue_plot.png", dpi=300)
plt.show()

# Generate stippling map for alpha = 0.10
alphatwo = 0.10
rejected, _ = fdrcorrection(pvals_flat, alpha=alphatwo)
raw_rejected = pvals_flat < 0.05

# Create boolean masks
fdr_mask_flat = rejected
raw_mask_flat = raw_rejected
difference_mask_flat = raw_mask_flat & (~fdr_mask_flat)

# Reshape to original grid shape
latitudes = masked_pval_cube.coord('latitude').points
longitudes = masked_pval_cube.coord('longitude').points
lon2d, lat2d = np.meshgrid(longitudes, latitudes)
grid_shape = masked_pval_cube.shape
difference_mask = np.full(grid_shape, False)
difference_mask_flat_full = np.full(grid_shape[0] * grid_shape[1], False)
valid_indices = ~np.isnan(masked_pvals.flatten())
difference_mask_flat_full[valid_indices] = difference_mask_flat
difference_mask = difference_mask_flat_full.reshape(grid_shape)

# Extract stippling coordinates
stipple_lats = lat2d[difference_mask]
stipple_lons = lon2d[difference_mask]

# Plot stippling map
fig = plt.figure(figsize=(14, 7))
ax = plt.axes(projection=ccrs.PlateCarree())
ax.set_global()
ax.scatter(stipple_lons, stipple_lats, color='crimson', s=1, alpha=0.6,
           transform=ccrs.PlateCarree(), label='Stippling')
ax.coastlines(resolution='50m', linewidth=1)
ax.add_feature(cfeature.LAND, facecolor='lightgray')
ax.add_feature(cfeature.OCEAN, facecolor='aliceblue')
ax.add_feature(cfeature.BORDERS, linestyle=':', linewidth=0.5)
gl = ax.gridlines(draw_labels=True, linewidth=0.5, color='gray', alpha=0.5, linestyle='--')
gl.top_labels = False
gl.right_labels = False
gl.xlabel_style = {'size': 10}
gl.ylabel_style = {'size': 10}
plt.title(f"Rejected Statistically Significant Grid Points\n{reanalysis.upper()} {variable.replace('_', ' ')} ({t1}–{t2})",
          fontsize=16, weight='bold')
plt.tight_layout()
plt.savefig("stippling_map.png", dpi=300)
plt.show()






