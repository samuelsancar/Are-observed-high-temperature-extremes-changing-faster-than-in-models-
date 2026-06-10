"""
13/11/2025
---------- 

ERA5 tail width trends vs CANARI model trends.
At each grid point, calculate the percentage of CANARI models with a tail width trend less than or equal to the ERA5 trend. 

To run: python CANARI_agreement.py
"""

import iris
import iris.analysis
import numpy as np
import matplotlib.pyplot as plt
import os
from Unused_code.BnD_error_function import binomial_uncertainty
from Lat_Land_mask_function import lat_land_mask

# Beginning and end year+1
t1 = 1958
t2 = 2023

lat_minimum = -60
lat_maximum = 90

# Load template
template = iris.load_cube("template_cube_1deg.nc")

# Load and mask ERA5
era5 = iris.load_cube(
    f"/data/users/samuel.sancar/Heat_extremes/Trends/Reanalysis_trends/era5_daily_maximum_2m_temperature_tailwidth_trends_{t1}-{t2-1}.nc")

era5_cube, _ = lat_land_mask(era5, template, lat_minimum, lat_maximum)
era5_data = era5_cube.data

# Load and mask CANARI members 
CANARI_files = iris.load(
    f"/data/users/samuel.sancar/Heat_extremes/Trends/CANARI_trends/*_tailwidth_trend_{t1}-{t2-1}.nc"
)

model_stack = np.array([lat_land_mask(cube, template, lat_minimum, lat_maximum)[0].data for cube in CANARI_files])

n_models = model_stack.shape[0]

# Agreement field
agreement = np.mean(model_stack <= era5_data, axis=0)

agreement_cube = era5_cube.copy()
agreement_cube.data = agreement

agreement_data = agreement_cube.data

# Valid grid points
valid = ~np.ma.getmaskarray(agreement_data)
n_gridpoints = np.count_nonzero(valid)

# Sampling uncertainty (binomial, ensemble size)
p_low_list = []
p_high_list = []

for i in range(101):
    k = (n_models / 100) * i if i > 0 else 0
    p_low, p_high, *_ = binomial_uncertainty(n_models, k)
    p_low_list.append(100 * p_low / n_models)
    p_high_list.append(100 * p_high / n_models)

# Threshold scan
thresholds = np.arange(0, 1.01, 0.01)
fraction_gridpoints = []
fraction_area_list = []

# Area weights (constant across thresholds)
if not agreement_cube.coord("latitude").has_bounds():
    agreement_cube.coord("latitude").guess_bounds()
if not agreement_cube.coord("longitude").has_bounds():
    agreement_cube.coord("longitude").guess_bounds()

area_weights = iris.analysis.cartography.area_weights(agreement_cube)
total_area = np.sum(area_weights[valid])

for thresh in thresholds:
    condition = (agreement_data >= thresh) & valid

    # Grid-point fraction
    fraction_gridpoints.append(100 * np.count_nonzero(condition) / n_gridpoints)

    # Area fraction
    area = np.sum(area_weights[condition])
    fraction_area_list.append(100 * area / total_area)

# ERA5 extreme coverage (agreement == 1)
extreme_condition = (agreement_data == 1) & valid
OHpcnt = 100 * np.count_nonzero(extreme_condition) / n_gridpoints
number = np.count_nonzero(extreme_condition)

print(
    f"ERA5 exceeds all models over {OHpcnt:.2f}% of land "
    f"({number} grid points)"
)

# Save outputs
outpath = "/data/users/samuel.sancar/Heat_extremes/Trends/Trend_agreement/"
os.makedirs(outpath, exist_ok=True)

np.save(os.path.join(outpath, "CANARI_ERA5_fraction_gridpoints.npy"), fraction_gridpoints)
np.save(os.path.join(outpath, "CANARI_ERA5_fraction_area.npy"), fraction_area_list)
np.save(os.path.join(outpath, "CANARI_ERA5_agreement_field.npy"),agreement_data.filled(np.nan))
np.save(os.path.join(outpath, "ERA5_extreme_fraction.npy"), OHpcnt)

# Plots
y = 100 - np.arange(101)

# Area agreement curve
plt.figure(figsize=(12, 7))
plt.plot(thresholds * 100, fraction_area_list, color="red", label="Model agreement (area)")
plt.plot(thresholds * 100, y, "--", color="gray", label="1:1 line")
plt.xlabel("Percentage of models ≤ ERA5 (%)", fontsize=14, weight="bold")
plt.ylabel("Percentage of land area (%)", fontsize=14, weight="bold")
plt.title("Model agreement with ERA5 (area-weighted)", fontsize=16, weight="bold")
plt.legend()
plt.grid(True, linestyle="--", alpha=0.6)
plt.tight_layout()
plt.show()

# Grid-point agreement curve with sampling uncertainty
plt.figure(figsize=(12, 7))
plt.errorbar(
    thresholds * 100,
    fraction_gridpoints,
    xerr=[p_low_list, p_high_list],
    fmt=".",
    capsize=4
)
plt.plot(thresholds * 100, y, "--", color="gray")
plt.xlabel("Percentage of models ≤ ERA5 (%)", fontsize=14, weight="bold")
plt.ylabel("Percentage of land grid points (%)", fontsize=14, weight="bold")
plt.title("Model agreement with ERA5 (grid points)", fontsize=16, weight="bold")
plt.grid(True, linestyle="--", alpha=0.6)
plt.tight_layout()
plt.show()




