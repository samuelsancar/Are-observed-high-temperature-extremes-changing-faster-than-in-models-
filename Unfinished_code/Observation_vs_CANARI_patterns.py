import iris
import numpy as np
from scipy.stats import pearsonr
import matplotlib.pyplot as plt
import seaborn as sns
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from Lat_Land_mask_function import lat_land_mask

# --- CONFIG ---
REANALYSIS = 'era5'
VARIABLE = 'daily_maximum_2m_temperature'
T1, T2 = 1958, 2023
LAT_MIN, LAT_MAX = -60, 90
DATADIR = '/data/users/samuel.sancar/Heat_extremes/Trends/Reanalysis_trends/'

# --- LOAD TEMPLATE AND MODEL DATA ---
template = iris.load_cube('template_cube_1deg.nc')
canari_files = iris.load('/data/users/samuel.sancar/Heat_extremes/Trends/CANARI_trends/*_tailwidth_trend.nc')
canari_regridded = [cube.regrid(template, iris.analysis.Linear()) for cube in canari_files]

# --- LOAD ERA5 TREND ---
era5_trend = iris.load_cube(f"{DATADIR}{REANALYSIS}_{VARIABLE}_tailwidth_trends_{T1}-{T2}.nc")
era5_trend = era5_trend.regrid(template, iris.analysis.Linear())

# --- APPLY LATITUDE AND LAND MASK ---
model_mask, _ = lat_land_mask(era5_trend, LAT_MIN, LAT_MAX)
model_mask_array = np.array(model_mask, dtype=bool)

# --- CLEAN DATA FUNCTION ---
def clean_data(cube):
    """Convert masked arrays to NaN and remove extreme outliers."""
    data = cube.data
    if hasattr(data, 'mask'):
        data = np.ma.filled(data, np.nan)
    data = np.where(np.abs(data) > 1e10, np.nan, data)
    return data

# --- PREP ERA5 DATA ---
era5_data = clean_data(era5_trend)
era5_flat = era5_data[model_mask_array].flatten()

# --- CORRELATIONS FOR EACH MEMBER ---
correlations = []
for cube in canari_regridded:
    member_data = clean_data(cube)
    member_flat = member_data[model_mask_array].flatten()
    valid = ~np.isnan(era5_flat) & ~np.isnan(member_flat)
    if np.sum(valid) > 0:
        corr, _ = pearsonr(era5_flat[valid], member_flat[valid])
        correlations.append(corr)
    else:
        correlations.append(np.nan)

# --- ENSEMBLE MEAN CORRELATION ---
ensemble_mean_full = np.nanmean([clean_data(cube) for cube in canari_regridded], axis=0)
ensemble_flat = ensemble_mean_full[model_mask_array].flatten()
valid = ~np.isnan(era5_flat) & ~np.isnan(ensemble_flat)
corr_mean, _ = pearsonr(era5_flat[valid], ensemble_flat[valid])

# --- HISTOGRAM OF CORRELATIONS ---
sns.histplot(correlations, bins=20, kde=True)
plt.axvline(corr_mean, color='red', linestyle='--', linewidth=2,
            label=f'Ensemble Mean Correlation = {corr_mean:.2f}')
plt.title('Pattern Correlation: ERA5 vs CANARI Members')
plt.xlabel('Correlation')
plt.ylabel('Frequency')
plt.legend()
plt.show()

import iris
import numpy as np
import pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt

# --- CONFIG ---
REANALYSIS = 'era5'
VARIABLE = 'daily_maximum_2m_temperature'
T1, T2 = 1958, 2023
LAT_MIN, LAT_MAX = -60, 90
DATADIR = '/data/users/samuel.sancar/Heat_extremes/Trends/Reanalysis_trends/'

# --- LOAD TEMPLATE AND MODEL DATA ---
template = iris.load_cube('template_cube_1deg.nc')
canari_files = iris.load('/data/users/samuel.sancar/Heat_extremes/Trends/CANARI_trends/*_tailwidth_trend.nc')
canari_regridded = [cube.regrid(template, iris.analysis.Linear()) for cube in canari_files]

# --- LOAD ERA5 TREND ---
era5_trend = iris.load_cube(f"{DATADIR}{REANALYSIS}_{VARIABLE}_tailwidth_trends_{T1}-{T2}.nc")
era5_trend = era5_trend.regrid(template, iris.analysis.Linear())

# --- APPLY LATITUDE AND LAND MASK ---
model_mask, _ = lat_land_mask(era5_trend, LAT_MIN, LAT_MAX)
model_mask_array = np.array(model_mask, dtype=bool)

# Apply mask to ERA5
era5_masked = era5_trend.data[model_mask_array]

# Apply mask to CANARI members
canari_masked = [cube.data[model_mask_array] for cube in canari_regridded]

# --- FLATTEN INTO MATRIX ---
n_members = len(canari_masked)
grid_size = era5_masked.size
data_matrix = np.zeros((n_members + 1, grid_size))  # +1 for ERA5

for i, member_data in enumerate(canari_masked):
    data_matrix[i, :] = member_data.flatten()
data_matrix[-1, :] = era5_masked.flatten()

# --- STANDARDISE ---
# --- HANDLE NaNs ---
valid_cols = ~np.isnan(data_matrix).any(axis=0)
data_matrix_clean = data_matrix[:, valid_cols]

# --- STANDARDISE ---
scaler = StandardScaler()
data_scaled = scaler.fit_transform(data_matrix_clean)

# --- CLUSTER USING KMeans ---
kmeans = KMeans(n_clusters=5, random_state=42)
clusters = kmeans.fit_predict(data_scaled)

# --- ADD LABELS ---
labels = [cube.attributes.get('model_id', f'Member_{i}') for i, cube in enumerate(canari_regridded)] + ['ERA5']
cluster_df = pd.DataFrame({'Model': labels, 'Cluster': clusters})
print(cluster_df)

# --- PCA FOR VISUALISATION ---
pca = PCA(n_components=2)
pca_coords = pca.fit_transform(data_scaled)

plt.figure(figsize=(8, 6))
for cluster in np.unique(clusters):
    idx = clusters == cluster
    plt.scatter(pca_coords[idx, 0], pca_coords[idx, 1], label=f'Cluster {cluster}')

# Highlight ERA5
plt.scatter(pca_coords[-1, 0], pca_coords[-1, 1], color='black', marker='*', s=200, label='ERA5')

# Add axis labels and title
plt.xlabel(f'PC1 ({pca.explained_variance_ratio_[0]*100:.1f}% variance)')
plt.ylabel(f'PC2 ({pca.explained_variance_ratio_[1]*100:.1f}% variance)')
plt.title('Model Clusters (PCA projection)')
plt.legend()
plt.show()