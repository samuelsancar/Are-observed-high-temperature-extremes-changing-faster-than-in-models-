"""
10/08/2025
------

 Function to calculate flase discovery rate corrected p-values using the Benjamin hochberg method.
 p-values are obtained using the Kendall-Tau non-parametric test.

 Returns: stipple_lats, stipple_lons, corrected_pvals, rejected

 """
 
import iris
import numpy as np
import scipy.stats as stats
import xarray as xr
from statsmodels.stats.multitest import fdrcorrection

def fdr_correction(k, alpha):

    # Apply FDR correction
    pvals_flat = k.values.flatten()
    rejected, corrected_pvals = fdrcorrection(pvals_flat, alpha=alpha)
    rejected_mask = rejected.reshape(k.shape)

    # Create significance mask
    significance_mask = xr.DataArray(data=rejected_mask, dims=k.dims, coords=k.coords, name="significance_mask")
    sig_mask_np = significance_mask.fillna(False).astype(bool).values

    # Extract stippling coordinates
    lon = k.coords['longitude'].values
    lat = k.coords['latitude'].values
    lon2d, lat2d = np.meshgrid(lon, lat)
    stipple_lats, stipple_lons = lat2d[sig_mask_np], lon2d[sig_mask_np]

    return stipple_lats, stipple_lons, corrected_pvals, rejected, sig_mask_np
