"""
10/08/2025
----------

Functions to calculate flase discovery rate corrected p-values using the Benjamin hochberg method.
p-values are obtained using the Kendall-Tau non-parametric test.

Returns: stipple_lats, stipple_lons, corrected_pvals, rejected
"""
 
import iris
import numpy as np
import scipy.stats as stats
import xarray as xr
from statsmodels.stats.multitest import fdrcorrection

def boolean_to_masked_cube(bool_field, template_cube, name, var_name):
    cube = template_cube.copy()
    cube.data = np.ma.masked_array(
        bool_field.astype(np.int8),
        mask=template_cube.data.mask
    )
    cube.long_name = name
    cube.var_name = var_name
    return cube

def fdr_from_masked_pvals(pval_cube, alpha):
    """
    Apply BH-FDR to a masked Iris cube of p-values.
    Returns a masked boolean array (True = significant).
    """
    pvals = pval_cube.data  # masked array

    valid = (~pvals.mask) & np.isfinite(pvals.data)
    pvals_flat = pvals.data[valid]

    rejected, _ = fdrcorrection(pvals_flat, alpha=alpha)

    sig = np.zeros(pvals.shape, dtype=bool)
    sig[valid] = rejected

    return sig


def stipple_from_cube(cube, lat2d, lon2d):
    data = cube.data
    mask = (~data.mask) & data.data.astype(bool)
    return lat2d[mask], lon2d[mask]

