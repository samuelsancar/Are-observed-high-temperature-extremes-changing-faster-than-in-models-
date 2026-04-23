"""
Code to calculate FDR correction for trends in tailwidths for JRA-55 reanalysis data.
"""

import iris
import matplotlib.cm as mpl_cm
import numpy as np
import os
import sys
from FDR_functions import fdr_from_masked_pvals, boolean_to_masked_cube
from Lat_Land_mask_function import lat_land_mask
from iris import Constraint

cmap_rb = mpl_cm.get_cmap('RdBu_r', 15)

# Beginning & end year+1
t1 = 1980
t2 = 2023

alpha = 0.1

lat_minimum = -60
lat_maximum = 80

reanalysis = 'jra55'
variable = '_minmax_surf.015_tmax.reg_tl319.'
dir = '/data/users/samuel.sancar/Heat_extremes/Trends/Reanalysis_trends/'

Stat_sig_name = f"{dir}{reanalysis}_{variable}_tailwidth_p_{t1}-{t2-1}.nc"
Mk_tau_name = f"{dir}{reanalysis}_{variable}_tailwidth_mktau_{t1}-{t2-1}.nc"

template = iris.load_cube('template_cube_1deg.nc')

def _var_or_name(name):
    return Constraint(cube_func=lambda c: (getattr(c, 'var_name', None) == name) or (c.name() == name))

# Load and mask p-values
Stat_sig = iris.load_cube(Stat_sig_name, _var_or_name('Tailwidth_trend_p-values'))
stat_sig_masked, no_of_gridpoints = lat_land_mask(Stat_sig, template, lat_minimum, lat_maximum)

# Apply FDR
fdr_sig = fdr_from_masked_pvals(stat_sig_masked, alpha)
n_rejected = fdr_sig.sum()

print(f"Number of significant grid points after FDR correction: "
      f"{n_rejected} out of {no_of_gridpoints} valid grid points")

# Save all FDR-significant trends
fdr_cube = boolean_to_masked_cube(fdr_sig, stat_sig_masked,name='fdr corrected trends',var_name='fdr_corrected_trends')

iris.save(fdr_cube,os.path.join(dir, f"{reanalysis}_{variable}_fdr_corrected_p_{t1}-{t2-1}.nc"))

# Positive-only trends
tau_cube = iris.load_cube(Mk_tau_name)
tau_masked, _ = lat_land_mask(tau_cube, template, lat_minimum, lat_maximum)

pos_sig = fdr_sig & (tau_masked.data > 0)

pos_fdr_cube = boolean_to_masked_cube(pos_sig, stat_sig_masked,name='fdr corrected positive trends',var_name='fdr_corrected_positive_trends')

iris.save(pos_fdr_cube,os.path.join(dir, f"{reanalysis}_{variable}_positive_fdr_corrected_p_{t1}-{t2-1}.nc"))