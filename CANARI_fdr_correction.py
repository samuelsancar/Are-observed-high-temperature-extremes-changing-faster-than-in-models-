"""
10/11/25
---------

Bejnamini Hochberg method applied to p-values CANARI members to correct for false discovery rate (alpha=0.1).

To run: python CANARI_fdr_correction.py
"""

import numpy as np
import os
import iris
from FDR_functions import fdr_from_masked_pvals, boolean_to_masked_cube
from Lat_Land_mask_function import lat_land_mask
from iris import Constraint

# Beginning & end year+1
t1 = 1958
t2 = 2023

alpha = 0.1

lat_minimum = -60
lat_maximum = 80

datadir = '/data/users/samuel.sancar/Heat_extremes/Trends/Trend_significance/CANARI/P-vals/'
outdir = "/data/users/samuel.sancar/Heat_extremes/Trends/Trend_significance/CANARI/FDR_corrections/"

def _var_or_name(name):
    return Constraint(cube_func=lambda c: (getattr(c, 'var_name', None) == name) or (c.name() == name))

template = iris.load_cube('template_cube_1deg.nc')

for i in np.arange(1,41): # 40 members
    
    print(f"Processing file: {i}")
    
     # Load p-value cube
    model_name = f'{datadir}{i}_p_{t1}-{t2-1}.nc'
    Stat_sig = iris.load_cube(model_name, _var_or_name('Tailwidth_trend_p-values'))
    stat_sig_masked, no_of_gridpoints = lat_land_mask(Stat_sig, template, lat_minimum, lat_maximum)

    # Apply FDR
    fdr_sig = fdr_from_masked_pvals(stat_sig_masked, alpha)
    n_rejected = fdr_sig.sum()

    print(f"Number of significant grid points after FDR correction: "
          f"{n_rejected} out of {no_of_gridpoints} valid grid points")

    # Save all FDR-significant trends
    fdr_cube = boolean_to_masked_cube(fdr_sig, stat_sig_masked,name='fdr corrected trends',var_name='fdr_corrected_trends')

    iris.save(fdr_cube,os.path.join(outdir, f"{i}_fdr_corrected_p_{t1}-{t2-1}.nc"))

    # Positive-only trends
    tau_cube = iris.load_cube(f'{datadir}{i}_mktau_{t1}-{t2-1}.nc')
    tau_masked, _ = lat_land_mask(tau_cube, template, lat_minimum, lat_maximum)

    pos_sig = fdr_sig & (tau_masked.data > 0)

    pos_fdr_cube = boolean_to_masked_cube(pos_sig, stat_sig_masked,name='fdr corrected positive trends',var_name='fdr_corrected_positive_trends')

    iris.save(pos_fdr_cube,os.path.join(outdir, f"{i}_positive_fdr_corrected_p_{t1}-{t2-1}.nc"))