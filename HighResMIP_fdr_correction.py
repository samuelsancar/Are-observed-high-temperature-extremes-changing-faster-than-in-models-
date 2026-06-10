"""
22/09/25
---------

Apply Benjamini Hochberg metthod to HighResMIP members to correct for false discovery rate (alpha=0.1).

To run: python3 HighResMIP_FDR_correction.py
"""

import numpy as np
import os
import iris
import pandas as pd
from FDR_functions import fdr_from_masked_pvals, boolean_to_masked_cube
from Lat_Land_mask_function import lat_land_mask
from iris import Constraint

# Beginning & end year+1
t1 = 1980
t2 = 2023

alpha = 0.1

lat_minimum = -60
lat_maximum = 80

datadir = f'/data/users/samuel.sancar/Heat_extremes/Trends/Trend_significance/P-vals/'
outdir = "/data/users/samuel.sancar/Heat_extremes/Trends/Trend_significance/FDR_corrections/"

def _var_or_name(name):
    return Constraint(cube_func=lambda c: (getattr(c, 'var_name', None) == name) or (c.name() == name))

template = iris.load_cube('template_cube_1deg.nc')

# Read in excel sheet containing model names
df = pd.read_excel("/data/users/samuel.sancar/Model_table.xlsx", sheet_name='Model Stack', engine="openpyxl")

# ensure consistent column names (strip spaces)
df.columns = df.columns.str.strip()

processed = 0

for row in np.arange(0,63):
    
    row_idx = row  
    row = df.iloc[row_idx]

    # Map columns to variables 
    Mode       = str(row["Coupled/SST"]).strip()
    model      = str(row["Model"]).strip()
    experiment = str(row["Base experiment"]).strip()
    member     = str(row["Member"]).strip()

    print(model, experiment, member)
    
    # Load p-value cube
    model_name = f'{datadir}{model}_{experiment}_{member}_p_{t1}-{t2-1}.nc'
    Stat_sig = iris.load_cube(model_name, _var_or_name('Tailwidth_trend_p-values'))
    stat_sig_masked, no_of_gridpoints = lat_land_mask(Stat_sig, template, lat_minimum, lat_maximum)

    # Apply FDR
    fdr_sig = fdr_from_masked_pvals(stat_sig_masked, alpha)
    n_rejected = fdr_sig.sum()

    print(f"Number of significant grid points after FDR correction: "
        f"{n_rejected} out of {no_of_gridpoints} valid grid points")

    # Save all FDR-significant trends
    fdr_cube = boolean_to_masked_cube(fdr_sig, stat_sig_masked,name='fdr corrected trends',var_name='fdr_corrected_trends')

    iris.save(fdr_cube,os.path.join(outdir, f"{model}_{experiment}_{member}_fdr_corrected_p_{t1}-{t2-1}.nc"))

    # Positive-only trends
    tau_cube = iris.load_cube(f'{datadir}{model}_{experiment}_{member}_mktau_{t1}-{t2-1}.nc')
    tau_masked, _ = lat_land_mask(tau_cube, template, lat_minimum, lat_maximum)

    pos_sig = fdr_sig & (tau_masked.data > 0)

    pos_fdr_cube = boolean_to_masked_cube(pos_sig, stat_sig_masked,name='fdr corrected positive trends',var_name='fdr_corrected_positive_trends')

    iris.save(pos_fdr_cube,os.path.join(outdir, f"{model}_{experiment}_{member}_positive_fdr_corrected_p_{t1}-{t2-1}.nc"))