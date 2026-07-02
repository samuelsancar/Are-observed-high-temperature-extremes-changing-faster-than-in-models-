## (C) Crown Copyright, Met Office. All rights reserved.
##
## This file is part of 'Are observed high temperature extremes changing faster than in models?' and is released under the BSD 3-Clause license.
## See LICENSE in the root of the repository for full licensing details.

# Are high temperature extremes changing faster than in models?

This repository contains the code required to reproduce all main and supplementary figures from Sancar et al 2026.

## Overview

Given the required data, this code calculates and analyses the tail widths trends of daily maximum temperatures (over the period 1958-2022) from: ERA5, JRA-55, MERRA2, BerkeleyEarth (not included in paper), 63 HighResMIP members (from multiple models), and 40 members from a single model initial condition large ensemble (SMILE) from the CANARI project (HadGEM3-GC31-MM model).

## Data and directories

### Reanalyses: ERA5, JRA-55, MERRA.
- Code is setup to calculate tail width trends from the raw imported files of daily maximum temperatures from each reanalysis. Each reanalysis must be stored in individual directories.
- Directories are set up so that data must be saved to the following convention: /data/users/samuel.sancar/Climate_models/[reanalysis].
  
### Observational Dataset: BerkeleyEarth (This dataset is not formally documented or peer-reviewed and is therefore not included in the main analysis).
- Directories are set up so that data must be saved to the following convention: /data/users/samuel.sancar/Climate_models/BerkeleyEarth.
  
### HighResMIP: Details are found in the Model_table.xlsx.
- Code is setup to calculate tail width trends from the raw imported files of daily maximum temperatures from each HighResMIP member.
- Excel spread sheet of model members and their details is supplied. Code utilises this spreadsheet to read in correct members from the correct directories. The spreadsheet can be modified to reflect the subset of model members available to the user.
- Directories are set up so that data must be saved to the following convention: /data/users/samuel.sancar/Climate_models/[Coupled|SST]/[Institution]/[member]. Therfore, each institution/model has it's own directory, with each member from each model having it's own directory within that (containing all the data). This can be changed
  
### SMILE: HadGEM3-GC31-MM, 40 members, from the CANARI project.
- Code is setup to calculate tail width trends from the raw imported files of daily maximum temperatures from each CANARI large ensemble member.
-  Directories are set up so that data must be saved to the following convention: /data/users/samuel.sancar/Climate_models/CANARI/[member]/.

Due to data volume and access restrictions, datasets are not included in this repository.  
Users should refer to the Data Availability Statement in Sancar et al. (2026) for where to access data.

Directory paths can be modified within the scripts to suit local configurations.

## Main workflow

### ERA5 analysis

ERA5_tailwidths.py -> ERA5_trend_calc.py -> ERA5_MKtest.py -> ERA5_fdr_correction.py -> ERA5_trend_plotter.py (figure 1a).

### JRA-55 analysis

JRA-55_tailwidths.py -> JRA-55_trend_calc.py -> JRA-55_MKtest.py -> JRA-55_fdr_correction.py -> JRA-55_trend_plotter.py (figure 1b) -> ERA5_JRA_significant_trends.py (figure 5a & 5b).

### MERRA2 analysis

MERRA_tailwidths.py -> MERRA_trend_calc.py -> MERRA_MKtest.py -> MERRA_fdr_correction.py -> MERRA_trend_plotter.py -> ERA5_MERRA_significant_trends.py (figure S3a & S3b)

### HighResMIP analysis

HighResMIP_tailwidths.py -> HighResMIP_stack.py -> HighResMIP_calc.py -> HighResMIP_MKtest.py -> HighResMIP_fdr_correction.py -> HighResMIP_trend_plotter.py (figure 1c, 1d, S2a & S2b)
-> HighResMIP_trends_histogram.py (figure 2a) -> HighResMIP_cumulative_trends.py (figure 2b) -> HighResMIP_mmm_tailwidths.py -> HighResMIP_mmm_trends_calc.py -> HighResMIP_std 
-> HighResMIP_mmm_trend_plotter (figure S1a) -> HighResMIP_mmm_trend_plotter_2.0.py (figure S1b).

### CANARI large ensemble analysis

CANARI_tailwidths.py -> CANARI_trend_calc.py -> CANARI_MKtest.py -> CANARI_trend_plotter.py -> CANARI_agreement.py -> CANARI_intercomparison.py (figure 3a-d) -> CANARI_region_selection.py (figure 3e & 3f)
CANARI_region_analysis.py (figure 4a-i).

## Extra analysis
- Code for extra analysis is available in the 'Extra_analysis' respository. Here we analyse and plot trends in the tail width of the BerkeleyEarth dataset.



