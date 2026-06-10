"""
15/05/2026
----------

Code to plot scatter plots for HighResMIP.
Scatter plot for each HighResMIP member comparing its tail width trend to the HighResMIP multi-model mean tail width trend.

To run: python HighResMIP_intercomparison_scatter_plot.py
"""

import iris
import numpy as np
import sys
from iris.cube import CubeList
import iris.plot as iplt
import cartopy.crs as ccrs
import matplotlib.cm as mpl_cm
import sys
import pandas as pd
import matplotlib.pyplot as plt
from Lat_Land_mask_function import lat_land_mask

from Save_figure_function import save_figure

cmap_rb = mpl_cm.get_cmap('coolwarm', 15)

# Beginning and end years
t1 = 1958
t2 = 2023

lat_minimum = -60
lat_maximum = 60

years = np.arange(t1, t2)

        
        
template = iris.load_cube('template_cube_1deg.nc')
df = pd.read_excel("/data/users/samuel.sancar/Model_table.xlsx", sheet_name=0, engine="openpyxl")

# ensure consistent column names (strip spaces)
df.columns = df.columns.str.strip()

all_slopes = []
all_corr = []

for i in np.arange(0,63):
    
    all_trend = []
    
    for row in np.arange(0, 63): 
        
        row_idx = row  
        print(row)
        row = df.iloc[row_idx]

        # Map columns to variables 
        Mode       = str(row["Coupled/SST"]).strip()
        Group      = str(row["Institution"]).strip()
        model      = str(row["Model"]).strip()
        experiment = str(row["Base experiment"]).strip()
        member     = str(row["Member"]).strip()

        print(Mode, Group, model, experiment, member)

        datadir = f'/data/users/samuel.sancar/Heat_extremes/Trends/{Mode}_model_trends/'
        filename = f"{model}_{experiment}_{member}_tailwidth_trend_1958-2022.nc"
        
        if row_idx == i:
            chosen = iris.load_cube(datadir+filename)
            chosen_p87p5 = iris.load_cube(datadir + filename.replace('tailwidth','p87p5'))
            i_model = model
            i_experiment = experiment
            i_member = member
        
        if row_idx != i:
            
            cube = iris.load_cube(datadir + filename)
            p87p5_trend = iris.load_cube(datadir + filename.replace('tailwidth','p87p5'))

            cube = iris.util.new_axis(cube)
            p87p5 = iris.util.new_axis(p87p5_trend)

            model_coord = iris.coords.DimCoord([row_idx],long_name='model',units='1')

            cube.add_dim_coord(model_coord, 0)
            p87p5.add_dim_coord(model_coord.copy(), 0)
            
            all_trend.append(cube)
            
    iris.util.unify_time_units(all_trend)
    iris.util.equalise_attributes(all_trend)
    
    trend_cube_list = iris.cube.CubeList(all_trend)
    trend_mm = trend_cube_list.concatenate_cube()
    trend_mmm = trend_mm.collapsed('model',iris.analysis.MEAN)
    trend_mmm = iris.util.squeeze(trend_mmm)

    # Guess bounds
    lon_coord = template.coord('longitude')
    if not lon_coord.has_bounds():
        lon_coord.guess_bounds()

    lat_coord = template.coord('latitude')
    if not lat_coord.has_bounds():
        lat_coord.guess_bounds()

    # Apply land mask
    masked_cube, _ = lat_land_mask(trend_mmm, template, lat_minimum, lat_maximum)
    chosen_masked, _ = lat_land_mask(chosen, template, lat_minimum, lat_maximum)
    chosen_p87p5_masked, _ = lat_land_mask(chosen_p87p5, template, lat_minimum, lat_maximum)
        
    x = masked_cube.data
    y = chosen_masked.data

    x = masked_cube.data.filled(np.nan)
    y = chosen_masked.data.filled(np.nan)

    # Filter out unphysical numbers arising from masking issues.
    mask = (
        np.isfinite(x) &
        np.isfinite(y) &
        (np.abs(x) < 5) &
        (np.abs(y) < 5)
    )

    x_valid = x[mask]
    y_valid = y[mask]

    # Correlation
    r = np.corrcoef(x_valid, y_valid)[0,1]
    print(f"Correlation: {r:.2f}")

    # Linear regression (slope + intercept)
    slope, intercept = np.polyfit(x_valid, y_valid, 1)
    print(f"Slope: {slope:.2f}, Intercept: {intercept:.2f}")

    plt.scatter(x_valid, y_valid, s=5)

    # Line of best fit
    x_fit = np.linspace(x_valid.min(), x_valid.max(), 100)
    y_fit = slope * x_fit + intercept
    
    all_slopes.append(slope)
    all_corr.append(r)
    #plt.plot(x_fit, y_fit, color='black', linewidth=2, label=f'Fit: y={slope:.2f}x + {intercept:.2f}')

    #plt.axhline(0, color='grey', linewidth=1)
    #plt.axvline(0, color='grey', linewidth=1)
    #plt.scatter(0, 0, color='red', s=40, zorder=5, label='(0,0)')

    #plt.title(f'HighResMIP multimodel mean\n and {i_model} {i_experiment} {i_member}')
    #plt.xlabel('HighResMIP multi-model mean tail width trend (°C/decade)')
    #plt.ylabel(f'{i_model} {i_experiment} {i_member} tail width trend (°C/decade)')

    #plt.legend()
    #save_figure(plt, f'HighResMIP_MMM_{i_model}_{i_experiment}_{i_member}_scatter_plot', 'Graphs')
    #plt.show()
    #plt.close()

all_slopes_array = np.array(all_slopes)
all_corr_array = np.array(all_corr)

plt.figure(figsize=(12, 6))
plt.hist(
    all_slopes_array, color="#4682B4",  
    alpha=0.85,
    edgecolor='white',
    linewidth=2
    , label='HRMIP members'
)

ax = plt.gca()
for spine in ax.spines.values():
    spine.set_linewidth(2)

plt.axvline(x = 0.67, color='black', linestyle='-', linewidth=2.5, label='ERA5')
plt.xlabel('Slope', fontsize=21)
plt.ylabel('Frequency', fontsize=21)
plt.xticks(fontsize=19)
plt.yticks(fontsize=19)
plt.grid(True, linestyle='--', alpha=0.5, linewidth=1.2)
plt.legend(fontsize=15, loc='upper right')
plt.title('d) Distribution of slopes from HRMIP members', fontsize=22, fontweight='bold')
plt.tight_layout()
save_figure(plt, 'Slopes_histogram', 'Graphs')
plt.show()

plt.figure(figsize=(12, 5))
plt.hist(
    all_corr_array, label = 'HighResMIP members'
)

plt.axvline(x = 0.19, color='black', linestyle='-', linewidth=2.5, label='ERA5')
plt.xlabel('Correlation')
plt.ylabel('Frequency')
plt.grid(True, linestyle='--', alpha=0.5, linewidth=1.2)
plt.legend(fontsize=13, loc='upper right')

plt.tight_layout()
save_figure(plt, 'Correlation_histogram', 'Graphs')
plt.show()

