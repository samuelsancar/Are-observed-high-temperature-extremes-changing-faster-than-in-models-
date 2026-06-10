"""
10/04/2026
----------

Code to plot the HighResMIP mean tail width trend (calculation of mmm is different to HighResMIP_mmm_trend_plotter.py).
Produces figure S1a in supplementary material.

To run: python HighResMIP_mmm_trend_plotter_2.0.py
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
from Map_features import map_features
from Save_figure_function import save_figure

cmap_rb = mpl_cm.get_cmap('coolwarm', 15)

# Beginning and end years
t1 = 1958
t2 = 2023

lat_minimum = -60
lat_maximum = 90

years = np.arange(t1, t2)

df = pd.read_excel("/data/users/samuel.sancar/Model_table.xlsx", sheet_name=0, engine="openpyxl")

# ensure consistent column names (strip spaces)
df.columns = df.columns.str.strip()

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

    template = iris.load_cube('template_cube_1deg.nc')
    datadir = f'/data/users/samuel.sancar/Heat_extremes/Trends/{Mode}_model_trends/'
    filename = f"{model}_{experiment}_{member}_tailwidth_trend_1958-2022.nc"
    
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

# Extract grid point locations for stippling
lon = cube.coord('longitude').points
lat = cube.coord('latitude').points
lon2d, lat2d = np.meshgrid(lon, lat)

# Plot trends
fig = plt.figure(figsize=(10, 6))
ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=0))

map_features('Detail', ax, lat_minimum, lat_maximum)

map = iplt.pcolormesh(masked_cube, cmap=cmap_rb, vmin=-0.15, vmax=0.15)

#cbar = plt.colorbar(map, ax = ax, orientation='horizontal', pad=0.05, aspect=50)
#cbar.ax.tick_params(labelsize=12)
#cbar.set_label('Trend in Tailwidth (°C per decade)', fontsize=14)
    
plt.title(f"a) HighResMIP MMM tail width trends", fontsize=17, weight='bold')
plt.tight_layout()
save_figure(plt, f'HighResMIP_MMM_tailwidths_2.0', 'Maps')
plt.show()