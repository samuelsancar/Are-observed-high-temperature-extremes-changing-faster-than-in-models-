"""

Script to calculate tailwidths of daily maximum 2m temperature from JRA-55 reanalysis data.

To run: python JRA-55_tailwidths.py

"""

import iris
import numpy as np
import iris.coord_systems as cs
import os

datadir = '/data/users/samuel.sancar/Heat_extremes/JRA-55_tailwidths/daily_maximum__2m_temperature/'
outdir = '/data/users/samuel.sancar/Heat_extremes/JRA-55_tailwidths/Yearly_tailwidths/'
template = iris.load_cube('template_cube_1deg.nc')

reanalysis = 'jra55'
variable = '_minmax_surf.015_tmax.reg_tl319.'

os.makedirs(outdir, exist_ok=True)

# Function to match coordinate systems and regrid
def match_and_regrid(tailwidth, template):
    def assign_default_coord_system(cube):
        default_cs = cs.GeogCS(6371229.0)
        for coord_name in ['latitude', 'longitude']:
            coord = cube.coord(coord_name)
            if coord.coord_system is None:
                coord.coord_system = default_cs

    assign_default_coord_system(tailwidth)
    assign_default_coord_system(template)

    for coord_name in ['latitude', 'longitude']:
        merged_coord = tailwidth.coord(coord_name)
        template_coord = template.coord(coord_name)
        if merged_coord.units != template_coord.units:
            merged_coord.units = template_coord.units

    return tailwidth.regrid(template, iris.analysis.Linear())

# Calculate tailwidths over years t1-t2
for y in np.arange(1958,2023):

    fn = reanalysis+variable+'_'+str(y)+'.nc'
    year_cubes = iris.load_cube(datadir+fn)

    #Calculate p99-p87.5 for this year:
    p99 = year_cubes.collapsed('time',iris.analysis.PERCENTILE,percent=99)
    p87p5 =year_cubes.collapsed('time',iris.analysis.PERCENTILE,percent=87.5)
    tail_width = p99-p87p5

    #regrid onto template:
    tail_width_rg = match_and_regrid(tail_width, template)
    p87p5_regrid = match_and_regrid(p87p5, template)

    out_file = os.path.join(outdir, f"{reanalysis}_{variable}_{y}_tail_width.nc")
    outfile_87p5 = os.path.join(outdir, f"{reanalysis}_{variable}_{y}_87p5.nc")
    iris.save(tail_width_rg, out_file)
    iris.save(p87p5_regrid, outfile_87p5)
    
    
    
    
# Extract trend data from masked cube
trend_data = model_mask.data

# Convert significance mask to p-values array
p_values = np.array(statsig_ds)

# Avoid log(0) by clipping p-values
p_values = np.clip(p_values, 1e-10, 1)

# Compute significance score
significance_score = -np.log10(p_values)

# Replace NaNs with zeros
trend_data = np.nan_to_num(trend_data, nan=0.0)
significance_score = np.nan_to_num(significance_score, nan=0.0)

# Compute weighted index
weighted_index = np.abs(trend_data) * significance_score

# Create Iris cube 
weighted_cube = model_mask.copy()
if weighted_index.shape != weighted_cube.shape:
    raise ValueError(f"Shape mismatch: weighted_index {weighted_index.shape} vs weighted_cube {weighted_cube.shape}")
weighted_cube.data = weighted_index

# Compute colour scale
try:
    vmax = np.nanpercentile(weighted_index, 99)
    if not np.isfinite(vmax):
        vmax = 1.0
except Exception:
    vmax = 1.0

# Plot Weighted Significance-Magnitude Index 
fig2 = plt.figure(figsize=(12, 6))
ax2 = plt.axes(projection=ccrs.Robinson(central_longitude=0))
ax2.set_extent([-180, 180, -60, 90], crs=ccrs.PlateCarree())

# Mask weighted index
weighted_masked, _ = lat_land_mask(weighted_cube, lat_minimum, lat_maximum)

iplt.pcolormesh(weighted_masked, cmap='coolwarm', vmin=0, vmax=vmax)

# Plot stippling 
if stipple_lons.size > 0 and stipple_lats.size > 0:
    ax2.scatter(stipple_lons, stipple_lats, color='black', s=0.4, alpha=0.5,
                transform=ccrs.PlateCarree(), label='Significant')

ax2.coastlines(resolution='50m', linewidth=0.5)
ax2.add_feature(cfeature.LAND, facecolor='lightgray')
ax2.add_feature(cfeature.OCEAN, facecolor='white')

# Colourbar
cbar2 = plt.colorbar(ax2.collections[0], orientation='horizontal', pad=0.05, aspect=50)
cbar2.set_label('Weighted Significance-Magnitude Index', fontsize=12)

plt.title(f"Weighted Significance-Magnitude Index for \n ERA5 tailwidth trends", fontsize=14, weight='bold')

plt.tight_layout()
save_figure(plt, 'ERA5_tailwidth_trends_weighted_index', 'Maps')
plt.show()


