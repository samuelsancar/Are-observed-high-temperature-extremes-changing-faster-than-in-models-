"""
    Function for map plotting details:
    
    Detailed: More detailed map.
    [anything]: Simpler map details.

    
"""

import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.cm as mpl_cm

def map_features(type, ax, lat_min, lat_max):
    
    if type == 'Detail':
        
        ax.set_extent([-180, 180, lat_min, lat_max], crs=ccrs.PlateCarree(central_longitude=0))
        ax.coastlines(resolution='50m', linewidth=0.5)
        ax.add_feature(cfeature.NaturalEarthFeature(category='cultural', name='admin_0_boundary_lines_land', facecolor='none', scale='50m', edgecolor='0', lw=.05))
        ax.add_feature(cfeature.LAKES.with_scale('110m'), edgecolor='0', facecolor='none', lw=.05)
        ax.add_feature(cfeature.STATES.with_scale('110m'), edgecolor='0', facecolor='none', lw=.2)
        ax.add_feature(cfeature.RIVERS.with_scale('110m'), edgecolor='tab:blue', facecolor='none', lw=.4)
        ax.add_feature(cfeature.LAND, facecolor='lightgray')
        ax.add_feature(cfeature.OCEAN, facecolor='white', zorder =0)
        
        for spine in ax.spines.values():
            spine.set_linewidth(1.5)   
            spine.set_edgecolor('black')
        
        gl = ax.gridlines(draw_labels=False, linewidth=0.5, color='gray', alpha=0.5, linestyle='--')
        gl.top_labels = False
        gl.right_labels = False
            
    else:
    
        ax.set_extent([-180, 0, 60, 90], crs=ccrs.PlateCarree())
        ax.coastlines()
        ax.add_feature(cfeature.LAND, facecolor='white')
        ax.add_feature(cfeature.OCEAN, facecolor='white', zorder = 0)

        gl = ax.gridlines(draw_labels=True, linewidth=0.5, color='gray', alpha=0.5, linestyle='--')
        gl.top_labels = False
        gl.right_labels = False