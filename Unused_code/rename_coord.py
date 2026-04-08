import iris
import numpy as np

# Load your cube
cube = iris.load_cube("/data/users/samuel.sancar/Heat_extremes/Tailwidthslist/HadGEM3/HadGEM3-GC31-LL_hist-1950_r1i1p1f1_tailwidth_1976.nc")

# Check shape
print("Cube shape:", cube.shape)

# If the last dimension is longitude, create a DimCoord
n_lon = cube.shape[-1]
lon_points = np.linspace(-180, 180, n_lon)  # or use actual values if known
longitude_coord = iris.coords.DimCoord(lon_points, standard_name='longitude', units='degrees')

# Add the longitude coordinate to the cube
cube.add_dim_coord(longitude_coord, cube.ndim - 1)

# Verify
print([c.name() for c in cube.coords()])

# Save updated cube
iris.save(cube, "/data/users/samuel.sancar/Heat_extremes/Tailwidthslist/HadGEM3/HadGEM3-GC31-LL_hist-1950_r1i1p1f1_tailwidth_1976.nc")