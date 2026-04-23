"""" 
Cube inspector:
"""

import iris

def cube_inspector(c, name):
    print(f"\n=== {name} summary ===")
    print(c.summary())
    print("Dim coords and dims:", [(coord.name(), dims) for coord, dims in 
                                  ((c.coord(cn), c.coord_dims(cn)) for cn in [co.name() for co in c.dim_coords])])
    print("All coords:", [(coord.name(), c.coord_dims(coord)) for coord in c.coords()])
    print("Shape:", c.shape)
    # Check time details
    if c.coords('time'):
        t = c.coord('time')
        print("time points:", t.points[:3], "...", t.points[-3:])
        print("time units:", t.units)
        print("time has bounds:", t.has_bounds())