"""
    Function to calculate trend using pearson r coefficient.

    returns: change per decade

"""

import iris

def calculate_trend(year_cube, data_cube):
    r = iris.analysis.stats.pearsonr(year_cube, data_cube, corr_coords='time')
    std_data = data_cube.collapsed('time', iris.analysis.STD_DEV)
    std_year = year_cube.collapsed('time', iris.analysis.STD_DEV)
    return r * std_data / std_year * 10  # decadal trend