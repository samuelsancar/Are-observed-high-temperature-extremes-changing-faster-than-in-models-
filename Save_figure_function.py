"""
Function to save figures with consistent format and dpi.
"""

import matplotlib.pyplot as plt
import os

def save_figure(plt, filename, location):
    folder = f"/home/users/samuel.sancar/Documents/Temperature_extremes/Plots/{location}/"
    os.makedirs(folder, exist_ok=True)  
    full_path = os.path.join(folder, f"{filename}.png")  

    plt.savefig(full_path, dpi=2000, format='png')
    print(f"Figure saved to: {full_path}")
