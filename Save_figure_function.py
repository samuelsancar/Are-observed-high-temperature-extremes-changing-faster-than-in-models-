"""
Function to save figures with consistent styling and formatting.
"""

import matplotlib.pyplot as plt
import os

def save_figure(plt, filename, location):
    folder = f"/home/users/samuel.sancar/Documents/Heat_Extremes/Plots/{location}/"
    os.makedirs(folder, exist_ok=True)  
    full_path = os.path.join(folder, f"{filename}.png")  

    plt.savefig(full_path, dpi=500, format='png')
    print(f"Figure saved to: {full_path}")
