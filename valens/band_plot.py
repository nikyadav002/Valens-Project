import os
import numpy as np
import matplotlib as mpl
mpl.use("agg")
mpl.rcParams["axes.unicode_minus"] = False
import matplotlib.pyplot as plt
from pymatgen.io.vasp import Vasprun, BSVasprun
from pymatgen.electronic_structure.plotter import BSPlotter

def plot_band_structure(vasprun_path, kpoints_path=None, output="band_structure.png",
                        ylim=None, figsize=(6, 5), dpi=400, font="Arial"):
    """
    Plots the electronic band structure from vasprun.xml.
    
    Args:
        vasprun_path (str): Path to vasprun.xml file.
        kpoints_path (str, optional): Path to KPOINTS file (for labels).
        output (str): Output filename.
        ylim (tuple, optional): Energy range (min, max).
        figsize (tuple): Figure size in inches.
        dpi (int): Resolution of the output image.
        font (str): Font family.
    """
    
    # --- Font configuration ---
    font_map = {
        "arial": "Arial",
        "helvetica": "Helvetica",
        "times": "Times New Roman",
        "times new roman": "Times New Roman",
    }
    font = font_map.get(font.lower(), "Arial")
    mpl.rcParams["font.family"] = font
    mpl.rcParams["axes.linewidth"] = 1.4
    mpl.rcParams["font.weight"] = "bold"
    mpl.rcParams["font.size"] = 12
    mpl.rcParams["xtick.major.width"] = 1.2
    mpl.rcParams["ytick.major.width"] = 1.2

    print(f"🔍 Reading {vasprun_path} ...")
    
    try:
        # Load VASP output
        # BSVasprun is optimized for band structures
        vr = BSVasprun(vasprun_path, parse_projected_eigen=False)
        bs = vr.get_band_structure(kpoints_filename=kpoints_path, line_mode=True)
    except Exception as e:
        raise ValueError(f"Failed to load band structure: {e}")

    # Use BSPlotter to get the data in a plot-friendly format
    bs_plotter = BSPlotter(bs)
    data = bs_plotter.bs_plot_data(zero_to_efermi=True)
    
    # Extract data
    distances = data['distances'] # List of lists (one per segment)
    energies = data['energy']     # List of dicts (one per segment), keys are '1' (up) and '-1' (down)
    ticks = data['ticks']         # Dict with 'distance' and 'label'
    vbm = data['vbm']             # VBM energy (should be 0 after zero_to_efermi=True)
    
    # Setup plot
    fig, ax = plt.subplots(figsize=figsize)
    
    # Colors
    color_vb = "#8e44ad" # Purple
    color_cb = "#2a9d8f" # Teal
    
    # Plot bands
    # Iterate over segments
    for i in range(len(distances)):
        d = distances[i]
        
        # Iterate over spin channels
        for spin in energies[i]:
            # energies[i][spin] is a list of arrays (one per band)
            for band in energies[i][spin]:
                # Determine color based on energy relative to VBM (0 eV)
                # We check the mean energy of the band segment
                if np.mean(band) <= 0:
                    c = color_vb
                else:
                    c = color_cb
                
                ax.plot(d, band, color=c, lw=1.5, alpha=1.0)

    # Setup X-axis (K-path)
    ax.set_xticks(ticks['distance'])
    # Clean up labels (remove formatting like $ if needed, but pymatgen usually does a good job)
    clean_labels = [l.replace("$\\mid$", "|") for l in ticks['label']]
    ax.set_xticklabels(clean_labels, fontsize=12, fontweight="bold")
    
    # Draw vertical lines at high-symmetry points
    for d in ticks['distance']:
        ax.axvline(d, color="k", lw=0.8, ls="-", alpha=0.3)
        
    # Draw VBM line (E=0)
    ax.axhline(0, color="k", lw=0.8, ls="--", alpha=0.5)

    # Setup Y-axis
    ax.set_ylabel("Energy (eV)", fontsize=14, fontweight="bold", labelpad=8)
    if ylim:
        ax.set_ylim(ylim)
    else:
        # Default zoom around gap
        ax.set_ylim(-4, 4)
        
    ax.set_xlim(distances[0][0], distances[-1][-1])

    plt.tight_layout()
    plt.savefig(output, dpi=dpi)
    plt.close(fig)
    print(f"✅ Band structure saved to {output}")
