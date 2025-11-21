#!/usr/bin/env python3
"""
Valens CLI Tool
===============

A post-processing tool for VASP outputs, designed to create publication-quality
plots with a modern aesthetic. Currently supports Density of States (DOS) plotting.

Features:
- Smart legend positioning
- Gradient fills for DOS peaks
- Custom font support
"""

import os
import sys
import argparse
import numpy as np
import matplotlib as mpl
mpl.use("agg")
mpl.rcParams["axes.unicode_minus"] = False
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.patches import Polygon
from matplotlib.ticker import AutoMinorLocator
from pymatgen.io.vasp.outputs import Vasprun

# ===============================================================
# Gradient fill aesthetic
# ===============================================================
def gradient_fill(x, y, ax=None, color=None, direction=1, **kwargs):
    """
    Fills the area under a curve with a vertical gradient.

    Args:
        x (array-like): X-axis data (Energy).
        y (array-like): Y-axis data (DOS).
        ax (matplotlib.axes.Axes, optional): The axes to plot on. Defaults to current axes.
        color (str, optional): The base color for the gradient.
        direction (int, optional): Direction of the gradient (unused currently).
        **kwargs: Additional arguments passed to ax.plot.

    Returns:
        matplotlib.lines.Line2D: The line object representing the curve.
    """
    if ax is None:
        ax = plt.gca()
    
    # Plot the main line
    line, = ax.plot(x, y, color=color, lw=2, **kwargs)
    
    # Determine fill color and alpha
    fill_color = line.get_color() if color is None else color
    alpha = line.get_alpha() or 1.0
    zorder = line.get_zorder()

    # Create a gradient image
    z = np.empty((100, 1, 4))
    rgb = mcolors.to_rgb(fill_color)
    z[:, :, :3] = rgb
    
    # Gradient alpha: starts from 0.3 to ensure visibility even for small peaks (e.g. CBM)
    z[:, :, -1] = np.linspace(0.3, alpha, 100)[:, None]
    
    xmin, xmax, ymin, ymax = x.min(), x.max(), 0, max(y.max(), 1e-6)

    # Display the gradient image
    im = ax.imshow(z, aspect="auto", extent=[xmin, xmax, ymin, ymax],
                   origin="lower", zorder=zorder)
    
    # Clip the gradient to the area under the curve
    xy = np.column_stack([x, y])
    xy = np.vstack([[xmin, ymin], xy, [xmax, ymin], [xmin, ymin]])
    clip = Polygon(xy, lw=0, facecolor="none", closed=True)
    ax.add_patch(clip)
    im.set_clip_path(clip)
    
    return line


# ===============================================================
# Data container
# ===============================================================
class ValensDos:
    """
    Simple container for Total DOS data.
    
    Attributes:
        energies (np.ndarray): Array of energy values (shifted by Fermi energy).
        total (np.ndarray): Total DOS values.
        efermi (float): Fermi energy.
    """
    def __init__(self, energies, total, efermi):
        self.energies = np.array(energies)
        self.total = np.array(total)
        self.efermi = float(efermi)


# ===============================================================
# Load DOS
# ===============================================================
def load_dos(vasprun, elements=None, **_):
    """
    Loads DOS data from a vasprun.xml file using pymatgen.

    Args:
        vasprun (str): Path to the vasprun.xml file or directory containing it.
        elements (list or dict, optional): Specific elements to extract PDOS for.

    Returns:
        tuple: (ValensDos object, dict of PDOS data)
    """
    print("🔍 Reading vasprun.xml ...")
    
    # Handle directory input
    if os.path.isdir(vasprun):
        vasprun = os.path.join(vasprun, "vasprun.xml")
    
    if not os.path.exists(vasprun):
        raise FileNotFoundError(f"{vasprun} not found")

    # Parse VASP output
    vr = Vasprun(vasprun)
    dos = vr.complete_dos
    efermi = getattr(dos, "efermi", vr.efermi)
    
    # Shift energies to set Fermi level at 0
    energies = dos.energies - efermi

    # Calculate total DOS (sum of spins)
    total = np.zeros_like(energies)
    for spin in dos.densities:
        total += dos.densities[spin]

    # Extract Projected DOS
    pdos = get_pdos(dos, elements)
    
    return ValensDos(energies, total, efermi), pdos


# ===============================================================
# Extract PDOS
# ===============================================================
def get_pdos(dos, elements=None):
    """
    Extracts Projected DOS (PDOS) for specified elements.

    Args:
        dos (pymatgen.electronic_structure.dos.CompleteDos): The complete DOS object.
        elements (list or dict, optional): Elements to extract. If None, extracts all.

    Returns:
        dict: A dictionary where keys are element symbols and values are dicts of orbital DOS.
    """
    structure = dos.structure
    symbols = [str(site.specie) for site in structure]
    
    # If no elements specified, use all unique elements in the structure
    if not elements:
        unique = sorted(set(symbols))
        elements = {el: () for el in unique}
    else:
        # Ensure elements is a dict if passed as list
        if isinstance(elements, list):
             elements = {el: () for el in elements}

    pdos = {}
    for el in elements:
        # Find all sites corresponding to this element
        el_sites = [s for s in structure if str(s.specie) == el]
        el_pdos = {}
        
        for site in el_sites:
            try:
                site_dos = dos.get_site_spd_dos(site)
            except Exception:
                continue
            
            # Sum up contributions from orbitals (s, p, d, f)
            for orb, orb_dos in site_dos.items():
                label = orb.name[0] # e.g., 's', 'p', 'd'
                dens = np.zeros_like(dos.energies)
                for spin in orb_dos.densities:
                    dens += orb_dos.densities[spin]
                el_pdos[label] = el_pdos.get(label, 0) + dens
        
        pdos[el] = el_pdos
    return pdos


# ===============================================================
# Plotting with smart legend & font control (Valens theme)
# ===============================================================
def plot_dos(dos, pdos, out="valens_dos.png",
             xlim=(-6, 6), ylim=None, figsize=(5, 4),
             dpi=400, legend_loc="auto", font="Arial"):
    """
    Plots the Total and Projected DOS with the Valens visual style.

    Args:
        dos (ValensDos): The total DOS data.
        pdos (dict): The projected DOS data.
        out (str): Output filename.
        xlim (tuple): Energy range (min, max).
        ylim (tuple, optional): DOS range (min, max).
        figsize (tuple): Figure size in inches.
        dpi (int): Resolution of the output image.
        legend_loc (str): Legend location strategy.
        font (str): Font family to use.
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
    mpl.rcParams["font.size"] = 10

    plt.style.use("default")
    fig, ax = plt.subplots(figsize=figsize)
    
    # Fermi level line
    ax.axvline(0, color="k", lw=0.8, ls="--", alpha=0.7)

    # Color palette for elements
    palette = ["#4b0082", "#e63946", "#2a9d8f", "#ffb703", "#6a994e", "#8e44ad", "#118ab2"]
    lines, labels = [], []

    # Plot PDOS for each element
    for i, (el, el_pdos) in enumerate(pdos.items()):
        c = palette[i % len(palette)]
        total_el = sum(el_pdos.values())
        
        # Apply gradient fill
        gradient_fill(dos.energies, total_el, ax=ax, color=c, alpha=0.9)
        
        line, = ax.plot(dos.energies, total_el, lw=1.5, color=c, label=el)
        lines.append(line)
        labels.append(el)

    # Plot Total DOS
    total_line, = ax.plot(dos.energies, dos.total, color="black", lw=1.4, label="Total DOS")
    lines.append(total_line)
    labels.append("Total DOS")

    # Axis settings
    ax.set_xlim(*xlim)
    if ylim:
        ax.set_ylim(*ylim)
    ax.set_xlabel("Energy (eV)", fontsize=12, weight="bold", labelpad=6)
    ax.set_ylabel("Density of States", fontsize=12, weight="bold", labelpad=6)
    ax.set_xticks(np.linspace(xlim[0], xlim[1], 5))
    ax.set_yticks([])

    # --- Smart legend visibility ---
    # Only show legend if PDOS peaks are significant enough in the visible range
    show_legend = False
    
    # Determine mask for visible x-range
    x_mask = (dos.energies >= xlim[0]) & (dos.energies <= xlim[1])
    
    if ylim:
        y_threshold = 0.10 * ylim[1]
        for el_pdos in pdos.values():
            total_el = sum(el_pdos.values())[x_mask]
            if np.max(total_el) >= y_threshold:
                show_legend = True
                break
    else:
        # If no ylim provided, default to showing legend if there is data
        if np.any(x_mask):
             show_legend = True

    if show_legend:
        # Check for overlap to decide legend position
        if np.any(x_mask):
            data_ymax = np.max(dos.total[x_mask])
            overlap = data_ymax > 0.8 * (ylim[1] if ylim else data_ymax)
        else:
            overlap = False
            
        n_items = len(labels)
        if overlap or n_items > 4:
            loc, ncol = "best", 2
        else:
            loc, ncol = "upper right", 1

        legend = ax.legend(
            lines, labels,
            frameon=False,
            fontsize=10,
            loc=loc,
            ncol=ncol,
            handlelength=1.5,
            columnspacing=0.8,
            handletextpad=0.6,
        )
        for text in legend.get_texts():
            text.set_fontweight("bold")
    else:
        print("ℹ️ Legend hidden (no PDOS exceeds 10% of visible range).")

    ax.xaxis.set_minor_locator(AutoMinorLocator(2))
    plt.tight_layout(pad=0.4)
    plt.savefig(out, dpi=dpi)
    plt.close(fig)
    print(f"📊 DOS plot saved as {out} with font: {font}")


# ===============================================================
# Main CLI
# ===============================================================
def main():
    """
    Main entry point for the CLI.
    """
    parser = argparse.ArgumentParser(description="Valens: VASP Post-Processing Tool")
    subparsers = parser.add_subparsers(dest="command", help="Available commands")

    # --- DOS Subcommand ---
    dos_parser = subparsers.add_parser("dos", help="Plot Density of States (DOS)")
    dos_parser.add_argument("vasprun", nargs="?", default=".", help="Path to vasprun.xml or directory containing it")
    dos_parser.add_argument("-e", "--elements", nargs="+", help="Specific elements to plot PDOS for")
    dos_parser.add_argument("--xlim", nargs=2, type=float, default=[-6, 6], help="Energy range (min max)")
    dos_parser.add_argument("--ylim", nargs=2, type=float, help="DOS range (min max)")
    dos_parser.add_argument("--scale", type=float, default=1.0, help="Scaling factor for Y-axis (DOS density)")
    dos_parser.add_argument("-o", "--output", default="valens_dos.png", help="Output filename")
    dos_parser.add_argument("--font", default="Arial", help="Font family")

    args = parser.parse_args()

    if args.command == "dos":
        try:
            dos_data, pdos_data = load_dos(args.vasprun, elements=args.elements)
            
            # Apply scaling
            if args.scale != 1.0:
                print(f"⚖️ Scaling DOS by factor of {args.scale}")
                dos_data.total /= args.scale
                for el in pdos_data:
                    for orb in pdos_data[el]:
                        pdos_data[el][orb] /= args.scale

            plot_dos(dos_data, pdos_data, out=args.output,
                     xlim=tuple(args.xlim), ylim=tuple(args.ylim) if args.ylim else None,
                     font=args.font)
        except Exception as e:
            print(f"❌ Error: {e}")
            sys.exit(1)
    else:
        parser.print_help()

if __name__ == "__main__":
    main()
