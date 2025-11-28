#!/usr/bin/env python3
"""
DOS Plotting Module
===================

Handles Density of States (DOS) plotting with gradient fills and smart legend.
"""

import os
import numpy as np
import matplotlib as mpl
mpl.use("agg")
mpl.rcParams["axes.unicode_minus"] = False
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.patches import Polygon
from matplotlib.ticker import AutoMinorLocator
from pymatgen.io.vasp import Vasprun


# ===============================================================
# Gradient fill aesthetic
# ===============================================================
def gradient_fill(x, y, ax=None, color=None, xlim=None, **kwargs):
    """
    Fills the area under a curve with a vertical gradient.

    Args:
        x (array-like): X-axis data (Energy).
        y (array-like): Y-axis data (DOS).
        ax (matplotlib.axes.Axes, optional): The axes to plot on. Defaults to current axes.
        color (str, optional): The base color for the gradient.
        xlim (tuple, optional): X-axis limits to restrict gradient fill.
        **kwargs: Additional arguments passed to ax.plot.

    Returns:
        matplotlib.lines.Line2D: The line object representing the curve.
    """
    if ax is None:
        ax = plt.gca()
    
    # Don't filter by xlim - use full data range for better appearance
    if len(x) == 0 or len(y) == 0:
        return None
    
    # Plot the main line
    line, = ax.plot(x, y, color=color, lw=2, **kwargs)
    
    # Determine fill color and alpha
    fill_color = line.get_color() if color is None else color
    alpha = line.get_alpha() or 1.0
    zorder = line.get_zorder()

    # Create a gradient image with more aggressive alpha
    z = np.empty((100, 1, 4))
    rgb = mcolors.to_rgb(fill_color)
    z[:, :, :3] = rgb
    
    
    # Gradient: transparent at bottom (y=0), opaque near the curve
    # This creates a gradient from bottom to top of the filled area
    z[:, :, -1] = np.linspace(0.05, 0.75, 100)[:, None]
    
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
    
    # Handle directory input
    if os.path.isdir(vasprun):
        vasprun = os.path.join(vasprun, "vasprun.xml")
    
    if not os.path.exists(vasprun):
        raise FileNotFoundError(f"{vasprun} not found")

    # Parse VASP output
    vr = Vasprun(vasprun)
    dos = vr.complete_dos
    efermi = getattr(dos, "efermi", vr.efermi)
    
    # Determine reference energy
    # For insulators/semiconductors, align VBM to 0 eV
    # For metals, align Fermi level to 0 eV
    try:
        # Try using BandStructure first (more robust)
        bs = vr.get_band_structure()
        if not bs.is_metal():
            efermi = bs.get_vbm()["energy"]
    except Exception:
        # Fallback to DOS-based detection
        try:
            cbm, vbm = dos.get_cbm_vbm()
            if cbm - vbm > 0.01:  # Band gap detected
                efermi = vbm
        except Exception:
            pass
    
    # Shift energies to set reference at 0
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
             dpi=400, legend_loc="auto", font="Arial",
             show_fermi=False, show_total=True, plotting_config=None,
             legend_cutoff=0.10, scale_factor=1.0):
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
        show_fermi (bool): Whether to draw a dashed line at the Fermi level (E=0).
        show_total (bool): Whether to plot the Total DOS.
        plotting_config (list): List of (Element, Orbital) tuples to plot.
        legend_cutoff (float): Threshold (as fraction) for showing items in legend (default: 0.10).
        scale_factor (float): Factor to scale the Y-axis limits (zoom in).
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

    plt.style.use("default")
    fig, ax = plt.subplots(figsize=figsize)
    
    # Fermi level line (optional)
    if show_fermi:
        ax.axvline(0, color="k", lw=0.8, ls="--", alpha=0.7)

    # Color palette for elements
    # Expanded color palette for better distinction
    # Reordered to maximize contrast between consecutive items
    palette = [
        "#4b0082", # Indigo
        "#e63946", # Red
        "#2a9d8f", # Teal
        "#ffb703", # Yellow
        "#0077b6", # Blue
        "#8e44ad", # Purple
        "#d62828", # Dark Red
        "#118ab2", # Light Blue
        "#f4a261", # Orange
        "#003049", # Dark Blue
        "#6a994e", # Green
        "#023e8a", # Royal Blue
        "#0096c7", # Cyan
        "#00b4d8", # Sky Blue
        "#48cae4", # Light Cyan
        "#90e0ef", # Pale Blue
        "#ade8f4", # Very Pale Blue
        "#caf0f8"  # White Blue
    ]
    lines, labels = [], []
    
    # Determine mask for visible x-range (used for scaling and legend)
    x_mask = (dos.energies >= xlim[0]) & (dos.energies <= xlim[1])

    # Determine what to plot
    if plotting_config:
        items_to_plot = plotting_config
    else:
        # Default: Plot all orbitals for each loaded element
        items_to_plot = []
        for el, el_pdos in pdos.items():
            for orb in el_pdos.keys():
                items_to_plot.append((el, orb))

    # Plot PDOS
    max_visible_y = 0  # Track maximum Y value in visible range
    
    for i, (el, orb) in enumerate(items_to_plot):
        if el not in pdos:
            continue
            
        # Assign unique color for each orbital contribution
        c = palette[i % len(palette)]
        
        if orb == 'total':
            y_data = sum(pdos[el].values())
            label = el
        else:
            if orb in pdos[el]:
                y_data = pdos[el][orb]
                label = f"{el}({orb})"
            else:
                continue
        
        # Check contribution in visible range
        visible_y = y_data[x_mask]
        if len(visible_y) > 0:
            max_y = np.max(visible_y)
            max_visible_y = max(max_visible_y, max_y)
            
            # Store for later threshold check
            line, = ax.plot(dos.energies, y_data, lw=1.5, color=c, label=label, alpha=0)
            lines.append((line, y_data, max_y, c, label))
    
    # Calculate threshold (legend_cutoff of max visible)
    threshold = legend_cutoff * max_visible_y
    
    # Filter legend items but keep all plot lines
    final_lines = []
    final_labels = []
    
    for line, y_data, max_y, c, label in lines:
        # Always plot the line (make it visible)
        line.set_alpha(1.0)
        
        # Always apply gradient fill for visible lines
        gradient_fill(dos.energies, y_data, ax=ax, color=c, alpha=0.9)
        
        # Only add to legend if above main threshold
        if max_y >= threshold:
            final_lines.append(line)
            final_labels.append(label)
    
    # Update lines and labels for legend creation
    lines = final_lines
    labels = final_labels

    # Plot Total DOS
    if show_total:
        ax.plot(dos.energies, dos.total, color="k", lw=1.2, label="Total DOS")
        gradient_fill(dos.energies, dos.total, ax=ax, color="k", alpha=0.15)
    
    # Auto-scale Y-axis based on visible range if ylim not provided
    if not ylim and len(lines) > 0:
        max_y_in_range = 0
        for line in lines:
            y_data = line.get_ydata()
            visible_y = y_data[x_mask]
            if len(visible_y) > 0:
                max_y_in_range = max(max_y_in_range, np.max(visible_y))
        
        if show_total:
            visible_total = dos.total[x_mask]
            if len(visible_total) > 0:
                max_y_in_range = max(max_y_in_range, np.max(visible_total))
        
        if max_y_in_range > 0:
            # Apply scaling factor to the limit (zoom in)
            limit = (max_y_in_range * 1.1) / scale_factor
            ax.set_ylim(0, limit)

    # Axis settings
    ax.set_xlim(*xlim)
    if ylim:
        ax.set_ylim(*ylim)
    ax.set_xlabel("Energy (eV)", fontsize=14, weight="bold", labelpad=6)
    ax.set_ylabel("Density of States", fontsize=14, weight="bold", labelpad=6)
    
    # Set x-ticks with 1 eV spacing
    xticks = np.arange(np.ceil(xlim[0]), np.floor(xlim[1]) + 1, 1)
    ax.set_xticks(xticks)
    # Format tick labels: show integers without .0
    tick_labels = [f'{int(x)}' if x == int(x) else f'{x}' for x in xticks]
    ax.set_xticklabels(tick_labels, fontweight="bold")
    ax.set_yticks([])

    # --- Smart legend visibility ---
    # Only show legend if there are items to display
    show_legend = len(lines) > 0

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
            fontsize=13,
            loc=loc,
            ncol=ncol,
            handlelength=1.5,
            columnspacing=0.8,
            handletextpad=0.6,
        )
        for text in legend.get_texts():
            text.set_fontweight("bold")

    ax.xaxis.set_minor_locator(AutoMinorLocator(2))
    plt.tight_layout(pad=0.4)
    plt.savefig(out, dpi=dpi)
    plt.close(fig)
