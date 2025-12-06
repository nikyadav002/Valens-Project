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
from pymatgen.electronic_structure.core import Spin


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
    
    # Gradient Logic based on relative height
    # We want opacity to be proportional to height relative to the max visible value (ymax_ref)
        
    # Create normalized alpha gradient (0 to 1)
    # We map y-values to alpha values. 
    # Since imshow fills a rectangle, we create a vertical gradient 
    # and clip it later.
    
    # Opacity range: 0.05 (at axis) to 0.95 (at max visible height)
    # This ensures "Darker colour at the top"
    min_alpha = 0.05
    max_alpha = 0.95
    
    # Create the gradient array (vertical)
    # 0 is bottom, 1 is top
    gradient_vector = np.linspace(min_alpha, max_alpha, 100)
    
    # IMPORTANT: Restore alpha scaling so total DOS (alpha=0.15) stays faint
    gradient_vector *= alpha
    
    # If data is negative (Spin Down), we want opaque at bottom (peak) and transparent at top (axis)
    if np.mean(y) < 0:
        gradient_vector = gradient_vector[::-1]
        
    z[:, :, -1] = gradient_vector[:, None]
    
    xmin, xmax = x.min(), x.max()
    
    # Determine extent based on LOCAL curve limits
    # User requested "each curve should have its own gradient"
    # So we scale from 0 to curve.max()
    
    local_ymax = max(y.max(), abs(y.min()))
    if local_ymax == 0: local_ymax = 1.0 # Prevent singular extent
    
    if np.mean(y) < 0:
        # Spin Down: Extent from -local_ymax to 0
        extent_ymin = -local_ymax
        extent_ymax = 0
    else:
        # Spin Up: Extent from 0 to local_ymax
        extent_ymin = 0
        extent_ymax = local_ymax
    
    # Display the gradient image
    # We use the explicit extent calculated above
    im = ax.imshow(z, aspect="auto", extent=[xmin, xmax, extent_ymin, extent_ymax],
                   origin="lower", zorder=zorder)
    
    # Clip the gradient to the area under the curve
    # We need to close the polygon at y=0
    xy = np.column_stack([x, y])
    
    # Construct polygon vertices:
    # Start at (xmin, 0), go along curve (x, y), end at (xmax, 0), close back to start
    verts = np.vstack([[x[0], 0], xy, [x[-1], 0], [x[0], 0]])
    
    clip = Polygon(verts, lw=0, facecolor="none", closed=True)
    ax.add_patch(clip)
    im.set_clip_path(clip)
    
    return line


# ===============================================================
# Data container
# ===============================================================
class ValyteDos:
    """
    Container for Total DOS data.
    
    Attributes:
        energies (np.ndarray): Array of energy values (shifted by Fermi energy).
        densities (dict): Dictionary of {Spin: np.ndarray} for total DOS.
        efermi (float): Fermi energy.
    """
    def __init__(self, energies, densities, efermi):
        self.energies = np.array(energies)
        self.densities = densities
        self.efermi = float(efermi)
        
    @property
    def total(self):
        """Returns the sum of all spin channels."""
        tot = np.zeros_like(self.energies)
        for spin in self.densities:
            tot += self.densities[spin]
        return tot

    @property
    def spin_up(self):
        return self.densities.get(Spin.up, np.zeros_like(self.energies))

    @property
    def spin_down(self):
        return self.densities.get(Spin.down, np.zeros_like(self.energies))


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
        tuple: (ValyteDos object, dict of PDOS data)
    """
    
    # Handle directory input
    if os.path.isdir(vasprun):
        vasprun = os.path.join(vasprun, "vasprun.xml")
    
    if not os.path.exists(vasprun):
        raise FileNotFoundError(f"{vasprun} not found")

    # Parse VASP output
    vr = Vasprun(vasprun)
    dos = vr.complete_dos
    
    # Get Fermi Energy
    efermi = dos.efermi
    
    # Attempt to align VBM to 0 for insulators/semiconductors
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

    # Extract Projected DOS
    pdos = get_pdos(dos, elements)
    
    return ValyteDos(energies, dos.densities, efermi), pdos


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
              pdos[element][orbital] = {Spin.up: array, Spin.down: array}
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
                
                # Initialize dictionaries for spins if not exists
                if label not in el_pdos:
                    el_pdos[label] = {}
                
                for spin in orb_dos.densities:
                    if spin not in el_pdos[label]:
                        el_pdos[label][spin] = np.zeros_like(dos.energies)
                    el_pdos[label][spin] += orb_dos.densities[spin]
        
        pdos[el] = el_pdos
    return pdos


# ===============================================================
# Plotting with smart legend & font control (Valyte theme)
# ===============================================================
def plot_dos(dos, pdos, out="valyte_dos.png",
             xlim=(-6, 6), ylim=None, figsize=(5, 4),
             dpi=400, legend_loc="auto", font="Arial",
             show_fermi=False, show_total=True, plotting_config=None,
             legend_cutoff=0.10, scale_factor=1.0):
    """
    Plots the Total and Projected DOS with the Valyte visual style.

    Args:
        dos (ValyteDos): The total DOS data.
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
    
    # Check if spin-polarized
    is_spin_polarized = Spin.down in dos.densities
    
    # Fermi level line (optional)
    if show_fermi:
        ax.axvline(0, color="k", lw=0.8, ls="--", alpha=0.7)
        
    # Zero line for spin polarized plots
    if is_spin_polarized:
        ax.axhline(0, color="k", lw=0.5, alpha=1.0)

    # Color palette for elements
    # Expanded color palette for better distinction
    # Reordered to maximize contrast between consecutive items
    palette = [
        "#4b0082", # Indigo
        "#0096c7", # Cyan
        "#e63946", # Red
        "#023e8a", # Royal Blue
        "#ffb703", # Yellow
        "#2a9d8f", # Teal
        "#8e44ad", # Purple
        "#118ab2", # Light Blue
        "#d62828", # Dark Red
        "#00b4d8", # Sky Blue
        "#f4a261", # Orange
        "#003049", # Dark Blue
        "#6a994e", # Green
        "#48cae4", # Light Cyan
        "#0077b6", # Blue
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
    min_visible_y = 0  # Track minimum Y value (for spin down)
    
    for i, (el, orb) in enumerate(items_to_plot):
        if el not in pdos:
            continue
            
        # Assign unique color for each orbital contribution
        c = palette[i % len(palette)]
        
        # Prepare data for plotting
        # y_data_up and y_data_down
        
        if orb == 'total':
            # Sum all orbitals for this element
            y_up = np.zeros_like(dos.energies)
            y_down = np.zeros_like(dos.energies)
            
            for o_data in pdos[el].values():
                y_up += o_data.get(Spin.up, np.zeros_like(dos.energies))
                y_down += o_data.get(Spin.down, np.zeros_like(dos.energies))
                
            label = el
        else:
            if orb in pdos[el]:
                y_up = pdos[el][orb].get(Spin.up, np.zeros_like(dos.energies))
                y_down = pdos[el][orb].get(Spin.down, np.zeros_like(dos.energies))
                label = f"{el}({orb})"
            else:
                continue
        
        # Invert spin down
        y_down = -y_down
        
        # Check contribution in visible range
        visible_y_up = y_up[x_mask]
        visible_y_down = y_down[x_mask]
        
        has_visible_data = False
        current_max_y = 0
        
        if len(visible_y_up) > 0:
            max_y = np.max(visible_y_up)
            max_visible_y = max(max_visible_y, max_y)
            current_max_y = max(current_max_y, max_y)
            if max_y > 1e-6: has_visible_data = True
            
        if is_spin_polarized and len(visible_y_down) > 0:
            min_y = np.min(visible_y_down)
            min_visible_y = min(min_visible_y, min_y)
            current_max_y = max(current_max_y, abs(min_y))
            if abs(min_y) > 1e-6: has_visible_data = True

        # Store for later threshold check and plotting
        # We plot a dummy line for the legend
        line, = ax.plot(dos.energies, y_up, lw=1.5, color=c, label=label, alpha=0)
        lines.append({
            'line': line,
            'y_up': y_up,
            'y_down': y_down,
            'max_y': current_max_y,
            'color': c,
            'label': label
        })
    
    # Calculate threshold (legend_cutoff of max visible)
    # Use the overall max absolute value found
    global_max = max(max_visible_y, abs(min_visible_y))
    threshold = legend_cutoff * global_max
    
    # Auto-scale Y-axis calculation (to determine ymax_ref)
    if ylim:
        pass # ymax_ref = ylim[1] (Unused)
    else:
        # Determine likely ymax based on logic later in the function
        pass

    # Filter legend items but keep all plot lines
    final_lines = []
    final_labels = []
    
    for item in lines:
        line = item['line']
        y_up = item['y_up']
        y_down = item['y_down']
        c = item['color']
        label = item['label']
        max_y = item['max_y']
        
        # Always plot the line (make it visible)
        line.set_alpha(1.0)
        
        # Apply gradient fill for visible lines
        # Spin Up
        gradient_fill(dos.energies, y_up, ax=ax, color=c, alpha=0.9)
        
        # Spin Down
        if is_spin_polarized:
             gradient_fill(dos.energies, y_down, ax=ax, color=c, alpha=0.9)
        
        # Only add to legend if above main threshold
        if max_y >= threshold:
            final_lines.append(line)
            final_labels.append(label)
    
    # Update lines and labels for legend creation
    lines = final_lines
    labels = final_labels

    # Plot Total DOS
    if show_total:
        y_total_up = dos.spin_up
        y_total_down = -dos.spin_down
        
        ax.plot(dos.energies, y_total_up, color="k", lw=1.2, label="Total DOS")
        # Use ymax_ref here too, assuming we want scaling relative to frame too
        # But Total DOS is often much larger. Maybe keep it separate max?
        # User said "Maximum value after scaling". So consistent reference is good.
        gradient_fill(dos.energies, y_total_up, ax=ax, color="k", alpha=0.15)
        
        if is_spin_polarized:
            ax.plot(dos.energies, y_total_down, color="k", lw=1.2)
            gradient_fill(dos.energies, y_total_down, ax=ax, color="k", alpha=0.15)
            
            # Update max/min range for auto-scaling
            visible_total_up = y_total_up[x_mask]
            visible_total_down = y_total_down[x_mask]
            if len(visible_total_up) > 0:
                max_visible_y = max(max_visible_y, np.max(visible_total_up))
            if len(visible_total_down) > 0:
                min_visible_y = min(min_visible_y, np.min(visible_total_down))
    
    # Auto-scale Y-axis based on visible range if ylim not provided
    if not ylim:
        if max_visible_y > 0 or min_visible_y < 0:
            # Apply scaling factor to the limit (zoom in)
            upper_limit = (max_visible_y * 1.1) / scale_factor
            lower_limit = (min_visible_y * 1.1) / scale_factor if is_spin_polarized else 0
            
            # If spin polarized, maybe make symmetric if user wants? 
            # For now, let's just use the data range.
            # But often symmetric is nicer. Let's stick to data range for now.
            
            ax.set_ylim(lower_limit, upper_limit)
    else:
        ax.set_ylim(*ylim)

    # Axis settings
    ax.set_xlim(*xlim)
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
        # Simplified overlap check for now
        loc, ncol = "upper right", 1
        
        # If spin polarized, upper right might cover spin up data
        # But usually it's fine.

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
