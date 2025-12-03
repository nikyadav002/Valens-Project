#!/usr/bin/env python3
"""
Valyte CLI Tool
===============

A post-processing tool for VASP outputs, designed to create publication-quality
plots with a modern aesthetic. Supports DOS and band structure plotting.

Features:
- DOS plotting with gradient fills
- Band structure plotting
- Smart legend positioning
- Custom font support
"""

import os
import sys
import argparse
import re
import warnings

# Suppress pymatgen warnings
warnings.filterwarnings("ignore", category=UserWarning, module="pymatgen")

from valyte.supercell import create_supercell
from valyte.band import generate_band_kpoints
from valyte.band_plot import plot_band_structure
from valyte.dos_plot import load_dos, plot_dos
from valyte.kpoints import generate_kpoints_interactive

def parse_element_selection(inputs):
    """
    Parses user input for elements and orbitals.
    
    Args:
        inputs (list): List of strings, e.g., ["Ag", "Bi(s)", "O(p)"]
        
    Returns:
        tuple: (elements_to_load, plotting_config)
            elements_to_load (list): Elements to extract from VASP data.
            plotting_config (list): List of (Element, Orbital) tuples.
    """
    if not inputs:
        return None, None

    elements_to_load = set()
    plotting_config = []
    
    # Regex to match "Element" or "Element(orbital)"
    pattern = re.compile(r"^([A-Za-z]+)(?:\(([spdf])\))?$")
    
    for item in inputs:
        match = pattern.match(item)
        if match:
            el = match.group(1)
            orb = match.group(2) # None if no orbital specified
            
            elements_to_load.add(el)
            if orb:
                plotting_config.append((el, orb))
            else:
                plotting_config.append((el, 'total'))
        else:
            print(f"⚠️ Warning: Could not parse '{item}'. Ignoring.")
            

    return list(elements_to_load), plotting_config


# ===============================================================
# Main CLI
# ===============================================================
def main():
    """
    Main entry point for the CLI.
    """
    parser = argparse.ArgumentParser(description="Valyte: VASP Post-Processing Tool")
    subparsers = parser.add_subparsers(dest="command", help="Available commands")

    # --- DOS Subcommand ---
    dos_parser = subparsers.add_parser("dos", help="Plot Density of States (DOS)")
    dos_parser.add_argument("filepath", nargs="?", help="Path to vasprun.xml or directory containing it (optional)")
    dos_parser.add_argument("--vasprun", help="Explicit path to vasprun.xml (alternative to positional argument)")
    dos_parser.add_argument("-e", "--elements", nargs="+", help="Elements/Orbitals to plot (e.g., 'Fe O' or 'Fe(d) O(p)')")
    dos_parser.add_argument("-o", "--output", default="valyte_dos.png", help="Output filename")
    dos_parser.add_argument("--xlim", nargs=2, type=float, default=[-6, 6], help="Energy range (min max)")
    dos_parser.add_argument("--ylim", nargs=2, type=float, help="DOS range (min max)")
    dos_parser.add_argument("--scale", type=float, default=1.0, help="Scaling factor for Y-axis (zoom in)")
    dos_parser.add_argument("--fermi", action="store_true", help="Draw dashed line at Fermi level (E=0)")
    dos_parser.add_argument("--pdos", action="store_true", help="Plot only Projected DOS (hide Total DOS)")
    dos_parser.add_argument("--legend-cutoff", type=float, default=0.10, help="Threshold for legend visibility (0.0-1.0)")
    dos_parser.add_argument("--font", default="Arial", help="Font family")

    # --- Supercell Subcommand ---
    supercell_parser = subparsers.add_parser("supercell", help="Create a supercell")
    supercell_parser.add_argument("nx", type=int, help="Supercell size x")
    supercell_parser.add_argument("ny", type=int, help="Supercell size y")
    supercell_parser.add_argument("nz", type=int, help="Supercell size z")
    supercell_parser.add_argument("-i", "--input", default="POSCAR", help="Input POSCAR file")
    supercell_parser.add_argument("-o", "--output", default="POSCAR_supercell", help="Output filename")

    # --- Band Structure Subcommand ---
    band_parser = subparsers.add_parser("band", help="Band structure utilities")
    band_subparsers = band_parser.add_subparsers(dest="band_command", help="Band commands")

    # Band Plotting (default if no subcommand)
    # Note: argparse doesn't easily support default subcommands, so we handle this in logic
    band_parser.add_argument("--vasprun", default=".", help="Path to vasprun.xml or directory")
    band_parser.add_argument("--kpoints", help="Path to KPOINTS file (for labels)")
    band_parser.add_argument("-o", "--output", default="valyte_band.png", help="Output filename")
    band_parser.add_argument("--ylim", nargs=2, type=float, help="Energy range (min max)")
    band_parser.add_argument("--font", default="Arial", help="Font family")

    # Band KPOINTS Generation
    kpt_gen_parser = band_subparsers.add_parser("kpt-gen", help="Generate KPOINTS for band structure")
    kpt_gen_parser.add_argument("-i", "--input", default="POSCAR", help="Input POSCAR file")
    kpt_gen_parser.add_argument("-n", "--npoints", type=int, default=40, help="Points per segment")
    kpt_gen_parser.add_argument("-o", "--output", default="KPOINTS", help="Output filename")

    # --- K-Point Generation (Interactive) ---
    subparsers.add_parser("kpt", help="Interactive K-Point Generation (SCF)")

    args = parser.parse_args()

    if args.command == "dos":
        # Resolve filepath: positional > flag > current dir
        target_path = args.filepath if args.filepath else args.vasprun
        if not target_path:
            target_path = "."
            
        elements, plotting_config = parse_element_selection(args.elements)
        
        try:
            dos_data, pdos_data = load_dos(target_path, elements)
            plot_dos(
                dos_data, pdos_data, 
                out=args.output,
                xlim=tuple(args.xlim),
                ylim=tuple(args.ylim) if args.ylim else None,
                font=args.font,
                show_fermi=args.fermi,
                show_total=not args.pdos,
                plotting_config=plotting_config,
                legend_cutoff=args.legend_cutoff,
                scale_factor=args.scale
            )
            # print(f"✅ DOS plot saved to {args.output}") # Silent mode
        except Exception as e:
            print(f"❌ Error: {e}")
            sys.exit(1)

    elif args.command == "supercell":
        try:
            create_supercell(
                poscar_path=args.input,
                nx=args.nx,
                ny=args.ny,
                nz=args.nz,
                output=args.output
            )
        except Exception as e:
            print(f"❌ Error: {e}")
            sys.exit(1)
            
    elif args.command == "kpt":
        try:
            generate_kpoints_interactive()
        except Exception as e:
            print(f"❌ Error: {e}")
            sys.exit(1)

    elif args.command == "band":
        if args.band_command == "kpt-gen":
            try:
                generate_band_kpoints(
                    poscar_path=args.input,
                    npoints=args.npoints,
                    output=args.output
                )
            except Exception as e:
                print(f"❌ Error: {e}")
                sys.exit(1)
        elif args.band_command == "plot" or args.band_command is None:
             # Default behavior for 'valyte band' is plotting
            try:
                # Determine input path: --vasprun > positional > current dir
                target_path = args.vasprun or args.filepath or "."
                if os.path.isdir(target_path):
                    target_path = os.path.join(target_path, "vasprun.xml")
                
                # Determine KPOINTS path
                kpoints_path = args.kpoints
                if not kpoints_path:
                    # Try to find KPOINTS in the same directory as vasprun.xml
                    base_dir = os.path.dirname(target_path)
                    potential_kpoints = os.path.join(base_dir, "KPOINTS")
                    if os.path.exists(potential_kpoints):
                        kpoints_path = potential_kpoints

                plot_band_structure(
                    vasprun_path=target_path,
                    kpoints_path=kpoints_path,
                    output=args.output,
                    ylim=tuple(args.ylim) if args.ylim else None,
                    font=args.font
                )
            except Exception:
                import traceback
                traceback.print_exc()
                sys.exit(1)
        else:
            band_parser.print_help()
    else:
        parser.print_help()

if __name__ == "__main__":
    main()
