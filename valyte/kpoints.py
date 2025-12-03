"""
Interactive K-Point Generation Module
=====================================

Handles interactive generation of KPOINTS files based on user input and structure.
"""

import os
import sys
import numpy as np
from pymatgen.core import Structure
from pymatgen.io.vasp.inputs import Kpoints

def generate_kpoints_interactive():
    """
    Interactively generates a KPOINTS file based on user input and POSCAR.
    """
    print("\n🔮 Valyte K-Point Generator")
    
    # Check for POSCAR
    poscar_path = "POSCAR"
    if not os.path.exists(poscar_path):
        # Try finding any POSCAR* file
        files = [f for f in os.listdir('.') if f.startswith('POSCAR')]
        if files:
            poscar_path = files[0]
            print(f"   Found structure: {poscar_path}")
        else:
            print("❌ POSCAR file not found in current directory.")
            return

    try:
        structure = Structure.from_file(poscar_path)
    except Exception as e:
        print(f"❌ Error reading structure: {e}")
        return

    # Scheme Selection
    print("\nSelect K-Mesh Scheme:")
    print("  1. Monkhorst-Pack")
    print("  2. Gamma (Default)")
    
    choice = input("   > ").strip()
    
    if choice == '1':
        scheme = 'MP'
    else:
        scheme = 'Gamma'

    # K-Spacing Input
    print("\nEnter K-Spacing (units of 2π/Å):")
    print("   (Typical values: 0.03 - 0.04)")
    
    try:
        kspacing_str = input("   > ").strip()
        if not kspacing_str:
            kspacing = 0.04 # Default
        else:
            kspacing = float(kspacing_str)
    except ValueError:
        print("❌ Invalid number. Exiting.")
        return

    if kspacing <= 0:
        print("ℹ️  Using Gamma-Only (1 1 1)")
        kpts = Kpoints.gamma_automatic((1, 1, 1))
        grid = (1, 1, 1)
    else:
        # Calculate grid based on spacing
        # Formula: N = |b| / (spacing * 2*pi)
        # pymatgen reciprocal lattice lengths include the 2*pi factor.
        
        recip_lattice = structure.lattice.reciprocal_lattice
        b_lengths = recip_lattice.abc
        
        # We multiply spacing by 2*pi because the input is a coefficient of 2*pi/A?
        # Or rather, standard convention (like VASP KSPACING) often implies 2*pi is involved in the density.
        # Empirically: N = |b| / (input * 2*pi) matches expected results.
        grid = [max(1, int(l / (kspacing * 2 * np.pi) + 0.5)) for l in b_lengths]
        
        # Create Kpoints object
        if scheme == 'MP':
            kpts = Kpoints.monkhorst_automatic(grid)
        else:
            kpts = Kpoints.gamma_automatic(grid)

    # Print Summary
    print("\n📊 Summary")
    print(f"   Structure: {structure.formula}")
    print(f"   Lattice:   a={structure.lattice.a:.2f}, b={structure.lattice.b:.2f}, c={structure.lattice.c:.2f} Å")
    print(f"   K-Mesh:    {grid[0]} x {grid[1]} x {grid[2]} ({scheme})")

    # Write KPOINTS
    output_file = "KPOINTS"
    kpts.write_file(output_file)
    print(f"\n✅ Generated {output_file}!")
