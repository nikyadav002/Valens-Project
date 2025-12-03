"""
Supercell generation module for Valyte.
"""

import os
from pymatgen.core import Structure


def create_supercell(poscar_path="POSCAR", nx=1, ny=1, nz=1, output="POSCAR_supercell"):
    """
    Creates a supercell from a POSCAR file.
    
    Args:
        poscar_path (str): Path to input POSCAR file.
        nx (int): Supercell size in x direction.
        ny (int): Supercell size in y direction.
        nz (int): Supercell size in z direction.
        output (str): Output filename for the supercell POSCAR.
    """
    
    if not os.path.exists(poscar_path):
        raise FileNotFoundError(f"{poscar_path} not found")
    
    # Read structure
    structure = Structure.from_file(poscar_path)
    
    # Create supercell
    supercell = structure.copy()
    supercell.make_supercell([nx, ny, nz])
    
    # Write output
    supercell.to(filename=output, fmt="poscar")
    
    supercell_atoms = len(supercell)
    print(f"✅ Supercell created: {output} ({supercell_atoms} atoms)")
