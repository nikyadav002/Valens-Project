"""
Band structure KPOINTS generation module for Valens.
"""

import os
from pymatgen.core import Structure
from pymatgen.symmetry.bandstructure import HighSymmKpath


def generate_band_kpoints(poscar_path="POSCAR", npoints=40, output="KPOINTS"):
    """
    Generates KPOINTS file in line-mode for band structure calculations.
    Uses SeeK-path method for high-symmetry path determination.
    
    Args:
        poscar_path (str): Path to input POSCAR file.
        npoints (int): Number of points per segment (default: 40).
        output (str): Output filename for KPOINTS.
    """
    
    if not os.path.exists(poscar_path):
        raise FileNotFoundError(f"{poscar_path} not found")
    
    # Read structure
    structure = Structure.from_file(poscar_path)
    
    # Get high-symmetry path using SeeK-path method
    kpath = HighSymmKpath(structure, path_type="setyawan_curtarolo")
    
    # Get the path
    path = kpath.kpath["path"]
    kpoints = kpath.kpath["kpoints"]
    
    # Write KPOINTS file
    with open(output, 'w') as f:
        f.write("k-points for band structure\n")
        f.write(f"{npoints}\n")
        f.write("Line-mode\n")
        f.write("Reciprocal\n")
        
        # Write each segment
        for segment in path:
            for i in range(len(segment) - 1):
                start = segment[i]
                end = segment[i + 1]
                
                start_coords = kpoints[start]
                end_coords = kpoints[end]
                
                f.write(f"  {start_coords[0]:.6f}  {start_coords[1]:.6f}  {start_coords[2]:.6f}  ! {start}\n")
                f.write(f"  {end_coords[0]:.6f}  {end_coords[1]:.6f}  {end_coords[2]:.6f}  ! {end}\n")
                f.write("\n")
    
    # Print success message
    path_str = ' → '.join([' - '.join(seg) for seg in path])
    print(f"✅ KPOINTS generated: {output} ({path_str}, {npoints} pts/seg)")
