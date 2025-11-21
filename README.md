<p align="center">
  <img src="valens/Logo.png" alt="Valens Logo" width="400"/>
</p>

# Valens

**Valens** is a comprehensive CLI tool for VASP workflows, providing both pre-processing and post-processing capabilities with a focus on clean, publication-quality outputs and modern aesthetics.

## Features

### Pre-processing
- **Supercell Creation**: Generate supercells from POSCAR files.
- **Band KPOINTS Generation**: Automatic high-symmetry path detection for band structure calculations.

### Post-processing
- **Smart Plotting**: Automatically handles total DOS and Projected DOS (PDOS).
- **Orbital-Resolved**: Plots individual orbitals (s, p, d, f) by default.
- **Adaptive Legend**: Intelligently hides the legend if PDOS contributions are low.
- **Gradient Fill**: Aesthetically pleasing gradient fills for DOS peaks.
- **Custom Fonts**: Supports Arial, Helvetica, Times New Roman.

> [!NOTE]
> Band structure plotting is coming soon!

## Installation

Clone the repository and install in editable mode:

```bash
git clone https://github.com/nikyadav002/Valens-Project
cd Valens-Project
pip install -e .
```

## Updating Valens

To update to the latest version:

```bash
cd Valens-Project
git pull
pip install -e .
```

## Usage

The main command is `valens`.

### Create Supercell

Generate a supercell from a POSCAR file:

```bash
valens supercell nx ny nz [options]
```

**Example:**
```bash
# Create a 2×2×2 supercell
valens supercell 2 2 2

# Specify input and output files
valens supercell 3 3 1 -i POSCAR_primitive -o POSCAR_3x3x1
```

**Options:**
- `-i`, `--input`: Input POSCAR file (default: `POSCAR`).
- `-o`, `--output`: Output filename (default: `POSCAR_supercell`).

### Generate Band Structure KPOINTS

Automatically generate KPOINTS file with high-symmetry path:

```bash
valens band kpt-gen [options]
```

**Example:**
```bash
# Generate KPOINTS with default settings (40 points per segment)
valens band kpt-gen

# Specify number of points and custom filenames
valens band kpt-gen -n 60 -i POSCAR_relaxed -o KPOINTS_band
```

**Options:**
- `-i`, `--input`: Input POSCAR file (default: `POSCAR`).
- `-n`, `--npoints`: Points per segment (default: `40`).
- `-o`, `--output`: Output filename (default: `KPOINTS`).

### Plot DOS

```bash
valens dos [path/to/vasprun.xml] [options]
```

You can provide the path as a positional argument, use the `--vasprun` flag, or omit it to use the current directory.

**Examples:**
```bash
# Plot all orbitals for all elements (Default)
valens dos

# Plot specific elements (Total PDOS)
valens dos -e Fe O

# Plot specific orbitals
valens dos -e "Fe(d)" "O(p)"

# Plot mixed (Fe Total and Fe d-orbital)
valens dos -e Fe "Fe(d)"
```

**Options:**
- `-e`, `--elements`: Specific elements or orbitals to plot.
    - Example: `-e Fe O` (Plots Total PDOS for Fe and O).
    - Example: `-e Fe(d) O(p)` (Plots Fe d-orbital and O p-orbital).
    - Example: `-e Fe Fe(d)` (Plots Fe Total and Fe d-orbital).
- `--xlim`: Energy range (default: `-6 6`).
- `--ylim`: DOS range (e.g., `--ylim 0 10`).
- `--scale`: Scaling factor for Y-axis (e.g., `--scale 3` divides DOS by 3).
- `--fermi`: Draw a dashed line at the Fermi level (E=0). Default is OFF.
- `--pdos`: Plot only Projected DOS (hide Total DOS).
- `--legend-cutoff`: Threshold for legend visibility (default: `0.10` = 10%).
- `-o`, `--output`: Output filename (default: `valens_dos.png`).
- `--font`: Font family (default: `Arial`).

**Example:**

```bash
valens dos ./vasp_data --xlim -5 5 -o my_dos.png
```
