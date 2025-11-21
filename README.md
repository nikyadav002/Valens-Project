# Valens

**Valens** is a CLI-based post-processing tool for VASP outputs, designed to create publication-quality Density of States (DOS) plots with a clean, modern aesthetic.

## Features

- **Smart Plotting**: Automatically handles total DOS and Projected DOS (PDOS).
- **Adaptive Legend**: Intelligently hides the legend if PDOS contributions are low, keeping the plot clean.
- **Gradient Fill**: aesthetically pleasing gradient fills for DOS peaks.
- **Custom Fonts**: Supports Arial, Helvetica, Times New Roman.

> [!NOTE]
> Band structure plotting and other features are coming soon!

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
- `-o`, `--output`: Output filename (default: `valens_dos.png`).
- `--font`: Font family (default: `Arial`).

**Example:**

```bash
valens dos ./vasp_data --xlim -5 5 -o my_dos.png
```
