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

## Usage

The main command is `valens`.

### Plot DOS

```bash
valens dos [path/to/vasprun.xml] [options]
```

You can provide the path as a positional argument, use the `--vasprun` flag, or omit it to use the current directory.

**Examples:**
```bash
# Use current directory
valens dos

# Use positional argument
valens dos ./vasp_data

# Use flag
valens dos --vasprun ./vasp_data/vasprun.xml
```

**Options:**
- `-e`, `--elements`: Specific elements to plot PDOS for (e.g., `-e Fe O`).
- `--xlim`: Energy range (default: `-6 6`).
- `--ylim`: DOS range (e.g., `--ylim 0 10`).
- `--scale`: Scaling factor for Y-axis (e.g., `--scale 3` divides DOS by 3).
- `-o`, `--output`: Output filename (default: `valens_dos.png`).
- `--font`: Font family (default: `Arial`).

**Example:**

```bash
valens dos ./vasp_data --xlim -5 5 -o my_dos.png
```
