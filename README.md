# Species Identifier

A professional, graph-based suite for identifying and counting chemical species (adsorbates and gas-phase molecules) from molecular dynamics trajectories. This tool is specifically designed for surface science applications involving complex oxides such as perovskites (e.g., SrTiO₃, LaSrFeO₃).

---

## Features

- **Surface-Inclusive Analysis**: Automatically distinguishes between gas-phase, physisorbed, and chemisorbed species using Z-height mapping and bond connectivity.
- **Generalized Logic**: Easily customizable for diverse chemical systems by defining element-specific parameters in a YAML configuration.
- **Bond Hysteresis**: Maintains stable species identification by applying a distance margin to existing bonds between frames to account for vibrational fluctuations.

---

## Installation

Using a **Conda environment** is recommended for managing scientific dependencies like ASE and NumPy.

### 1. Clone the repository
```bash
git clone [https://github.com/bl981219/species_identifier.git](https://github.com/bl981219/species_identifier.git)
cd species_identifier

---

## Installation

It is recommended to use a virtual environment to manage dependencies.

### 1. Clone the repository

```bash
git clone https://github.com/bl981219/species_identifier.git
cd species_identifier
```

### 2. Create and Activate Conda Environment

```bash
conda create -n species_env python=3.9 -y
conda activate species_env
```

### 3. Install Dependencies and Package

Install the required libraries from your 'requirements.txt' and then install the package in standard mode (or editable mode `-e` for development).

```bash
pip install -r requirements.txt
pip install .
```

---

## Quickstart Example

A sample trajectory and configuration file are provided in the `examples/` directory.

```bash
cd examples

species-analyze \
  --config config.yaml \
  --xdatcar XDATCAR \
  --out example_results.csv
```

This command will:

- Process the sample trajectory
- Identify species in each frame
- Output `example_results.csv` with time-resolved species counts

---

## Configuration

The analysis is controlled by `config.yaml`. This allows customization for different chemical systems without modifying the source code.

### Key Parameters

| Parameter | Type | Description |
|---------|------|-------------|
| adsorbate_elements | List | Elements tracked as molecules (e.g., `["C","H","O"]`) |
| lattice_elements | List | Elements forming the static slab surface (e.g., `["O"]`) |
| cutoffs | Dict | Bond distance thresholds (Å). Keys must be alphabetically ordered (e.g., `C-H: 1.2`) |
| species | Dict | Mapping from atom counts to species names (e.g., `"1,4,0": "CH4"`) |

---

## Repository Structure

The project follows a professional `src/` layout to ensure a clean namespace and prevent collisions.

```
.
├── pyproject.toml
├── README.md
├── config.yaml
├── requirements.txt
├── src/
│   └── species_identifier/
│       ├── __init__.py
│       └── analyzer.py
└── examples/
    ├── XDATCAR
    ├── config.yaml
    └── example_results.csv
```

### Description

| File | Purpose |
|------|---------|
| pyproject.toml | Build system and CLI entry points |
| README.md | Documentation |
| config.yaml | Configuration template |
| requirements.txt | Dependency list |
| analyzer.py | Core analysis logic |
| examples/ | Example files |

---

## Usage

Run the tool on your own data:

```bash
species-analyze \
  --config your_config.yaml \
  --xdatcar your_XDATCAR \
  --out your_results.csv
```

---

## Requirements

Dependencies are automatically installed with:

```bash
pip install .
```

### Required Packages

- Python ≥ 3.9
- NumPy — Numerical array operations
- ASE (Atomic Simulation Environment) — Trajectory handling
- PyYAML — Configuration management

---

## Output Format

The output CSV contains species counts for each frame.

Columns represent chemical species identified during the simulation.

### Example

```
Frame,CH4,CO2,H2O
0,2,1,0
1,1,2,1
```