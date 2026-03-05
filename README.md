# Species Identifier

A graph-based suite for identifying and counting chemical species (adsorbates and gas-phase molecules) from molecular dynamics trajectories. This tool is specifically designed for surface science applications involving complex oxides such as perovskites (e.g., SrTiO<sub>3</sub>).

---

## Features

- **Surface-Inclusive Analysis**: Automatically distinguishes between gas-phase and chemisorbed species using Z-height mapping and bond connectivity.
- **Generalized Logic**: Easily customizable for diverse chemical systems by defining element-specific parameters in a YAML configuration.
- **Bond Hysteresis**: Maintains stable species identification by applying a distance margin to existing bonds between frames to account for vibrational fluctuations.

---

## Installation

```bash
git clone https://github.com/bl981219/species_identifier.git
cd species_identifier
pip install .
```

---

## Quickstart Example

A sample trajectory and configuration file are provided in the `examples/` directory.

```bash
cd examples
species-analyze --config config.yaml --xdatcar XDATCAR --out example_results.csv
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
| adsorbates | List | Elements tracked as molecules (e.g., `["C","H","O"]`). Order determines `species_map` keys. |
| lattice_non_metals | List | Elements forming the static slab surface (e.g., `["O"]`) |
| cutoffs | Dict | Bond distance thresholds (Å). Keys are element pairs (e.g., `C-H: 1.2` or `H-C: 1.2`). |
| species_map | Dict | Mapping from atom counts to species names. Keys correspond to `adsorbates` order. |

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

## Output Format

The output CSV contains species counts for each frame.

Columns represent chemical species identified during the simulation.

### Example

```
Frame,CH4,CO2,H2O,...
0,2,1,0,...
1,1,2,1,...
```