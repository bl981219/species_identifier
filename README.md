# Species Analyzer

A strict, graph-based tool to identify and count chemical species (adsorbates and gas-phase molecules) from molecular dynamics trajectories (e.g., VASP `XDATCAR`).

## Features
* **Surface-Inclusive:** Distinguishes between gas-phase, physisorbed, and chemisorbed species by mapping atom Z-heights and bond connectivity to a static surface model.
* **Customizable:** Uses a `config.yaml` file so you can easily swap bond cutoffs, tracking elements, and molecule definitions without editing raw Python code.

## Installation
1. Clone the repository.
2. Create a virtual environment and install dependencies:
   ```bash
   python -m venv venv
   source venv/bin/activate  # On Windows use: venv\Scripts\activate
   pip install -r requirements.txt

## Usage
Ensure you have your trajectory (`XDATCAR`) in the same directory, then run:
   ```bash
   python analyzer.py --config config.yaml --xdatcar XDATCAR --out results.csv

## Quickstart Example
Want to test the analyzer right away? We have provided a sample trajectory in the `examples/` directory.

Navigate to the examples folder and run the script:
```bash
cd examples
python ../analyzer.py --config config.yaml --xdatcar XDATCAR --out example_results.csv