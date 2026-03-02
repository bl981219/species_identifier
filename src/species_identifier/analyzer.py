import argparse
import csv
import sys
import yaml
import logging
import numpy as np
from ase.io import read
from collections import defaultdict
from typing import Dict, List, Set, FrozenSet

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")

def parse_args():
    parser = argparse.ArgumentParser(description="Species Identifier for MD Trajectories.")
    parser.add_argument("--config", required=True, help="Path to config.yaml")
    parser.add_argument("--xdatcar", required=True, help="Path to XDATCAR or trajectory")
    parser.add_argument("--out", required=True, help="Output CSV file path")
    return parser.parse_args()

def load_config(path: str) -> dict:
    with open(path, "r") as f:
        return yaml.safe_load(f)

def get_pairwise_cutoffs(cutoffs_dict: Dict[str, float]) -> Dict[FrozenSet[str], float]:
    """Parse 'A-B': 1.2 into a dictionary keyed by frozenset to safely support homonuclear bonds."""
    parsed = {}
    for pair_str, dist in cutoffs_dict.items():
        elements = tuple(pair_str.split('-'))
        if len(elements) == 2:
            parsed[frozenset(elements)] = float(dist)
    return parsed

def build_graph(atoms, pairwise_cutoffs, hysteresis_margin, previous_bonds):
    """
    Build adjacency list for the molecule graph.
    Applies hysteresis margin to previously existing bonds.
    """
    symbols = atoms.get_chemical_symbols()
    n_atoms = len(atoms)
    
    # Use ASE's built-in robust PBC distance calculator
    distances = atoms.get_all_distances(mic=True)
    
    adj = {i:[] for i in range(n_atoms)}
    current_bonds = set()
    
    for i in range(n_atoms):
        for j in range(i + 1, n_atoms):
            pair = frozenset([symbols[i], symbols[j]])
            if pair in pairwise_cutoffs:
                cutoff = pairwise_cutoffs[pair]
                
                # Apply Bond Hysteresis to stabilize molecular identification against vibrations
                if frozenset([i, j]) in previous_bonds:
                    cutoff += hysteresis_margin
                
                if distances[i, j] <= cutoff:
                    adj[i].append(j)
                    adj[j].append(i)
                    current_bonds.add(frozenset([i, j]))
                    
    return adj, current_bonds

def find_molecules(adj: dict, valid_indices: set) -> List[List[int]]:
    """Find connected components representing molecules among valid adsorbate indices."""
    visited = set()
    molecules =[]
    
    for i in valid_indices:
        if i not in visited:
            comp =[]
            queue = [i]
            visited.add(i)
            while queue:
                curr = queue.pop(0)
                comp.append(curr)
                for neighbor in adj[curr]:
                    if neighbor in valid_indices and neighbor not in visited:
                        visited.add(neighbor)
                        queue.append(neighbor)
            molecules.append(comp)
    return molecules

def get_species_name(molecule_indices: List[int], symbols: List[str], adsorbate_elements: List[str], species_map: dict) -> str:
    """Determine the name of the species based on atom counts."""
    counts = {el: 0 for el in adsorbate_elements}
    for idx in molecule_indices:
        sym = symbols[idx]
        if sym in counts:
            counts[sym] += 1
            
    # Lookup by exactly formatted count mapping (e.g., "1,4,0")
    count_key = ",".join(str(counts[el]) for el in adsorbate_elements)
    if count_key in species_map:
        return species_map[count_key]
    
    # Fallback to smart empirical formula (e.g., "C1H4" -> "CH4")
    formula = "".join(f"{el}{counts[el] if counts[el] > 1 else ''}" for el in adsorbate_elements if counts[el] > 0)
    return formula if formula else "Unknown"

def analyze_trajectory(traj, config: dict):
    adsorbate_elements = config.get("adsorbate_elements", [])
    lattice_elements = config.get("lattice_elements",[])
    cutoffs = get_pairwise_cutoffs(config.get("cutoffs", {}))
    species_map = config.get("species", {})
    hysteresis_margin = config.get("hysteresis_margin", 0.15)
    dynamic_indices = config.get("dynamic_indices", False)
    
    z_threshold = config.get("surface_z_max", None)
    if z_threshold is None:
        logging.warning("No 'surface_z_max' specified in config. Defaulting to 15.0 Å.")
        z_threshold = 15.0
    
    results =[]
    previous_bonds = set()
    all_species_keys = set(species_map.values())
    surface_indices, adsorbate_indices = set(), set()
    
    for frame_idx, atoms in enumerate(traj):
        symbols = atoms.get_chemical_symbols()
        positions = atoms.get_positions()
        
        # 1. Identify surface vs adsorbate atoms 
        # By default, compute once for frame 0 to avoid atoms incorrectly swapping roles mid-simulation
        if dynamic_indices or frame_idx == 0:
            surface_indices.clear()
            adsorbate_indices.clear()
            for i, (sym, pos) in enumerate(zip(symbols, positions)):
                if sym in lattice_elements and pos[2] <= z_threshold:
                    surface_indices.add(i)
                elif sym in adsorbate_elements:
                    adsorbate_indices.add(i)
                    
        # 2. Build connectivity graph
        adj, current_bonds = build_graph(atoms, cutoffs, hysteresis_margin, previous_bonds)
        previous_bonds = current_bonds
        
        # 3. Find connected components strictly amongst adsorbates
        molecules = find_molecules(adj, adsorbate_indices)
        
        # 4. Classify each molecule as gas-phase or chemisorbed
        frame_counts = defaultdict(int)
        
        for mol in molecules:
            # Graph logic: A molecule is chemisorbed if any of its atoms neighbors a surface atom
            is_chemisorbed = any(neighbor in surface_indices for idx in mol for neighbor in adj[idx])
            species_name = get_species_name(mol, symbols, adsorbate_elements, species_map)
            
            # Prefix with '*' if chemisorbed to distinguish state
            if is_chemisorbed:
                species_name = f"*{species_name}"
            
            frame_counts[species_name] += 1
            all_species_keys.add(species_name)
            
        results.append({"Frame": frame_idx, **frame_counts})
        
    return results, sorted(list(all_species_keys))

def main():
    args = parse_args()
    try:
        config = load_config(args.config)
    except Exception as e:
        logging.error(f"Failed to load config: {e}")
        sys.exit(1)
        
    logging.info(f"Loading trajectory from {args.xdatcar}...")
    try:
        traj = read(args.xdatcar, index=":")
        logging.info(f"Successfully loaded {len(traj)} frames.")
    except Exception as e:
        logging.error(f"Error reading trajectory: {e}")
        sys.exit(1)
        
    logging.info("Starting species identification...")
    results, species_columns = analyze_trajectory(traj, config)
    
    logging.info(f"Writing time-resolved results to {args.out}...")
    with open(args.out, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["Frame"] + species_columns)
        writer.writeheader()
        for row in results:
            complete_row = {col: row.get(col, 0) for col in writer.fieldnames}
            complete_row["Frame"] = row["Frame"]
            writer.writerow(complete_row)
            
    logging.info("Analysis complete.")

if __name__ == "__main__":
    main()