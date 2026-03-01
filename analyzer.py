#!/usr/bin/env python3
"""
Strict Graph-Based Species Counting (Surface-Inclusive)
"""

import csv
import sys
import argparse
import yaml
import numpy as np
from collections import Counter
from pathlib import Path
from ase.io import read
from ase import Atoms

# ==========================================
# 1) CONFIGURATION LOADER
# ==========================================
def load_config(config_path):
    with open(config_path, 'r') as f:
        config = yaml.safe_load(f)
    
    # Parse cutoffs into the tuple format your code expects: ("C", "H"): 1.19
    parsed_cutoffs = {}
    for bond, dist in config['cutoffs'].items():
        elem1, elem2 = bond.split('-')
        parsed_cutoffs[tuple(sorted((elem1, elem2)))] = float(dist)
        
    # Parse species map into tuples: (1, 4, 0): "CH4"
    parsed_species = {}
    for counts, name in config['species'].items():
        c_tuple = tuple(map(int, counts.split(',')))
        parsed_species[c_tuple] = name
        
    return config['parameters'], parsed_cutoffs, parsed_species

# ==========================================
# 2) HELPERS
# ==========================================
def get_bond_cutoff(s1: str, s2: str, metals: set, base_cutoffs: dict) -> float:
    k1 = 'M' if s1 in metals else s1
    k2 = 'M' if s2 in metals else s2
    return float(base_cutoffs.get(tuple(sorted((k1, k2))), 0.0))

def count_key_from_syms(mol_syms: list) -> tuple[int, int, int]:
    # This remains specific to C, H, O based on your logic, 
    # but the YAML map makes it easier to extend later.
    c = Counter(mol_syms)
    return (int(c.get('C', 0)), int(c.get('H', 0)), int(c.get('O', 0)))

def build_adjacency_and_bonds(sub_syms, reactive_indices, dists, prev_bonds, metals, is_lattice_o_map, base_cutoffs, hysteresis):
    adj = {int(g): set() for g in reactive_indices}
    new_bonds = set()
    n_sub = len(reactive_indices)

    for i in range(n_sub):
        gi = int(reactive_indices[i])
        si = sub_syms[i]
        i_is_metal = si in metals
        i_is_latO = is_lattice_o_map.get(gi, False)

        for j in range(i + 1, n_sub):
            gj = int(reactive_indices[j])
            sj = sub_syms[j]
            j_is_metal = sj in metals
            j_is_latO = is_lattice_o_map.get(gj, False)

            if (i_is_metal and j_is_latO) or (j_is_metal and i_is_latO):
                continue

            r_on = get_bond_cutoff(si, sj, metals, base_cutoffs)
            if r_on <= 0.0:
                continue

            pair = (gi, gj) if gi < gj else (gj, gi)
            threshold = r_on + (hysteresis if pair in prev_bonds else 0.0)

            if dists[i, j] < threshold:
                adj[gi].add(gj)
                adj[gj].add(gi)
                new_bonds.add(pair)

    return adj, new_bonds

def get_connected_components(adj: dict[int, set[int]], nodes: list[int]) -> list[list[int]]:
    visited = set()
    components = []
    for node in nodes:
        if node in visited:
            continue
        component = []
        stack = [node]
        visited.add(node)
        while stack:
            curr = stack.pop()
            component.append(curr)
            for neighbor in adj.get(curr, []):
                if neighbor not in visited:
                    visited.add(neighbor)
                    stack.append(neighbor)
        components.append(component)
    return components

# ==========================================
# 3) FRAME ANALYSIS
# ==========================================
def analyze_frame(atoms, metals, prev_bonds, surface_z, lattice_o_indices, 
                  is_lattice_o_map, reactive_indices, frame_idx, params, base_cutoffs, species_map):
    
    atoms.pbc = list(params['pbc_analysis'])
    syms = np.array(atoms.get_chemical_symbols())
    pos = atoms.get_positions(wrap=True)

    ads_counts = Counter()
    gas_counts = Counter()

    if len(reactive_indices) == 0:
        return ads_counts, gas_counts, set()

    sub_pos = pos[reactive_indices]
    sub_syms = syms[reactive_indices]
    sub_atoms = Atoms(symbols=sub_syms, positions=sub_pos, cell=atoms.cell, pbc=atoms.pbc)
    dists = sub_atoms.get_all_distances(mic=True)

    adj, new_bonds = build_adjacency_and_bonds(
        sub_syms, reactive_indices, dists, prev_bonds, metals, 
        is_lattice_o_map, base_cutoffs, params['bond_hysteresis']
    )

    molecules = get_connected_components(adj, reactive_indices)

    other_label = "Other"
    for mol_idxs in molecules:
        m_atoms = [i for i in mol_idxs if syms[i] in metals]
        lat_o_atoms = [i for i in mol_idxs if i in lattice_o_indices]
        
        is_adsorbed = False
        species_indices = []
        
        if len(m_atoms) > 0 or len(lat_o_atoms) > 0:
            is_adsorbed = True
            species_indices = [i for i in mol_idxs if (syms[i] not in metals) and (i not in lattice_o_indices)]
        else:
            is_adsorbed = False 
            species_indices = mol_idxs
                
        if not species_indices:
            continue

        final_syms = syms[species_indices]
        key = count_key_from_syms(final_syms) 
        name = species_map.get(key, other_label)
        
        if is_adsorbed:
            ads_counts[f"{name}*"] += 1
        else:
            gas_counts[f"{name}_gas"] += 1

    return ads_counts, gas_counts, new_bonds

# ==========================================
# 4) MAIN
# ==========================================
def main():
    parser = argparse.ArgumentParser(description="Strict Graph-Based Species Counting.")
    parser.add_argument("--config", default="config.yaml", help="Path to config.yaml")
    parser.add_argument("--xdatcar", default="XDATCAR", help="Path to XDATCAR")
    parser.add_argument("--out", default="species_conservation.csv", help="Output CSV name")
    args = parser.parse_args()

    # Load Configuration
    params, base_cutoffs, species_map = load_config(args.config)
    
    # Setup Columns dynamically based on config
    species_names = list(dict.fromkeys(species_map.values())) + ["Other"]
    cols_ads = [f"{v}*" for v in species_names]
    cols_gas = [f"{v}_gas" for v in species_names]
    all_cols = ["Frame"] + cols_ads + cols_gas

    print(f"Reading {args.xdatcar}...")
    traj = read(args.xdatcar, index=':') 
    n_frames = len(traj)
    
    # Extract properties directly from the first frame of the XDATCAR
    frame_zero = traj[0]
    frame_zero.pbc = list(params['pbc_analysis'])
    f0_syms = np.array(frame_zero.get_chemical_symbols())
    f0_pos = frame_zero.get_positions(wrap=True)

    # Dynamically find metals (anything not in the adsorbate list) from Frame 0
    ads_elements = set(params['adsorbate_elements'])
    metals_list = sorted({s for s in f0_syms if s not in ads_elements})
    metals_set = set(metals_list)
    if not metals_set:
        sys.exit("Error: No metals detected in frame 0.")
    print(f"Detected metals: {metals_list}")
# Calculate static surface properties
    is_metal_arr = np.array([s in metals_set for s in f0_syms])
    surface_z = float(np.percentile(f0_pos[is_metal_arr, 2], params['surface_metal_z_percentile']))
    
    # 1. Identify Lattice Surface Atoms (Generalized from Lattice O)
    lattice_elements = set(params['lattice_elements'])
    is_lattice_arr = np.isin(f0_syms, list(lattice_elements))
    
    # Find lattice atoms that act as the physical surface
    lattice_o_indices = set(np.where(is_lattice_arr & (f0_pos[:, 2] < (surface_z + params['lattice_z_tolerance'])))[0])
    is_lattice_o_map = {idx: (idx in lattice_o_indices) for idx in range(len(f0_syms))}

    # 2. Build the Reactive Mask Dynamically
    reactive_mask = np.zeros(len(f0_syms), dtype=bool)
    
    # Add Adsorbates (Dynamic elements, configurable depth)
    ads_elements = set(params['adsorbate_elements'])
    is_ads = np.isin(f0_syms, list(ads_elements))
    reactive_mask |= is_ads & (f0_pos[:, 2] > (surface_z - params['adsorbate_depth']))
    
    # Add Metals (Configurable depth)
    reactive_mask |= is_metal_arr & (f0_pos[:, 2] > (surface_z - params['metal_depth']))
    
    # Add Subsurface Lattice Atoms (Configurable depth)
    reactive_mask |= is_lattice_arr & (f0_pos[:, 2] > (surface_z - params['lattice_depth']))
    
    reactive_indices = np.where(reactive_mask)[0].astype(int)

    prev_bonds = set()
    print(f"Analyzing {n_frames} frames...")

    with open(args.out, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(all_cols)

        for i, atoms in enumerate(traj):
            ads, gas, new_bonds = analyze_frame(
                atoms, metals_set, prev_bonds, surface_z, lattice_o_indices,
                is_lattice_o_map, reactive_indices, i + 1, params, base_cutoffs, species_map
            )
            prev_bonds = new_bonds

            row = [i + 1]
            row.extend([ads.get(k, 0) for k in cols_ads])
            row.extend([gas.get(k, 0) for k in cols_gas])
            writer.writerow(row)
            print(f"Processed {i + 1}/{n_frames}...", end='\r')

    print(f"\n[Done] Wrote {args.out}.")

if __name__ == "__main__":
    main()