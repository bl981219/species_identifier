#!/usr/bin/env python3
import csv
import sys
import yaml
import argparse
import numpy as np
from collections import Counter
from pathlib import Path
from ase.io import read
from ase.neighborlist import neighbor_list
import logging

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(message)s")

def get_bond_cutoff(s1, s2, metals, cutoffs):
    pair = frozenset([s1, s2])
    if pair in cutoffs: return cutoffs[pair]
    k1 = 'M' if s1 in metals else s1
    k2 = 'M' if s2 in metals else s2
    pair_m = frozenset([k1, k2])
    return cutoffs.get(pair_m, 0.0)

def get_species_name(mol_syms, adsorbates, species_map):
    c = Counter(mol_syms)
    key = ",".join(str(c.get(el, 0)) for el in adsorbates)
    if key in species_map: return species_map[key]
    formula = "".join(f"{el}{c[el] if c[el] > 1 else ''}" for el in adsorbates if c[el] > 0)
    return formula if formula else "Unknown"

def build_graph(atoms, sub_syms, reactive_indices, prev_bonds, metals, is_lat_nm, cutoffs, config):
    max_cut = max(cutoffs.values()) + config['bond_hysteresis'] if cutoffs else 0.0
    i_list, j_list, d_list = neighbor_list('ijd', atoms, max_cut)
    adj = {int(g): set() for g in reactive_indices}
    new_bonds = set()
    g_to_l = {g: l for l, g in enumerate(reactive_indices)}
    r_set = set(reactive_indices)
    
    for i, j, d in zip(i_list, j_list, d_list):
        if i >= j or i not in r_set or j not in r_set: continue
        si, sj = sub_syms[g_to_l[i]], sub_syms[g_to_l[j]]
        if config.get('break_surface_bonds', True):
            if (si in metals and is_lat_nm.get(j)) or (sj in metals and is_lat_nm.get(i)): continue
        r_on = get_bond_cutoff(si, sj, metals, cutoffs)
        if r_on <= 0: continue
        pair = (i, j) if i < j else (j, i)
        if d < (r_on + (config['bond_hysteresis'] if pair in prev_bonds else 0.0)):
            adj[i].add(j); adj[j].add(i); new_bonds.add(pair)
    return adj, new_bonds

def get_components(adj, nodes):
    visited, components = set(), []
    for node in nodes:
        if node in visited: continue
        comp, stack = [], [node]
        visited.add(node);
        while stack:
            curr = stack.pop(); comp.append(curr)
            for n in adj.get(curr, []):
                if n not in visited: visited.add(n); stack.append(n)
        components.append(comp)
    return components

def main():
    parser = argparse.ArgumentParser(description="Graph-Based Species Identification Tool")
    parser.add_argument("--config", required=True, help="Path to config.yaml")
    parser.add_argument("--xdatcar", default="XDATCAR", help="MD Trajectory file")
    parser.add_argument("--out", default="species_results.csv", help="Output CSV name")
    args = parser.parse_args()

    with open(args.config, 'r') as f: config = yaml.safe_load(f)
    config['parsed_cutoffs'] = {frozenset(k.split('-')): float(v) for k, v in config.get('cutoffs', {}).items()}

    logging.info(f"Loading trajectory from {args.xdatcar}...")
    traj = read(args.xdatcar, index=':')
    
    # Initialize Surface from Frame 0
    f0 = traj[0]
    f0.pbc = list(config.get('pbc', [True, True, False]))
    syms, pos = np.array(f0.get_chemical_symbols()), f0.get_positions(wrap=True)
    adsorbates, lattice_nms = config.get('adsorbates', []), config.get('lattice_non_metals', [])
    metals = set(syms) - set(adsorbates) - set(lattice_nms)
    m_z = pos[[s in metals for s in syms], 2]
    surface_z = float(np.percentile(m_z, config.get('surface_metal_z_percentile', 98.0)))
    
    lat_nm_indices = {i for i, (s, z) in enumerate(zip(syms, pos[:, 2])) if s in lattice_nms and z < (surface_z + 0.7)}
    is_lat_nm = {i: (i in lat_nm_indices) for i in range(len(syms))}
    depths = config.get('depth_cutoffs', {})
    reactive = [i for i, (s, z) in enumerate(zip(syms, pos[:, 2])) if 
                (s in metals and z > surface_z + depths.get('surface_metals', -2.0)) or
                (i in lat_nm_indices and z > surface_z + depths.get('lattice_non_metals', -2.0)) or
                (s in adsorbates and z > surface_z + depths.get('adsorbates', -0.5))]
    reactive = np.array(reactive)

    logging.info(f"Surface established at Z={surface_z:.2f} Å. Analyzing {len(traj)} frames...")
    prev_bonds, results, all_keys = set(), [], set()
    
    for i, atoms in enumerate(traj):
        frame_counts = Counter()
        cur_syms = np.array(atoms.get_chemical_symbols())
        adj, prev_bonds = build_graph(atoms, cur_syms[reactive], reactive, prev_bonds, metals, is_lat_nm, config['parsed_cutoffs'], config)
        
        for mol in get_components(adj, reactive):
            is_ads = any(cur_syms[j] in metals or j in lat_nm_indices for j in mol)
            # Subtract metals and lattice non-metals to find species
            fragment = [j for j in mol if cur_syms[j] not in metals and j not in lat_nm_indices]
            if not fragment: continue
            
            name = get_species_name(cur_syms[fragment], adsorbates, config.get('species_map', {}))
            label = f"{name}*" if is_ads else f"{name}_gas"
            frame_counts[label] += 1
            all_keys.add(label)
        
        results.append({"Frame": i + 1, **frame_counts})
        if (i+1) % 10 == 0: sys.stdout.write(f"\rProgress: {i+1}/{len(traj)} frames")

    with open(args.out, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=["Frame"] + sorted(all_keys))
        writer.writeheader()
        for r in results: writer.writerow(r)
    logging.info(f"\n[Done] Results saved to {args.out}")

if __name__ == "__main__": main()