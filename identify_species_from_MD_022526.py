#!/usr/bin/env python3
"""
Strict Graph-Based Species Counting (Surface-Inclusive)
-------------------------------------------------------
METHODOLOGY:
1. Define Reactive Atoms (STATIC): C, H, plus Surface Metals (M) and Surface Oxygens (O) from Frame 0.
2. Build Adjacency Graph:
   - Connect atoms based on cutoffs.
   - CONSTRAINT: Do not bond Metals to Lattice Oxygens (keeps slab disjoint).
3. Partition Graph into Connected Components (Clusters).
4. Identify Species by Subtraction:
   - Cluster with Metal? -> Adsorbed. Remove Metal -> Identify Remainder.
   - Cluster with Lattice O? -> Adsorbed. Remove Lattice O -> Identify Remainder.
   - Floating Cluster? -> Check Z-height -> Identify.
"""

import csv
import sys
import argparse
import numpy as np
from collections import Counter
from pathlib import Path
from ase.io import read
from ase import Atoms

# ==========================================
# 1) PARAMETERS
# ==========================================

# BASE_CUTOFFS = {
#     ("H", "H"): 1.00,
#     ("C", "H"): 1.60,
#     ("C", "C"): 1.90,
#     ("H", "O"): 1.55,  # FIXED: 'H' before 'O' (and bumped to 1.55 for vibration)
#     ("O", "O"): 1.75,
#     ("C", "O"): 1.85,
#     ("C", "M"): 2.20,  # FIXED: 'C' before 'M'
#     ("M", "O"): 2.20,
#     ("H", "M"): 2.10,  # FIXED: 'H' before 'M'
# }
BASE_CUTOFFS = {
    ("H", "H"): 0.84,
    ("C", "H"): 1.19,
    ("C", "C"): 1.64,
    ("H", "O"): 1.08,  # FIXED: 'H' before 'O' (and bumped to 1.55 for vibration)
    ("O", "O"): 1.58,
    ("C", "O"): 1.53,
    ("C", "M"): 2.20,  # FIXED: 'C' before 'M'
    ("M", "O"): 2.10,
    ("H", "M"): 1.80,  # FIXED: 'H' before 'M'
}
# Bond hysteresis margin in Å:
# If a bond existed in the previous frame, keep it until d < (r_on + BOND_HYSTERESIS)
BOND_HYSTERESIS = 0.15  # typical 0.10–0.20 Å


PBC_ANALYSIS = (True, True, False)
SURFACE_METAL_Z_PERCENTILE = 98.0
# ADSORPTION_Z_THRESHOLD = 1.5
LATTICE_O_TOLERANCE = 0.7 

SPECIES_MAP = {
    (1, 4, 0): "CH4",   (1, 3, 0): "CH3",   (0, 0, 2): "O2",
    (0, 2, 1): "H2O",   (1, 0, 1): "CO",    (1, 0, 2): "CO2",
    (0, 2, 0): "H2",    (0, 1, 1): "OH",    (1, 2, 1): "CH2O",
    (1, 3, 1): "CH3O",  (1, 1, 1): "HCO",   (0, 0, 1): "O",
    (0, 1, 0): "H",
    (2, 4, 0): "C2H4",  (2, 6, 0): "C2H6",
    (2, 4, 1): "C2H4O", (2, 3, 1): "C2H3O",  (1, 1, 0): "CH",
    (1, 2, 0): "CH2"
}

OTHER_LABEL = "Other"

SPECIES_NAMES = list(dict.fromkeys(SPECIES_MAP.values())) + [OTHER_LABEL]
COLS_ADS = [f"{v}*" for v in SPECIES_NAMES]
COLS_GAS = [f"{v}_gas" for v in SPECIES_NAMES]
ALL_COLS = ["Frame"] + COLS_ADS + COLS_GAS


# ==========================================
# 2) HELPERS
# ==========================================

def get_bond_cutoff(s1: str, s2: str, metals: set) -> float:
    k1 = 'M' if s1 in metals else s1
    k2 = 'M' if s2 in metals else s2
    return float(BASE_CUTOFFS.get(tuple(sorted((k1, k2))), 0.0))

def count_key_from_syms(mol_syms: list) -> tuple[int, int, int]:
    c = Counter(mol_syms)
    return (int(c.get('C', 0)), int(c.get('H', 0)), int(c.get('O', 0)))

def build_adjacency_and_bonds(
    sub_syms: np.ndarray,
    reactive_indices: np.ndarray,
    dists: np.ndarray,
    prev_bonds: set[tuple[int, int]],
    metals: set,
    is_lattice_o_map: dict
):
    adj = {int(g): set() for g in reactive_indices}
    new_bonds: set[tuple[int, int]] = set()
    n_sub = len(reactive_indices)

    for i in range(n_sub):
        gi = int(reactive_indices[i]) # g is indices
        si = sub_syms[i] # s is element symbol
        
        i_is_metal = si in metals
        i_is_latO = is_lattice_o_map.get(gi, False)

        for j in range(i + 1, n_sub):
            gj = int(reactive_indices[j])
            sj = sub_syms[j]
            
            j_is_metal = sj in metals
            j_is_latO = is_lattice_o_map.get(gj, False)

            if (i_is_metal and j_is_latO) or (j_is_metal and i_is_latO): # checks if the pair of atoms consists of one Metal and one Lattice Oxygen.
                continue

            r_on = get_bond_cutoff(si, sj, metals)
            if r_on <= 0.0:
                # print(f'This bond is not defined somehow: {si}-{sj}')
                continue

            pair = (gi, gj) if gi < gj else (gj, gi)
            threshold = r_on + (BOND_HYSTERESIS if pair in prev_bonds else 0.0)

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
# 3) FRAME ANALYSIS (Optimized)
# ==========================================

def analyze_frame(atoms: Atoms, metals: set, prev_bonds: set[tuple[int, int]], 
                  surface_z: float, lattice_o_indices: set, 
                  is_lattice_o_map: dict, reactive_indices: np.ndarray,frame_idx: int):
    """
    Analyzes a single frame using pre-calculated static surface parameters.
    """
    atoms.pbc = list(PBC_ANALYSIS)
    syms = np.array(atoms.get_chemical_symbols())
    pos = atoms.get_positions(wrap=True)

    ads_counts = Counter()
    gas_counts = Counter()

    if len(reactive_indices) == 0:
        return ads_counts, gas_counts, set()

    # --- 1. Build Distance Matrix (MIC) ONLY for reactive atoms ---
    sub_pos = pos[reactive_indices]
    sub_syms = syms[reactive_indices]
    sub_atoms = Atoms(symbols=sub_syms, positions=sub_pos, cell=atoms.cell, pbc=atoms.pbc)
    dists = sub_atoms.get_all_distances(mic=True)

    # --- 2. Build Graph ---
    adj, new_bonds = build_adjacency_and_bonds(
        sub_syms=sub_syms,
        reactive_indices=reactive_indices,
        dists=dists,
        prev_bonds=prev_bonds,
        metals=metals,
        is_lattice_o_map=is_lattice_o_map
    )

    # --- 3. Find Molecules (Connected Components) ---
    molecules = get_connected_components(adj, reactive_indices)

    # --- 4. Classify Molecules (Subtraction Logic) ---
    for mol_idxs in molecules:
        m_atoms = [i for i in mol_idxs if syms[i] in metals]
        lat_o_atoms = [i for i in mol_idxs if i in lattice_o_indices]
        
        is_adsorbed = False
        species_indices = []
        
        # Rule A & B Combined: Chemisorbed on Metal AND/OR Lattice Oxygen
        if len(m_atoms) > 0 or len(lat_o_atoms) > 0:
            is_adsorbed = True
            # Subtraction: Keep atoms that are NOT metals AND NOT lattice oxygens
            species_indices = [i for i in mol_idxs if (syms[i] not in metals) and (i not in lattice_o_indices)]
            
        # # Rule C: Floating / Gas Phase (Z-height check)
        # else:
        #     species_indices = mol_idxs
        #     min_z = np.min(pos[mol_idxs, 2])
            
        #     # If the lowest atom is near the surface, it is physisorbed/trapped
        #     if (min_z - surface_z) < ADSORPTION_Z_THRESHOLD:
        #         is_adsorbed = True 
        #     else:
        #         is_adsorbed = False # True gas phase

        else:
            # If no surface atoms are in the cluster, it MUST be gas phase
            is_adsorbed = False 
            species_indices = mol_idxs
                
        # Safety catch for completely empty clusters (e.g., pure metal/oxide clusters)
        if not species_indices:
            continue

        # Finally, identify and log the molecule
        final_syms = syms[species_indices]
        key = count_key_from_syms(final_syms) 
        name = SPECIES_MAP.get(key, OTHER_LABEL)
        
        if is_adsorbed:
            ads_counts[f"{name}*"] += 1
        else:
            gas_counts[f"{name}_gas"] += 1
            if name == 'H': # --- DEBUGGING SECTION FOR UNCOORDINATED HYDROGEN ---
                h_indices_str = ", ".join(map(str, species_indices))
                print(f"Frame {frame_idx} | H_gas detected at indices: {h_indices_str}")
                
                # Get the global index of this single H atom
                h_idx = species_indices[0] 
                
                all_dists = atoms.get_distances(h_idx, range(len(atoms)), mic=True)
                
                # Ignore distance to itself (which is 0.0)
                all_dists[h_idx] = np.inf 
                
                # Find the absolute closest atom in the whole box
                closest_idx = int(np.argmin(all_dists))
                closest_dist = all_dists[closest_idx]
                closest_elem = syms[closest_idx]
                
                # Check if this closest atom was actually in our graph
                in_graph = "YES" if closest_idx in reactive_indices else "NO (MASKED OUT!)"
                
                print(f"   -> Closest atom: Index {closest_idx:05d}, Element {closest_elem}, Distance {closest_dist:.2f} Å")
                print(f"   -> Was closest atom included in the reactive graph? {in_graph}")

                # Log spatial coordinates
                p = pos[h_idx]
                print(f"   -> Pos: {p[0]:8.3f}, {p[1]:8.3f}, {p[2]:8.3f}\n")

    return ads_counts, gas_counts, new_bonds


# ==========================================
# 4) MAIN
# ==========================================

def main():
    parser = argparse.ArgumentParser(description="Strict Graph-Based Species Counting (Surface-Inclusive).")
    parser.add_argument("--xdatcar", default="XDATCAR", help="Path to XDATCAR")
    parser.add_argument("--poscar_init", default="POSCAR_initial", help="Path to POSCAR_initial")
    parser.add_argument("--out", default="species_conservation.csv", help="Output CSV name")
    args = parser.parse_args()

    try:
        script_dir = Path(__file__).resolve().parent
    except NameError:
        script_dir = Path.cwd()

    f_xdatcar = (script_dir / args.xdatcar).resolve()
    f_poscar_init = (script_dir / args.poscar_init).resolve()
    f_output = (script_dir / args.out).resolve()

    if not f_poscar_init.exists():
        sys.exit(f"Error: {f_poscar_init} not found.")
    if not f_xdatcar.exists():
        sys.exit(f"Error: {f_xdatcar} not found.")

    # --- Step A: Read Initial State & Identify Metals ---
    print(f"Reading {f_poscar_init}...")
    init_atoms = read(str(f_poscar_init))
    init_syms = np.array(init_atoms.get_chemical_symbols())

    metals_list = sorted({s for s in init_syms if s not in ('C', 'H', 'O')})
    metals_set = set(metals_list)
    if not metals_set:
        sys.exit("Error: No metals detected in POSCAR_initial.")
    print(f"Detected metals: {metals_list}")

    # --- Step B: Pre-calculate Static Surface Properties from First Frame ---
    print(f"Reading {f_xdatcar}...")
    
    # Note: Changed to ':' to process the entire trajectory instead of just ':100'
    traj = read(str(f_xdatcar), index=':') 
    n_frames = len(traj)
    
    # Extract properties from the very first frame of the trajectory
    frame_zero = traj[0]
    frame_zero.pbc = list(PBC_ANALYSIS)
    f0_syms = np.array(frame_zero.get_chemical_symbols())
    f0_pos = frame_zero.get_positions(wrap=True)

    is_metal_arr = np.array([s in metals_set for s in f0_syms])
    surface_z = float(np.percentile(f0_pos[is_metal_arr, 2], SURFACE_METAL_Z_PERCENTILE))
    print(f"Calculated surface Z-height from frame 0: {surface_z:.2f} Å")

    is_O_arr = (f0_syms == 'O')
    lattice_o_indices = set(np.where(is_O_arr & (f0_pos[:, 2] < (surface_z + LATTICE_O_TOLERANCE)))[0])
    print(f"Identified {len(lattice_o_indices)} lattice O atoms in frame 0.")
    is_lattice_o_map = {idx: (idx in lattice_o_indices) for idx in range(len(f0_syms))}

    z_cutoff = surface_z - 0.5
    is_C = (f0_syms == 'C') & (f0_pos[:, 2] > z_cutoff)
    is_H = (f0_syms == 'H') & (f0_pos[:, 2] > z_cutoff)
    is_valid_M = is_metal_arr & (f0_pos[:, 2] > z_cutoff - 2) # include subsurface metal
    is_valid_O = is_O_arr & (f0_pos[:, 2] > surface_z - 2) # include subsurface oxygen into analysis
    
    reactive_mask = is_C | is_H | is_valid_M | is_valid_O
    reactive_indices = np.where(reactive_mask)[0].astype(int)

    # anaylze the reactive atoms in frame 0
    print('the highest Z of the reactive oxygens in frame 0 is: {:.2f} Å'.format(np.max(f0_pos[is_valid_O, 2]) if np.any(is_valid_O) else 0.0))
    print('the lowest Z of the reactive oxygens in frame 0 is: {:.2f} Å'.format(np.min(f0_pos[is_valid_O, 2]) if np.any(is_valid_O) else 0.0))
    print(f"--- Static Surface Model Initialized ---")
    print(f"Calculated Surface Z-Height: {surface_z:.2f} Å")
    print("Number of reactive atoms: C={}, H={}, M={}, O={}".format(
        np.sum(is_C), np.sum(is_H), np.sum(is_valid_M), np.sum(is_valid_O)
    ))

    # --- Step C: Process Trajectory ---
    prev_bonds = set()

    print(f"Analyzing {n_frames} frames...")
    with open(f_output, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(ALL_COLS)

        for i, atoms in enumerate(traj):
            # Pass all the static structural data into the analyzer
            ads, gas, new_bonds = analyze_frame(
                atoms=atoms, 
                metals=metals_set, 
                prev_bonds=prev_bonds,
                surface_z=surface_z,
                lattice_o_indices=lattice_o_indices,
                is_lattice_o_map=is_lattice_o_map,
                reactive_indices=reactive_indices,
                frame_idx=i + 1
            )
            prev_bonds = new_bonds

            row = [i + 1]
            row.extend([ads.get(k, 0) for k in COLS_ADS])
            row.extend([gas.get(k, 0) for k in COLS_GAS])
            writer.writerow(row)

            if (i + 1) % 1 == 0:
                print(f"Processed {i + 1}/{n_frames}...", end='\r')

    print(f"\n[Done] Wrote {f_output.name}.")

if __name__ == "__main__":
    main()