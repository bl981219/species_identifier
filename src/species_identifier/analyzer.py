#!/usr/bin/env python3
import csv
import sys
import os
import yaml
import argparse
import logging
from collections import Counter
from pathlib import Path
from typing import List, Dict, Set, Tuple, Any, Optional
import numpy as np
import multiprocessing as mp
from ase.io import read
from ase.atoms import Atoms
from ase.neighborlist import neighbor_list
from tqdm import tqdm

# Module-level logger
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())

def analyze_frame_standalone(
    atoms: Atoms, 
    frame_idx: int, 
    config_state: Dict[str, Any],
    prev_bonds: Optional[Set[Tuple[int, int]]] = None
) -> Tuple[int, Dict[str, int], List[Dict[str, Any]], Dict[int, str], Set[Tuple[int, int]]]:
    """Standalone function for parallel processing of a single frame."""
    atoms.pbc = config_state['pbc']
    syms = np.array(atoms.get_chemical_symbols())
    pos = atoms.get_positions(wrap=True)
    
    reactive_indices = config_state['reactive_indices']
    max_search_radius = config_state['max_search_radius']
    cutoff_matrix = config_state['cutoff_matrix']
    hysteresis = config_state['hysteresis']
    metals_set = config_state['metals_set']
    is_lattice_nm_map = config_state['is_lattice_nm_map']
    species_map = config_state['species_map']
    track_distances = config_state.get('track_distances', [])
    
    # 1. Component Atoms (Reactive subset)
    sub_syms = syms[reactive_indices]
    sub_pos = pos[reactive_indices]
    sub_atoms = Atoms(symbols=sub_syms, positions=sub_pos, cell=atoms.cell, pbc=atoms.pbc)
    
    # 2. Build Graph using neighbor_list
    ii, jj, dd = neighbor_list('ijd', sub_atoms, max_search_radius)
    
    adj = {int(g): set() for g in reactive_indices}
    new_bonds = set()
    
    for idx in range(len(ii)):
        i, j, dist = ii[idx], jj[idx], dd[idx]
        if i >= j: continue 
        
        gi, gj = int(reactive_indices[i]), int(reactive_indices[j])
        si, sj = sub_syms[i], sub_syms[j]
        
        # CONSTRAINT: Do not bond Metals to Lattice Non-Metals
        if (si in metals_set and is_lattice_nm_map.get(gj)) or (sj in metals_set and is_lattice_nm_map.get(gi)):
            continue
        
        pair = (gi, gj)
        r_limit = cutoff_matrix.get((si, sj), 0.0)
        if r_limit <= 0: continue
        
        threshold = r_limit + (hysteresis if (prev_bonds and pair in prev_bonds) else 0.0)
        
        if dist < threshold:
            adj[gi].add(gj)
            adj[gj].add(gi)
            new_bonds.add(pair)
    
    # 3. Find Connected Components
    components = []
    visited = set()
    for node in reactive_indices:
        if node not in visited:
            comp, stack = [], [int(node)]
            visited.add(int(node))
            while stack:
                curr = stack.pop()
                comp.append(curr)
                for n in adj.get(curr, []):
                    if n not in visited:
                        visited.add(n)
                        stack.append(n)
            components.append(comp)
    
    # 4. Classify by Subtraction
    frame_counts = Counter()
    bond_distances = []
    species_mapping = {}
    
    for mol_idxs in components:
        m_atoms = [i for i in mol_idxs if syms[i] in metals_set]
        lat_nm_atoms = [i for i in mol_idxs if is_lattice_nm_map.get(i)]
        
        # A component is "adsorbed" if it contains a metal or lattice non-metal
        is_adsorbed = (len(m_atoms) > 0 or len(lat_nm_atoms) > 0)
        
        # Species identification: REMOVE metals and REMOVE lattice non-metals
        species_indices = [i for i in mol_idxs if (syms[i] not in metals_set) and (not is_lattice_nm_map.get(i))]
        
        if not species_indices:
            continue
            
        final_syms = syms[species_indices]
        sc = Counter(final_syms)
        # Construct key using the adsorbates_list for consistency
        adsorbate_elements = config_state['adsorbates_list']
        key = ",".join(str(sc.get(el, 0)) for el in adsorbate_elements)
        
        name = species_map.get(key)
        if not name:
            # Fallback to basic chemical formula
            formula = "".join(f"{el}{sc[el] if sc[el] > 1 else ''}" for el in adsorbate_elements if sc[el] > 0)
            name = formula if formula else "Other"
        
        # Generic distance tracking for diatomics or specific requested species
        if name in track_distances and len(species_indices) == 2:
            dist = atoms.get_distance(species_indices[0], species_indices[1], mic=True)
            bond_distances.append({"Species": name, "Distance": dist})
        
        species_name = f"{name}{'*' if is_adsorbed else '_gas'}"
        frame_counts[species_name] += 1
        for idx in mol_idxs:
            species_mapping[int(idx)] = species_name
            
    return frame_idx + 1, dict(frame_counts), bond_distances, species_mapping, new_bonds

class TrajectoryAnalyzer:
    def __init__(self, config_path: str):
        self.config_path = Path(config_path)
        if not self.config_path.exists():
            raise FileNotFoundError(f"Config file not found: {config_path}")
            
        with open(self.config_path, 'r') as f:
            self.config = yaml.safe_load(f)
        
        self._validate_config()
        
        self.cutoffs_raw = self.config.get('cutoffs', {})
        self.hysteresis = self.config.get('bond_hysteresis', 0.15)
        self.adsorbates_list = self.config.get('adsorbates', ['C', 'H', 'O'])
        self.adsorbates_set = set(self.adsorbates_list)
        self.lattice_non_metals = self.config.get('lattice_non_metals', ['O'])
        self.species_map = self.config.get('species_map', {})
        self.pbc = self.config.get('pbc', [True, True, False])
        self.lattice_tol = self.config.get('lattice_non_metal_tolerance', 0.7)
        self.m_percentile = self.config.get('surface_metal_z_percentile', 98.0)
        self.track_distances = self.config.get('track_distances', ["O2"])
        
        self.depth_cutoffs = self.config.get('depth_cutoffs', {
            'adsorbates': -0.5,
            'surface_metals': -2.0,
            'lattice_non_metals': -2.0
        })
        
        self.metals_set = set()
        self.is_lattice_nm_map = {}
        self.reactive_indices = np.array([])
        self.surface_z = 0.0
        self.cutoff_matrix = {}
        self.max_search_radius = 0.0

    def _validate_config(self):
        """Basic validation of config keys and values."""
        required_keys = ['adsorbates', 'cutoffs', 'species_map']
        missing = [k for k in required_keys if k not in self.config]
        if missing:
            raise ValueError(f"Missing required config keys: {missing}")
        
        if not isinstance(self.config['cutoffs'], dict):
            raise ValueError("'cutoffs' must be a dictionary.")

    def _prepare_cutoff_matrix(self, elements: List[str]):
        """Build a fast lookup table for bond cutoffs between all element pairs."""
        self.cutoff_matrix = {}
        # Normalize cutoffs_raw keys to sorted order
        normalized_cutoffs = {}
        for k, v in self.cutoffs_raw.items():
            pair = tuple(sorted(k.split('-')))
            normalized_cutoffs[pair] = float(v)
            
        all_elems = set(elements) | self.metals_set | {'M'}
        max_cutoff = 0.0
        
        for e1 in all_elems:
            for e2 in all_elems:
                # 1. Try direct match
                pair = tuple(sorted([str(e1), str(e2)]))
                r_on = normalized_cutoffs.get(pair, 0.0)
                
                # 2. Try Metal generalization
                if r_on <= 0:
                    k1 = 'M' if e1 in self.metals_set else e1
                    k2 = 'M' if e2 in self.metals_set else e2
                    pair_m = tuple(sorted([str(k1), str(k2)]))
                    r_on = normalized_cutoffs.get(pair_m, 0.0)
                
                if r_on > 0:
                    self.cutoff_matrix[(e1, e2)] = float(r_on)
                    max_cutoff = max(max_cutoff, float(r_on))
        
        self.max_search_radius = max_cutoff + self.hysteresis

    def setup_static_surface(self, atoms: Atoms):
        """Setup static surface parameters from the initial frame."""
        syms = np.array(atoms.get_chemical_symbols())
        pos = atoms.get_positions(wrap=True)
        
        # 1. Identify Metals
        self.metals_list = sorted({s for s in syms if s not in self.adsorbates_set and s not in self.lattice_non_metals})
        self.metals_set = set(self.metals_list)
        logger.info(f"Detected metals: {self.metals_list}")
        
        # 2. Surface Z and Lattice Non-Metals
        is_metal_arr = np.array([s in self.metals_set for s in syms])
        if not any(is_metal_arr):
             self.surface_z = 0.0
             logger.warning("No metal atoms found; using Z=0.0 as surface reference.")
        else:
             self.surface_z = float(np.percentile(pos[is_metal_arr, 2], self.m_percentile))
        logger.info(f"Calculated surface Z-height: {self.surface_z:.2f} Å")
        
        is_lat_nm_arr = np.array([s in self.lattice_non_metals for s in syms])
        lat_nm_indices = set(np.where(is_lat_nm_arr & (pos[:, 2] < (self.surface_z + self.lattice_tol)))[0])
        self.is_lattice_nm_map = {idx: (idx in lat_nm_indices) for idx in range(len(syms))}
        logger.info(f"Identified {len(lat_nm_indices)} lattice non-metal atoms.")
        
        # 3. Define Reactive Indices
        # Use depth_cutoffs from config
        z_ads = self.surface_z + self.depth_cutoffs.get('adsorbates', -0.5)
        z_metal = self.surface_z + self.depth_cutoffs.get('surface_metals', -2.0)
        z_lat_nm = self.surface_z + self.depth_cutoffs.get('lattice_non_metals', -2.0)
        
        is_ads = np.array([s in self.adsorbates_set for s in syms]) & (pos[:, 2] > z_ads)
        is_valid_M = is_metal_arr & (pos[:, 2] > z_metal)
        is_valid_LNM = is_lat_nm_arr & (pos[:, 2] > z_lat_nm)
        
        reactive_mask = is_ads | is_valid_M | is_valid_LNM
        self.reactive_indices = np.where(reactive_mask)[0].astype(int)
        logger.info(f"Initialized {len(self.reactive_indices)} reactive atoms.")
        
        # 4. Prepare optimized cutoffs
        self._prepare_cutoff_matrix(list(set(syms)))

    def get_config_state(self) -> Dict[str, Any]:
        """Bundle serializable state for workers."""
        return {
            'pbc': self.pbc,
            'reactive_indices': self.reactive_indices,
            'max_search_radius': self.max_search_radius,
            'cutoff_matrix': self.cutoff_matrix,
            'hysteresis': self.hysteresis,
            'metals_set': self.metals_set,
            'is_lattice_nm_map': self.is_lattice_nm_map,
            'species_map': self.species_map,
            'adsorbates_set': self.adsorbates_set,
            'adsorbates_list': self.adsorbates_list,
            'track_distances': self.track_distances
        }

    def analyze(self, xdatcar_path: str, output_path: str, n_procs: int = 1, mapping_path: str = None, enable_hysteresis: bool = True):
        """Perform trajectory analysis and save results."""
        xdatcar = Path(xdatcar_path)
        if not xdatcar.is_file():
            raise FileNotFoundError(f"Trajectory file not found: {xdatcar_path}")
            
        logger.info(f"Reading trajectory from {xdatcar_path}...")
        traj = read(xdatcar_path, index=':')
        if not traj:
            raise ValueError(f"Trajectory file {xdatcar_path} is empty.")
            
        self.setup_static_surface(traj[0])
        config_state = self.get_config_state()
        
        results = []
        logger.info(f"Analyzing {len(traj)} frames using {n_procs} processes...")
        
        if n_procs > 1:
            logger.warning("Multiprocessing enabled: Bond hysteresis is DISABLED.")
            with mp.Pool(n_procs) as pool:
                args = [(atoms, i, config_state, None) for i, atoms in enumerate(traj)]
                results = list(tqdm(pool.starmap(analyze_frame_standalone, args), total=len(traj), desc="Analyzing"))
        else:
            prev_bonds = set() if enable_hysteresis else None
            for i, atoms in enumerate(tqdm(traj, desc="Analyzing")):
                res = analyze_frame_standalone(atoms, i, config_state, prev_bonds)
                results.append(res)
                if enable_hysteresis:
                    prev_bonds = res[4]
        
        results.sort(key=lambda x: x[0])
        
        # Decomposed output logic
        output_dir = Path(output_path).parent
        self._write_counts_csv(results, output_path)
        self._write_bond_distances_csv(results, output_dir / "bond_distances.csv")
        
        if mapping_path is None:
            mapping_path = output_dir / "species_mapping.csv"
        self._write_mapping_csv(results, mapping_path)
        
        logger.info(f"Analysis complete. Results: {output_path}")

    def _write_counts_csv(self, results, path):
        all_species = set()
        for _, counts, _, _, _ in results:
            all_species.update(counts.keys())
        headers = ["Frame"] + sorted(list(all_species))
        
        with open(path, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=headers)
            writer.writeheader()
            for frame, counts, _, _, _ in results:
                row = {"Frame": frame}
                row.update(counts)
                writer.writerow(row)

    def _write_bond_distances_csv(self, results, path):
        with open(path, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(["Frame", "Species", "Distance_A"])
            for frame, _, bond_dists, _, _ in results:
                for entry in bond_dists:
                    writer.writerow([frame, entry['Species'], entry['Distance']])

    def _write_mapping_csv(self, results, path):
        mapping_headers = ["Frame"] + [str(idx) for idx in self.reactive_indices]
        with open(path, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=mapping_headers)
            writer.writeheader()
            for frame, _, _, mapping, _ in results:
                row = {"Frame": frame}
                for idx in self.reactive_indices:
                    row[str(idx)] = mapping.get(int(idx), "Bulk")
                writer.writerow(row)

def main():
    parser = argparse.ArgumentParser(description="Species identification from MD trajectories.")
    parser.add_argument("--config", required=True, help="YAML configuration file.")
    parser.add_argument("--xdatcar", required=True, help="XDATCAR trajectory file.")
    parser.add_argument("--out", required=True, help="Output CSV for counts.")
    parser.add_argument("--nprocs", type=int, default=1, help="Number of processes.")
    
    # Use BooleanOptionalAction for better flag handling (Python 3.9+)
    if hasattr(argparse, 'BooleanOptionalAction'):
        parser.add_argument("--hysteresis", action=argparse.BooleanOptionalAction, default=True, help="Enable/disable bond hysteresis.")
    else:
        parser.add_argument("--hysteresis", action="store_true", default=True)
        parser.add_argument("--no-hysteresis", action="store_false", dest="hysteresis")
        
    args = parser.parse_args()
    
    # Configure logging here instead of module level
    logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
    
    try:
        analyzer = TrajectoryAnalyzer(args.config)
        analyzer.analyze(args.xdatcar, args.out, n_procs=args.nprocs, enable_hysteresis=args.hysteresis)
    except Exception as e:
        logger.error(f"Execution failed: {e}", exc_info=True)
        sys.exit(1)

if __name__ == "__main__":
    main()
