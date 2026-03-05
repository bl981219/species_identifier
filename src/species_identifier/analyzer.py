#!/usr/bin/env python3
import csv
import sys
import yaml
import argparse
import logging
import numpy as np
from collections import Counter
from pathlib import Path
from typing import List, Dict, Set, Tuple, Any, Optional, FrozenSet
from ase.io import read
from ase.atoms import Atoms
from ase.neighborlist import neighbor_list
from tqdm import tqdm

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")

class TrajectoryAnalyzer:
    def __init__(self, config_path: str):
        self.config_path = Path(config_path)
        if not self.config_path.exists():
            raise FileNotFoundError(f"Config file not found: {config_path}")
            
        with open(self.config_path, 'r') as f:
            self.config: Dict[str, Any] = yaml.safe_load(f)
        
        # Parse cutoffs into frozensets for easy matching
        self.cutoffs: Dict[FrozenSet[str], float] = {
            frozenset(k.split('-')): float(v) for k, v in self.config.get('cutoffs', {}).items()
        }
        self.max_cutoff: float = max(self.cutoffs.values()) if self.cutoffs else 0.0
        self.hysteresis: float = self.config.get('bond_hysteresis', 0.15)
        self.adsorbates: List[str] = self.config.get('adsorbates', [])
        self.lattice_nms: List[str] = self.config.get('lattice_non_metals', [])
        self.species_map: Dict[str, str] = self.config.get('species_map', {})
        self.break_surface: bool = self.config.get('break_surface_bonds', True)
        self.pbc: List[bool] = self.config.get('pbc', [True, True, False])
        
        # State tracking
        self.prev_bonds: Set[Tuple[int, int]] = set()
        self.metals: Set[str] = set()
        self.surface_z: float = 0.0
        self.is_lat_nm: Dict[int, bool] = {}
        self.reactive_indices: List[int] = []

    def get_bond_cutoff(self, s1: str, s2: str, pair_in_prev: bool) -> float:
        """Determines bond distance cutoff with optional hysteresis."""
        pair = frozenset([s1, s2])
        r_on = self.cutoffs.get(pair, 0.0)
        
        # Fallback to metal-generic keys (M-X)
        if r_on <= 0:
            k1 = 'M' if s1 in self.metals else s1
            k2 = 'M' if s2 in self.metals else s2
            r_on = self.cutoffs.get(frozenset([k1, k2]), 0.0)
            
        if r_on <= 0: return 0.0
        return r_on + (self.hysteresis if pair_in_prev else 0.0)

    def get_species_name(self, mol_syms: List[str]) -> str:
        """Translates atom counts into a chemical name or formula."""
        c = Counter(mol_syms)
        # Key generated based on the order of self.adsorbates in config
        key = ",".join(str(c.get(el, 0)) for el in self.adsorbates)
        if key in self.species_map:
            return self.species_map[key]
        
        # Fallback to basic chemical formula
        formula = "".join(f"{el}{c[el] if c[el] > 1 else ''}" for el in self.adsorbates if c[el] > 0)
        return formula if formula else "Unknown"

    def setup_surface_mask(self, atoms: Atoms) -> np.ndarray:
        """Identifies metal surface height and defines reactive atom indices."""
        syms = np.array(atoms.get_chemical_symbols())
        pos = atoms.get_positions(wrap=True)
        self.metals = set(syms) - set(self.adsorbates) - set(self.lattice_nms)
        
        m_mask = [s in self.metals for s in syms]
        m_z = pos[m_mask, 2]
        
        if len(m_z) == 0:
            logging.warning("No metal atoms found; using Z=0.0 as surface reference.")
            self.surface_z = 0.0
        else:
            self.surface_z = float(np.percentile(m_z, self.config.get('surface_metal_z_percentile', 98.0)))
        
        lat_nm_tol = self.config.get('lattice_non_metal_tolerance', 0.7)
        lat_nm_indices = {i for i, (s, z) in enumerate(zip(syms, pos[:, 2])) 
                         if s in self.lattice_nms and z < (self.surface_z + lat_nm_tol)}
        
        self.is_lat_nm = {i: (i in lat_nm_indices) for i in range(len(syms))}
        depths = self.config.get('depth_cutoffs', {})
        
        self.reactive_indices = [i for i, (s, z) in enumerate(zip(syms, pos[:, 2])) if 
            (s in self.metals and z > self.surface_z + depths.get('surface_metals', -2.0)) or
            (i in lat_nm_indices and z > self.surface_z + depths.get('lattice_non_metals', -2.0)) or
            (s in self.adsorbates and z > self.surface_z + depths.get('adsorbates', -0.5))]
        
        return np.array(self.reactive_indices)

    def build_graph(self, atoms: Atoms, reactive_indices: np.ndarray) -> Dict[int, Set[int]]:
        """Constructs an adjacency dictionary for reactive atoms."""
        syms = np.array(atoms.get_chemical_symbols())
        i_list, j_list, d_list = neighbor_list('ijd', atoms, self.max_cutoff + self.hysteresis)
        
        adj: Dict[int, Set[int]] = {int(g): set() for g in reactive_indices}
        new_bonds: Set[Tuple[int, int]] = set()
        r_set = set(reactive_indices)
        
        for i, j, d in zip(i_list, j_list, d_list):
            if i >= j or i not in r_set or j not in r_set: continue
            
            si, sj = syms[i], syms[j]
            if self.break_surface:
                # Prevent artificial bonding between metal substrate and lattice oxygens
                if (si in self.metals and self.is_lat_nm.get(j)) or (sj in self.metals and self.is_lat_nm.get(i)):
                    continue
            
            pair = (i, j)
            r_limit = self.get_bond_cutoff(si, sj, pair in self.prev_bonds)
            
            if d < r_limit:
                adj[i].add(j)
                adj[j].add(i)
                new_bonds.add(pair)
        
        self.prev_bonds = new_bonds
        return adj

    def get_components(self, adj: Dict[int, Set[int]], nodes: List[int]) -> List[List[int]]:
        """Finds connected components in the graph using BFS."""
        visited: Set[int] = set()
        components: List[List[int]] = []
        for node in nodes:
            if node in visited: continue
            comp, stack = [], [node]
            visited.add(node)
            while stack:
                curr = stack.pop()
                comp.append(curr)
                for n in adj.get(curr, []):
                    if n not in visited:
                        visited.add(n)
                        stack.append(n)
            components.append(comp)
        return components

    def analyze(self, xdatcar_path: str, output_path: str):
        """Main analysis loop across all trajectory frames."""
        xdatcar = Path(xdatcar_path)
        if not xdatcar.exists():
            raise FileNotFoundError(f"Trajectory file not found: {xdatcar_path}")
            
        logging.info(f"Loading trajectory from {xdatcar_path}...")
        traj = read(xdatcar_path, index=':')
        
        # Initial surface establishment
        reactive = self.setup_surface_mask(traj[0])
        all_labels: Set[str] = set()
        results: List[Dict[str, Any]] = []

        logging.info(f"Surface established at Z={self.surface_z:.2f} Å. Analyzing {len(traj)} frames...")
        
        for i, atoms in enumerate(tqdm(traj, desc="Processing Frames", unit="frame")):
            atoms.pbc = self.pbc
            if self.config.get('dynamic_surface', False):
                reactive = self.setup_surface_mask(atoms)
                
            adj = self.build_graph(atoms, reactive)
            frame_counts: Counter = Counter()
            cur_syms = np.array(atoms.get_chemical_symbols())
            
            for mol in self.get_components(adj, list(reactive)):
                # Distinguish if connected to slab (adsorbed) or not (gas)
                is_ads = any(cur_syms[j] in self.metals or self.is_lat_nm.get(j) for j in mol)
                fragment = [j for j in mol if cur_syms[j] not in self.metals and not self.is_lat_nm.get(j)]
                
                if not fragment: continue
                
                name = self.get_species_name(list(cur_syms[fragment]))
                label = f"{name}*" if is_ads else f"{name}_gas"
                frame_counts[label] += 1
                all_labels.add(label)
            
            results.append({"Frame": i + 1, **frame_counts})

        # Save to CSV
        with open(output_path, 'w', newline='') as f:
            fieldnames = ["Frame"] + sorted(list(all_labels))
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            for r in results:
                writer.writerow(r)
        
        logging.info(f"Analysis complete. Results saved to {output_path}")

def main():
    parser = argparse.ArgumentParser(description="Graph-Based Species Identification Tool")
    parser.add_argument("--config", required=True, help="Path to config.yaml")
    parser.add_argument("--xdatcar", default="XDATCAR", help="MD Trajectory file")
    parser.add_argument("--out", default="species_results.csv", help="Output CSV name")
    args = parser.parse_args()

    try:
        analyzer = TrajectoryAnalyzer(args.config)
        analyzer.analyze(args.xdatcar, args.out)
    except Exception as e:
        logging.error(f"Analysis failed: {e}", exc_info=True)
        sys.exit(1)

if __name__ == "__main__":
    main()