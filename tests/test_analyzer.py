import pytest
import numpy as np
import yaml
from ase.atoms import Atoms
from species_identifier.analyzer import TrajectoryAnalyzer, analyze_frame_standalone

@pytest.fixture
def basic_config(tmp_path):
    config = {
        'adsorbates': ['C', 'H', 'O'],
        'lattice_non_metals': ['O'],
        'cutoffs': {
            'H-H': 1.0,
            'O-O': 1.6,
            'M-O': 2.1
        },
        'species_map': {
            '0,0,2': 'O2',
            '0,0,1': 'O'
        },
        'depth_cutoffs': {
            'adsorbates': -0.5,
            'surface_metals': -2.0,
            'lattice_non_metals': -2.0
        },
        'track_distances': ['O2']
    }
    config_path = tmp_path / "config.yaml"
    with open(config_path, 'w') as f:
        yaml.dump(config, f)
    return str(config_path)

def test_cutoff_normalization(basic_config):
    analyzer = TrajectoryAnalyzer(basic_config)
    analyzer._prepare_cutoff_matrix(['H', 'O', 'Pt'])
    # H-H should be 1.0
    assert analyzer.cutoff_matrix.get(('H', 'H')) == 1.0
    # O-O should be 1.6
    assert analyzer.cutoff_matrix.get(('O', 'O')) == 1.6

def test_analyze_frame_basic(basic_config):
    analyzer = TrajectoryAnalyzer(basic_config)
    
    # Pt at z=0, O2 at z=3.0
    atoms = Atoms('Pt3O2', positions=[
        [0, 0, 0], [1, 1, 0], [2, 2, 0], # Pt
        [0, 0, 3], [0, 0, 4.2] # O2 (1.2 A apart)
    ], cell=[10, 10, 10], pbc=[True, True, True])
    
    analyzer.setup_static_surface(atoms)
    config_state = analyzer.get_config_state()
    
    frame_idx, counts, bond_dists, mapping, bonds = analyze_frame_standalone(atoms, 0, config_state)
    
    assert counts['O2_gas'] == 1
    assert len(bond_dists) == 1
    assert bond_dists[0]['Species'] == 'O2'
    assert np.isclose(bond_dists[0]['Distance'], 1.2)

def test_adsorbed_species(basic_config):
    analyzer = TrajectoryAnalyzer(basic_config)
    
    # O adsorbed on Pt
    atoms = Atoms('PtO', positions=[
        [0, 0, 0], # Pt
        [0, 0, 1.5] # O (bonded to Pt)
    ], cell=[10, 10, 10], pbc=[True, True, True])
    
    analyzer.setup_static_surface(atoms)
    config_state = analyzer.get_config_state()
    
    frame_idx, counts, bond_dists, mapping, bonds = analyze_frame_standalone(atoms, 0, config_state)
    
    assert counts['O*'] == 1
