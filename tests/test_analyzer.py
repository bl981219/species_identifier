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
    
    result = analyze_frame_standalone(atoms, 0, config_state)
    
    assert result.counts['O2_gas'] == 1
    assert len(result.bond_distances) == 1
    assert result.bond_distances[0]['Species'] == 'O2'
    assert np.isclose(result.bond_distances[0]['Distance'], 1.2)

def test_adsorbed_species(basic_config):
    analyzer = TrajectoryAnalyzer(basic_config)
    
    # O adsorbed on Pt
    atoms = Atoms('PtO', positions=[
        [0, 0, 0], # Pt
        [0, 0, 1.5] # O (bonded to Pt)
    ], cell=[10, 10, 10], pbc=[True, True, True])
    
    analyzer.setup_static_surface(atoms)
    config_state = analyzer.get_config_state()
    
    result = analyze_frame_standalone(atoms, 0, config_state)
    
    assert result.counts['O*'] == 1

def test_validate_config_missing_key(tmp_path):
    config = {'adsorbates': ['O']} # Missing cutoffs, species_map
    config_path = tmp_path / "bad_config.yaml"
    with open(config_path, 'w') as f:
        yaml.dump(config, f)
    
    with pytest.raises(ValueError, match="Missing required config keys"):
        TrajectoryAnalyzer(str(config_path))

def test_validate_config_invalid_values(tmp_path):
    # Case 1: Non-positive cutoff
    config = {
        'adsorbates': ['O'],
        'cutoffs': {'O-O': -1.0},
        'species_map': {'0,0,1': 'O'}
    }
    config_path = tmp_path / "bad_config_1.yaml"
    with open(config_path, 'w') as f:
        yaml.dump(config, f)
    with pytest.raises(ValueError, match="must be positive"):
        TrajectoryAnalyzer(str(config_path))

    # Case 2: Positive depth cutoff
    config = {
        'adsorbates': ['O'],
        'cutoffs': {'O-O': 1.0},
        'species_map': {'0,0,1': 'O'},
        'depth_cutoffs': {'adsorbates': 0.5}
    }
    config_path = tmp_path / "bad_config_2.yaml"
    with open(config_path, 'w') as f:
        yaml.dump(config, f)
    with pytest.raises(ValueError, match="must be ≤ 0"):
        TrajectoryAnalyzer(str(config_path))

def test_unknown_species_fallback(basic_config):
    analyzer = TrajectoryAnalyzer(basic_config)
    # H3 molecule (not in species_map)
    atoms = Atoms('Pt3H3', positions=[
        [0, 0, 0], [1, 1, 0], [2, 2, 0],
        [0, 0, 5], [0, 0, 5.5], [0.5, 0, 5.5]
    ], cell=[10, 10, 10], pbc=[True, True, True])
    
    analyzer.setup_static_surface(atoms)
    config_state = analyzer.get_config_state()
    result = analyze_frame_standalone(atoms, 0, config_state)
    
    # Formula fallback for H3 should be 'H3'
    assert result.counts['H3_gas'] == 1

def test_hysteresis_stabilizes_bond(basic_config):
    analyzer = TrajectoryAnalyzer(basic_config)
    # O2 molecule, distance is 1.65 (cutoff is 1.6)
    # Frame 0: distance 1.5 (bonded)
    # Frame 1: distance 1.65 (should be bonded with hysteresis 0.15)
    
    atoms0 = Atoms('Pt3O2', positions=[
        [0, 0, 0], [1, 1, 0], [2, 2, 0],
        [0, 0, 5], [0, 0, 6.5] # distance 1.5
    ], cell=[10, 10, 10], pbc=[True, True, True])
    
    analyzer.setup_static_surface(atoms0)
    config_state = analyzer.get_config_state()
    
    res0 = analyze_frame_standalone(atoms0, 0, config_state)
    assert res0.counts['O2_gas'] == 1
    
    atoms1 = Atoms('Pt3O2', positions=[
        [0, 0, 0], [1, 1, 0], [2, 2, 0],
        [0, 0, 5], [0, 0, 6.65] # distance 1.65
    ], cell=[10, 10, 10], pbc=[True, True, True])
    
    # With hysteresis (prev_bonds passed)
    res1_with = analyze_frame_standalone(atoms1, 1, config_state, prev_bonds=res0.new_bonds)
    assert res1_with.counts['O2_gas'] == 1
    
    # Without hysteresis
    res1_without = analyze_frame_standalone(atoms1, 1, config_state, prev_bonds=None)
    # It should split into 2 O atoms
    assert res1_without.counts['O_gas'] == 2

def test_break_surface_bonds_toggle(basic_config):
    # O at z=1.5 (lattice NM) and Pt at z=0 (metal)
    # distance 1.5 < 2.1 (M-O cutoff)
    
    with open(basic_config, 'r') as f:
        config = yaml.safe_load(f)
    config['lattice_non_metal_tolerance'] = 2.0
    with open(basic_config, 'w') as f:
        yaml.dump(config, f)

    atoms = Atoms('PtO', positions=[[0, 0, 0], [0, 0, 1.5]], cell=[10, 10, 10], pbc=[True, True, True])
    
    # Case 1: break_surface_bonds=True (default)
    analyzer_on = TrajectoryAnalyzer(basic_config)
    analyzer_on.setup_static_surface(atoms)
    res_on = analyze_frame_standalone(atoms, 0, analyzer_on.get_config_state())
    assert (0, 1) not in res_on.new_bonds and (1, 0) not in res_on.new_bonds
    
    # Case 2: break_surface_bonds=False
    config['break_surface_bonds'] = False
    with open(basic_config, 'w') as f:
        yaml.dump(config, f)
        
    analyzer_off = TrajectoryAnalyzer(basic_config)
    analyzer_off.setup_static_surface(atoms)
    res_off = analyze_frame_standalone(atoms, 0, analyzer_off.get_config_state())
    assert (0, 1) in res_off.new_bonds or (1, 0) in res_off.new_bonds
