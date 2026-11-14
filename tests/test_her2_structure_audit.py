import numpy as np
import pytest
from biotite.structure import AtomArray
from her2_structure_audit import contacts


def example():
    a = AtomArray(4)
    a.chain_id = np.array(['H','L','T','T'])
    a.res_id = np.array([1,1,10,11])
    a.element = np.array(['C','C','C','H'])
    a.coord = np.array([[0,0,0],[0.1,0,0],[4,0,0],[0.2,0,0]])
    return a


def test_interface_excludes_vh_vl_and_hydrogen():
    r = contacts(example(), ['H','L'], ['T'])
    assert r['heavy_atom_contact_pairs'] == 2
    assert r['heavy_atom_clash_pairs_below_1p5A'] == 0
    assert r['target_contact_residues'] == ['T:10']


def test_chain_validation():
    with pytest.raises(ValueError):
        contacts(example(), ['H'], ['H'])
    with pytest.raises(ValueError):
        contacts(example(), ['H'], ['X'])
