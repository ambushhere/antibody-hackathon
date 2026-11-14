import copy
import json
from pathlib import Path

import numpy as np
import pytest
from biotite.structure import AtomArray

from scripts import audit_structural_replicates as audit


@pytest.fixture
def request_inputs():
    _, _, _, queries, _, _ = audit.frozen_queries()
    request = audit.read_json(audit.RECEIPTS / 'request-seed42.json')
    return request, queries


def test_saved_requests_preserve_all_frozen_inputs_and_parameters(request_inputs):
    _, queries = request_inputs
    jobs = {seed: audit.submission(seed, queries, audit.RECEIPTS)[0] for seed in audit.SEEDS}
    assert len({job['id'] for job in jobs.values()}) == 3


@pytest.mark.parametrize('change, error', [
    ('sequence', 'exact T/H/L FASTA'),
    ('template', 'exact T/H/L FASTA'),
    ('seed', 'protocol mismatch'),
    ('duplicate', 'duplicate or order mismatch'),
])
def test_request_drift_fails(request_inputs, change, error):
    original, queries = request_inputs
    request = copy.deepcopy(original)
    if change == 'sequence':
        request['targets'][0]['content'] = request['targets'][0]['content'].replace('CHQL', 'AHQL', 1)
    elif change == 'template':
        request['targets'][0]['template'] = 'target.cif'
    elif change == 'seed':
        request['params']['seed'] = 43
    else:
        request['targets'].append(request['targets'][0])
    with pytest.raises(ValueError, match=error):
        audit.validate_request(request, 42, queries)


def test_results_coverage_rejects_missing_query_or_wrong_seed():
    row = {'id': 'wt', 'status': 'ok', 'samples': 1, 'msa': False, 'msa_depth': 1, 'seed': 42}
    audit.validate_rows([row], ['wt'], 42)
    with pytest.raises(ValueError, match='missing'):
        audit.validate_rows([row], ['wt', 'strong'], 42)
    with pytest.raises(ValueError, match='seed or MSA mismatch'):
        audit.validate_rows([row], ['wt'], 44)


def test_token_axes_check_full_sequence_and_file_order():
    atoms = AtomArray(4)
    atoms.atom_name[:] = 'CA'
    atoms.chain_id[:] = ['Y', 'X', 'X', 'Z']
    atoms.res_id[:] = [1, 1, 2, 1]
    atoms.res_name[:] = ['ALA', 'CYS', 'TYR', 'MET']
    roles = {'Y': 'L', 'X': 'T', 'Z': 'H'}
    expected = {'T': 'CY', 'H': 'M', 'L': 'A'}
    sidecar = {'model': audit.MODEL, 'n_tokens': 4, 'token_order': audit.TOKEN_ORDER}
    assert list(audit.token_axes(atoms, roles, expected, sidecar, 'fixture')) == ['Y', 'X', 'X', 'Z']
    with pytest.raises(ValueError, match='full CIF sequence'):
        audit.token_axes(atoms, roles, {**expected, 'T': 'WY'}, sidecar, 'fixture')
    with pytest.raises(ValueError, match='token order'):
        audit.token_axes(atoms, roles, expected, {**sidecar, 'token_order': 'input order'}, 'fixture')


def test_confidence_index_rank_is_independent_of_cif_order():
    roles = {'Y': 'L', 'X': 'T', 'Z': 'H'}
    tokens = np.array(['Y', 'X', 'Z'])
    provider = {
        'chains_ptm': {'0': 0.90, '1': 0.39, '2': 0.95},
        'pair_chains_iptm': [[0.90, 0.31, 0.91], [0.31, 0.39, 0.24], [0.91, 0.24, 0.95]],
    }
    arrays = {'pae': np.diag([1., 6., 2.]), 'pde': np.diag([0.2, 4., 0.3])}
    order, evidence = audit.chain_mean_error_order(roles, tokens, arrays, provider, 'fixture')
    assert order == ['H', 'T', 'L']
    assert evidence['status'].startswith('rank inference')
    scores, field = audit.indexed_iptm(provider, {}, order, 'fixture')
    assert field == 'pair_chains_iptm'
    assert [scores[k] for k in ('vh_vl_iptm', 'iptm_H_T', 'iptm_L_T')] == [0.91, 0.31, 0.24]
    arrays['pde'] = np.diag([0.3, 4., 0.2])
    with pytest.raises(ValueError, match='do not resolve a unique'):
        audit.chain_mean_error_order(roles, tokens, arrays, provider, 'fixture')


def test_pae_array_rejects_size_mismatch_and_nonfinite_values(tmp_path):
    path = tmp_path / 'confidence.npz'
    sidecar = {'arrays': {'pae': {'shape': [2, 2], 'dtype': 'float32', 'units': 'angstrom'}}}
    small = np.array([[0, 1], [2, 0]], dtype=np.float32)
    np.savez(path, pae=small)
    np.testing.assert_array_equal(audit.error_arrays(path, sidecar, 2)['pae'], small)
    with pytest.raises(ValueError, match='invalid pae'):
        audit.error_arrays(path, sidecar, 3)
    small[0, 1] = np.nan
    np.savez(path, pae=small)
    with pytest.raises(ValueError, match='invalid pae'):
        audit.error_arrays(path, sidecar, 2)


def test_download_receipt_binds_hashes_to_job_and_file(tmp_path):
    folder = tmp_path / 'seed-42'
    folder.mkdir()
    artifact = folder / 'control_wt.cif'
    artifact.write_text('original')
    path = tmp_path / 'receipt.json'
    path.write_text(json.dumps({'seed': 42, 'job_id': 'job42', 'files': [{
        'path': 'seed-42/control_wt.cif', 'sha256': audit.file_hash(artifact), 'bytes': 8,
    }]}))
    hashes = {artifact.name: audit.file_hash(artifact)}
    audit.verify_download_hashes(path, 42, 'job42', folder, hashes)
    with pytest.raises(ValueError, match='seed or job mismatch'):
        audit.verify_download_hashes(path, 42, 'job-other', folder, hashes)
    artifact.write_text('modified')
    hashes[artifact.name] = audit.file_hash(artifact)
    with pytest.raises(ValueError, match='SHA-256 or size mismatch'):
        audit.verify_download_hashes(path, 42, 'job42', folder, hashes)


def test_analogous_rf3_fixture_reuses_exact_chain_and_geometry_logic():
    folder = Path('results/campaign/japanfold-rf3-domain-iv')
    reference = Path('/tmp/her2-ablang/assets/1N8Z.cif')
    if not reference.is_file() or not folder.is_dir():
        pytest.skip('Local 1N8Z and extracted RF3 archive required')
    manifest, panel, canonical, queries, maps, _ = audit.frozen_queries()
    name = 'control_wt'
    refs = audit.read_reference(reference, manifest, canonical, queries[name], maps)
    atoms, roles = audit.validate_prediction(folder / f'{name}.cif', queries[name], protenix=False)
    sidecar = audit.read_json(folder / f'{name}_pae.json')
    tokens = audit.token_axes(atoms, roles, queries[name], {**sidecar, 'model': audit.MODEL}, name)
    arrays = audit.error_arrays(folder / f'{name}_pae.npz', sidecar, len(tokens))
    summary = audit.read_json(folder / f'{name}_summary_confidences.json')
    assert audit.rf3_summary_order(roles, tokens, arrays, summary, name) == ['H', 'L', 'T']
    mapped, geometry = audit.geometry_row(atoms, roles, refs, maps['domain-iv'], True)
    assert (len(refs[4]), len(refs[5]), geometry['wt_contact_count']) == (19, 109, 14)
    assert audit.footprint_pae(arrays['pae'], mapped, roles, refs[4], refs[3])
    assert not audit.antibody_mapping('denovo_fullseq_k1_1', panel)


def test_full_opendde_artifact_audit_when_local():
    reference = Path('/tmp/her2-ablang/assets/1N8Z.cif')
    if not reference.is_file() or any(not (audit.CAMPAIGN / f'seed-{seed}').is_dir() for seed in audit.SEEDS):
        pytest.skip('Local 1N8Z and all three OpenDDE archives required')
    result = audit.audit(reference)
    rows = result['predictions']
    assert len(rows) == 21
    assert {(row['seed'], row['query_id']) for row in rows} == {
        (seed, name) for seed in audit.SEEDS for name in audit.PANEL
    }
    assert all(row['target_fit_ca_count'] == 109 for row in rows)
    assert all('antibody_ca_rmsd_in_target_frame_A' not in row for row in rows if row['query_id'].startswith('denovo'))
    assert [row['wt_recovered'] for row in result['controls_by_seed']] == [2, 0, 5]
