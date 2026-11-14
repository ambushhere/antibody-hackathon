"""Audit the frozen HER2 context panel from local artifacts on CPU."""
import argparse
import csv
import gzip
import hashlib
import itertools
import json
import re
from pathlib import Path

import numpy as np
from biotite.sequence import ProteinSequence
from biotite.structure import get_residues
from biotite.structure.io.pdbx import CIFFile, get_structure

from scripts.prepare_her2_context import canonical_sequence, require, sha256
from src.her2_structure_audit import contacts, remap_target_residues

CONFIG = Path('configs/her2-context')
CAMPAIGN = Path('results/campaign')
PROTENIX = Path('campaign-20261114/context-panel')
RF3 = Path('japanfold-rf3-domain-iv')
ARCHIVE = Path('results/archives/structure-runs-manifest.json')
RECEIPTS = Path('results/provider-jobs')
ARMS = {
    'split': {'C': range(511, 603), 'D': range(613, 630)},
    'domain-iv': {'T': range(511, 653)},
    'ecd': {'T': range(23, 653)},
}
PANEL = (
    'control_wt', 'control_spr_strong', 'control_spr_weak', 'control_facs_low',
    'denovo_fullseq_k1_1', 'opt_WTAVGFFAMD', 'opt_YTGDGFPAMN',
)
SEED = 42
SAMPLE = 1
IMAGE_ID = 'sha256:d0d1d4a8a911b8001a0868aa30b9f88f598ce12236cdd0b75a4b3465e7a3da58'
REFERENCE_CHAINS = {'H': 'B', 'L': 'A'}


def read_json(path):
    return json.loads(path.read_text())


def file_hash(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def exact_coverage(actual, expected, label, ordered=False):
    actual, expected = list(actual), list(expected)
    missing, extra = set(expected) - set(actual), set(actual) - set(expected)
    require(not missing and not extra and len(actual) == len(expected)
            and (actual == expected or not ordered),
            f'{label}: missing {sorted(missing)}, extra {sorted(extra)}, '
            f'duplicate or order mismatch if counts differ ({len(actual)}/{len(expected)}) '
            'or required order differs')


def validate_file(path, expected_hash, expected_bytes=None):
    require(path.is_file(), f'Missing artifact {path}')
    require(expected_bytes is None or path.stat().st_size == expected_bytes, f'Wrong size for {path}')
    actual = file_hash(path)
    require(actual == expected_hash, f'SHA-256 mismatch for {path}')
    return actual


def validate_archive():
    archive = read_json(ARCHIVE)
    archive_hash = validate_file(ARCHIVE.parent / archive['archive'], archive['sha256'], archive['bytes'])
    hashes = {}
    for prefix in (PROTENIX, RF3):
        members = [
            item for item in archive['files']
            if Path(item['path']).is_relative_to(prefix)
        ]
        entries = {Path(item['path']): item for item in members}
        require(len(entries) == len(members), f'{prefix} has duplicate archive manifest paths')
        actual = {p.relative_to(CAMPAIGN) for p in (CAMPAIGN / prefix).rglob('*') if p.is_file()}
        # The duplicate download zip is outside the archive and is not consumed.
        actual.discard(RF3 / 'results.zip')
        exact_coverage(actual, entries, f'{prefix} raw files')
        for path, item in sorted(entries.items()):
            hashes[str(path)] = validate_file(CAMPAIGN / path, item['sha256'], item['bytes'])
    return {'manifest_sha256': file_hash(ARCHIVE), 'archive_sha256': archive_hash, 'files': hashes}


def validate_inputs():
    manifest = read_json(CONFIG / 'manifest.json')
    staged = {}
    for line in (CAMPAIGN / PROTENIX / 'input.sha256').read_text().splitlines():
        digest, path = line.split(maxsplit=1)
        name = Path(path).name
        require(name not in staged, f'Duplicate staged source {name}')
        staged[name] = digest
    exact_coverage(staged, [*manifest['artifact_sha256'], 'manifest.json'], 'staged input files')
    for name, digest in staged.items():
        require(manifest['artifact_sha256'].get(name, digest) == digest, f'Manifest disagrees for {name}')
        validate_file(CONFIG / name, digest)
    canonical = canonical_sequence((CONFIG / 'P04626.fasta').read_bytes())
    require(sha256(canonical.encode()) == manifest['canonical_sequence_sha256'],
            'Unexpected canonical sequence hash')
    exact_coverage(manifest['arms'], ARMS, 'manifest arms')
    panel = {row['query_id']: row for row in manifest['panel']}
    exact_coverage([row['query_id'] for row in manifest['panel']], PANEL, 'manifest panel')
    queries, mappings = {}, {}
    for arm, ranges in ARMS.items():
        entry = manifest['arms'][arm]
        exact_coverage(entry['target_chains'], ranges, f'{arm} target maps')
        mapping = {chain: m['query_position_to_canonical'] for chain, m in entry['target_chains'].items()}
        targets = {}
        for chain, positions in mapping.items():
            details = entry['target_chains'][chain]
            require(positions == list(ranges[chain]), f'{arm}/{chain} canonical numbering changed')
            targets[chain] = ''.join(canonical[p - 1] for p in positions)
            require(len(positions) == details['length']
                    and sha256(targets[chain].encode()) == details['sequence_sha256'],
                    f'{arm}/{chain} target sequence hash or length mismatch')
        expected_order = [*targets, 'H', 'L']
        items = read_json(CONFIG / entry['queries'])['queries']
        exact_coverage([q['query_id'] for q in items], panel, f'{arm} queries', ordered=True)
        queries[arm], mappings[arm] = {}, mapping
        for q in items:
            name = q['query_id']
            exact_coverage([c['chain_id'] for c in q['chains']], expected_order,
                           f'{arm}/{name} input chain order', ordered=True)
            expected = dict(targets)
            for chain, kind in [('H', 'vh'), ('L', 'vl')]:
                expected[chain] = panel[name][f'{kind}_sequence']
                require(sha256(expected[chain].encode()) == panel[name][f'{kind}_sha256'],
                        f'{name} {kind} sequence hash mismatch')
            require(all(c['sequence'] == expected[c['chain_id']] and c['msa_path'] is None
                        for c in q['chains']), f'{arm}/{name} query sequence or MSA mismatch')
            queries[arm][name] = expected
    require(manifest['sifts']['mature_to_canonical_offset'] == 22, 'Unexpected 1N8Z SIFTS offset')
    return manifest, panel, canonical, queries, mappings, staged


def validate_protenix_run(arm, queries, entry):
    folder = CAMPAIGN / PROTENIX / arm
    run = read_json(folder / 'result.json')
    provenance = run['provenance']
    config = read_json(folder / 'config.json')
    require(run['task'] == 'predict' and run['schema_version'] == 'foundry.cli.v1'
            and provenance['module'] == {'name': 'apheris-protenix', 'version': '1.3.0'}
            and config == {'schema_version': 'foundry.cli.v1', **provenance['resolved_config']},
            f'{arm} runner metadata mismatch')
    require(config['seeds'] == [SEED] and config['num_diffusion_samples'] == 1
            and config['write_detailed_scores'] is True and config['query_failures'] == 'fail'
            and config['missing_msas'] == 'fail' and config['weights'] == '/app/weights/protenix-v2.pt'
            and len(config['template_paths']) == 1
            and Path(config['template_paths'][0]).name == entry['template']
            and len(provenance['inputs']) == 1
            and Path(provenance['inputs'][0]['path']).name == entry['queries'],
            f'{arm} model, seed, sample or source mismatch')
    predictions = run['result']['predictions']
    exact_coverage([p['query_id'] for p in predictions], queries, f'{arm} prediction receipt', ordered=True)
    suffixes = {'path': '_model.cif', 'scores_path': '_scores.json',
                'detailed_scores_path': '_detailed_scores.cbor.gz'}
    expected_paths = {
        key: [f'{name}/seed_{SEED}/{name}_seed_{SEED}_sample_{SAMPLE}{suffix}' for name in queries]
        for key, suffix in suffixes.items()
    }
    for i, pred in enumerate(predictions):
        require((pred['query_index'], pred['seed'], pred['sample']) == (i, SEED, SAMPLE),
                f"{arm}/{pred['query_id']} query index, seed or sample mismatch")
    for key, suffix in suffixes.items():
        exact_coverage([p[key] for p in predictions], expected_paths[key], f'{arm} receipt {key}', ordered=True)
        exact_coverage([str(p.relative_to(folder)) for p in folder.rglob(f'*{suffix}')],
                       expected_paths[key], f'{arm} {suffix}')
    log = (folder.with_suffix('.log')).read_text()
    assigned = re.findall(r'Assigned templates to (\d+) entities across (\d+) queries', log)
    resolved = set(re.findall(r'Identified (\d+) resolved residue\(s\) from', log))
    lengths = set(re.findall(r'Read \d+ chains? with a total of (\d+) residues from `_pdbx_poly_seq_scheme`', log))
    expected_resolved = set().union(*(set(m['query_position_to_canonical'])
                                     - {p + 22 for p in m['unresolved_deposited_mature']}
                                     - {p + 22 for p in m['beyond_deposit_mature']}
                                     for m in entry['target_chains'].values()))
    total = sum(m['length'] for m in entry['target_chains'].values())
    require(assigned == [(str(len(entry['target_chains']) * len(queries)), str(len(queries)))]
            and resolved == {str(len(expected_resolved))} and lengths == {str(total)}
            and f'Seed set to {SEED}' in log and f'{len(queries)}/{len(queries)}' in log
            and 'Loaded constraint feature: #atom contact:0 #contact:0 #pocket:0' in log,
            f'{arm} template assignment, template coverage or completion log mismatch')
    return {p['query_id']: {key: folder / p[key] for key in suffixes} for p in predictions}, {
        'module': provenance['module'], 'config': config, 'template_assignments': int(assigned[0][0]),
        'template_resolved_residues_in_log': len(expected_resolved),
        'input_queries': entry['queries'], 'input_template': entry['template'], 'queries_completed': len(queries),
    }


def validate_rf3_run(queries):
    folder = CAMPAIGN / RF3
    job_path, result_path, usage_path = [RECEIPTS / f'{RF3.name}{suffix}.json'
                                        for suffix in ('', '-results', '-usage')]
    job, receipt, usage = (read_json(path) for path in (job_path, result_path, usage_path))
    require(job['model'] == usage['model'] == 'rf3' and job['kind'] == 'predict'
            and job['status'] == usage['status'] == 'succeeded' and job['error'] is None
            and job['id'] == receipt['job_id'] == usage['job']
            and job['results_ready'] is True and receipt['ready'] is True
            and job['done'] == job['total'] == len(queries)
            and job['params'] == {'diffusion_samples': 1, 'output_format': 'cif', 'seed': SEED,
                                  'use_msa_server': False, 'write_pae': True},
            'RF3 provider model, job completion, seed or sample mismatch')
    for items, key, label in [(job['targets'], 'id', 'job targets'), (receipt['rows'], 'id', 'results receipt'),
                              (usage['targets'], 'id', 'usage targets')]:
        exact_coverage([item[key] for item in items], queries, f'RF3 {label}')
    require(all(row['state'] == 'done' for row in job['targets'])
            and all(row['attempts'] == 1 for row in usage['targets']), 'RF3 target completion or attempt mismatch')
    local_rows = read_json(folder / 'results.json')
    require(local_rows == receipt['rows'], 'RF3 archive results differ from the provider receipt')
    structures = {name: [f'{name}.cif'] for name in queries}
    require(receipt['structures'] == structures, 'RF3 receipt structure mapping mismatch')
    expected = {'structure': [f'{name}.cif' for name in queries],
                'pae': [f'{name}_pae.npz' for name in queries],
                'pae_sidecar': [f'{name}_pae.json' for name in queries]}
    exact_coverage({item['type'] for item in receipt['artifacts']}, expected, 'RF3 receipt artifact types')
    unique_artifacts = {item['path']: item for item in receipt['artifacts']}
    require(all(item == unique_artifacts[item['path']] for item in receipt['artifacts']),
            'RF3 receipt repeats artifact paths with conflicting metadata')
    duplicates = sorted(path for path in unique_artifacts
                        if sum(item['path'] == path for item in receipt['artifacts']) > 1)
    require(set(duplicates) <= set(expected['pae_sidecar']), 'RF3 repeats non-sidecar artifact paths')
    for kind, files in expected.items():
        exact_coverage([item['path'] for item in unique_artifacts.values() if item['type'] == kind],
                       files, f'RF3 {kind} receipt files')
    for kind in ('pae', 'pae_sidecar'):
        paths = list(dict.fromkeys(receipt[kind])) if kind == 'pae_sidecar' else receipt[kind]
        exact_coverage(paths, expected[kind], f'RF3 {kind} receipt list')
    for item in receipt['artifacts']:
        if item['type'] == 'structure':
            require(structures.get(item['target']) == [item['path']], 'RF3 artifact target mismatch')
    for suffix in ('.cif', '_summary_confidences.json', '_pae.npz', '_pae.json'):
        exact_coverage([p.name for p in folder.glob(f'*{suffix}')],
                       [name + suffix for name in queries], f'RF3 {suffix}')
    require(all(row['status'] == 'ok' and row['samples'] == 1 and row['msa'] is False and row['msa_depth'] == 1
                for row in local_rows), 'RF3 target failed or used a different sample/MSA policy')
    return {row['id']: row for row in local_rows}, {
        'model': job['model'], 'id': job['id'], 'params': job['params'], 'queries_completed': job['done'],
        'receipt_sha256': {path.name: file_hash(path) for path in (job_path, result_path, usage_path)},
        'duplicate_receipt_sidecar_paths': duplicates,
        'template': 'none supplied according to protocol; receipt has no template parameter',
        'query_sequences': 'validated from outputs against domain-iv queries; submission payload unavailable',
        'model_revision': 'not exposed in receipts',
    }


def chain_sequence(atoms):
    _, names = get_residues(atoms)
    return ''.join(ProteinSequence.convert_letter_3to1(name) for name in names)


def sequence_roles(observed, expected, path, require_query_order=False):
    require(len(set(expected.values())) == len(expected), f'{path}: ambiguous expected chain sequences')
    by_sequence = {seq: role for role, seq in expected.items()}
    require(len(observed) == len(expected) and sorted(observed.values()) == sorted(expected.values()),
            f'{path}: chain sequence mismatch')
    roles = {chain: by_sequence[seq] for chain, seq in observed.items()}
    if require_query_order:
        require(list(roles) == list(roles.values()) == list(expected), f'{path}: chain IDs or order mismatch')
    return roles


def validate_prediction(path, expected, protenix):
    cif = CIFFile.read(path)
    site = cif.block['atom_site']
    require(set(site['pdbx_PDB_model_num'].as_array()) == {'1'}, f'{path}: expected one structure model')
    require(np.array_equal(site['label_asym_id'].as_array(), site['auth_asym_id'].as_array())
            and np.array_equal(site['label_seq_id'].as_array(), site['auth_seq_id'].as_array()),
            f'{path}: author and label numbering differ')
    atoms = get_structure(cif, model=1, use_author_fields=False, extra_fields=['b_factor'])
    require(len(atoms) == len(site['id'].as_array()) and not atoms.hetero.any()
            and not np.any(atoms.ins_code != '') and np.isfinite(atoms.coord).all(),
            f'{path}: unexpected coordinate records')
    chains = tuple(dict.fromkeys(atoms.chain_id))
    seqs = {}
    for chain in chains:
        selected = atoms[atoms.chain_id == chain]
        ids, _ = get_residues(selected)
        require(np.array_equal(ids, np.arange(1, len(ids) + 1))
                and np.array_equal(selected.res_id[selected.atom_name == 'CA'], ids),
                f'{path}: missing, duplicate or nonsequential residues in {chain}')
        seqs[str(chain)] = chain_sequence(selected)
    require(list(dict.fromkeys(zip(atoms.chain_id, atoms.res_id))) ==
            [(chain, index) for chain, seq in seqs.items() for index in range(1, len(seq) + 1)]
            and len(set(zip(atoms.chain_id, atoms.res_id, atoms.atom_name))) == len(atoms),
            f'{path}: unexpected atom order or duplicate atom keys')
    roles = sequence_roles(seqs, expected, path, require_query_order=protenix)
    return atoms, roles


def canonical_target(atoms, roles, mapping):
    targets = {chain: mapping[role] for chain, role in roles.items() if role in mapping}
    # Models have already been checked for complete 1-based label numbering.
    return remap_target_residues(atoms, targets)


def ca_coordinates(atoms):
    alpha = atoms[atoms.atom_name == 'CA']
    return {int(pos): xyz for pos, xyz in zip(alpha.res_id, alpha.coord)}


def rmsd(left, right):
    return float(np.sqrt(np.mean(np.sum((left - right) ** 2, axis=1))))


def target_transform(moving, fixed):
    moving_mean, fixed_mean = moving.mean(axis=0), fixed.mean(axis=0)
    u, _, vh = np.linalg.svd((moving - moving_mean).T @ (fixed - fixed_mean))
    rotation = u @ np.diag([1, 1, np.linalg.det(u @ vh)]) @ vh
    return rotation, fixed_mean - moving_mean @ rotation


def read_reference(path, manifest, canonical, wt, mappings):
    source = manifest['sources']['pdb_cif']
    digest = validate_file(path, source['sha256'])
    atoms = get_structure(CIFFile.read(path), model=1)
    atoms = atoms[~atoms.hetero]
    target = atoms[atoms.chain_id == 'C']
    target.res_id += 22
    ref_target_ids, names = get_residues(target)
    require(all(ProteinSequence.convert_letter_3to1(name) == canonical[pos - 1]
                for pos, name in zip(ref_target_ids, names)) and not np.any(target.ins_code != ''),
            '1N8Z target sequence disagrees with SIFTS numbering')
    require(len(target) == manifest['audit']['rcsb_target_polymer_atoms']
            and len(ref_target_ids) == manifest['audit']['rcsb_resolved_target_residues'],
            '1N8Z target coordinate count mismatch')
    split = set().union(*map(set, mappings['split'].values()))
    ca = ca_coordinates(target)
    require(split <= ca.keys(), 'Incomplete reference C-alpha coverage of split target')
    antibody = {}
    for role, chain in REFERENCE_CHAINS.items():
        selected = atoms[(atoms.chain_id == chain) & (atoms.res_id <= len(wt[role]))]
        ids, _ = get_residues(selected)
        require(list(ids) == list(range(1, len(wt[role]) + 1)) and not np.any(selected.ins_code != '')
                and chain_sequence(selected) == wt[role] and ca_coordinates(selected).keys() == set(ids),
                f'1N8Z {chain} variable sequence, numbering or C-alpha coverage mismatch')
        antibody[role] = selected
    fab = atoms[np.isin(atoms.chain_id, list(REFERENCE_CHAINS.values()))]
    geometry = contacts(target + fab, list(REFERENCE_CHAINS.values()), ['C'])
    fv = contacts(target + antibody['H'] + antibody['L'], list(REFERENCE_CHAINS.values()), ['C'])
    footprint = {int(value.split(':')[1]) for value in geometry['target_contact_residues']}
    require(footprint == set(manifest['audit']['footprint_canonical_positions'])
            and geometry['target_contact_residues'] == fv['target_contact_residues'] and footprint <= split,
            'Deposited Fab and Fv footprint mismatch')
    antibody_sites = {role: sorted(int(value.split(':')[1]) for value in fv['antibody_contact_residues']
                                   if value.split(':')[0] == chain)
                      for role, chain in REFERENCE_CHAINS.items()}
    return target, ca, antibody, antibody_sites, footprint, sorted(split), digest


def coordinate_keys(atoms):
    keys = list(zip(atoms.res_id, atoms.res_name, atoms.atom_name))
    require(len(set(keys)) == len(keys), 'Duplicate template or reference atom keys')
    return {key: xyz for key, xyz in zip(keys, atoms.coord)}


def audit_templates(manifest, mappings, targets, reference):
    deposit = coordinate_keys(reference[reference.atom_name != 'OXT'])
    summaries, key_sets = {}, {}
    for arm, entry in manifest['arms'].items():
        template = CONFIG / entry['template']
        cif = CIFFile.read(template)
        atoms = get_structure(cif, model=1, use_author_fields=False)
        chains = entry['target_chains']
        expected = targets[arm]['control_wt']
        exact_coverage(dict.fromkeys(atoms.chain_id), chains, f'{arm} template chain IDs', ordered=True)
        poly = cif.block['entity_poly']
        seqs = dict(zip(poly['entity_id'].as_array(), poly['pdbx_seq_one_letter_code_can'].as_array()))
        asym = cif.block['struct_asym']
        entity_by_chain = dict(zip(asym['id'].as_array(), asym['entity_id'].as_array()))
        remapped = atoms.copy()
        for chain, mapping in mappings[arm].items():
            mask = atoms.chain_id == chain
            ids, names = get_residues(atoms[mask])
            expected_ids = [i for i, pos in enumerate(mapping, start=1)
                            if pos - 22 not in chains[chain]['unresolved_deposited_mature']
                            and pos - 22 not in chains[chain]['beyond_deposit_mature']]
            require(''.join(seqs[entity_by_chain[chain]].split()) == expected[chain]
                    and list(ids) == expected_ids
                    and all(ProteinSequence.convert_letter_3to1(name) == expected[chain][i - 1]
                            for i, name in zip(ids, names)),
                    f'{arm}/{chain} template polymer sequence or coordinate index/gaps differ')
            remapped.res_id[mask] = np.asarray(mapping)[atoms.res_id[mask] - 1]
        keys = coordinate_keys(remapped)
        require(np.array_equal(remapped.res_id - 22, cif.block['atom_site']['auth_seq_id'].as_array(int)),
                f'{arm} template author and canonical positions differ')
        modeled = set().union(*map(set, mappings[arm].values()))
        require(keys.keys() == {key for key in deposit if key[0] in modeled}
                and all(np.array_equal(xyz, deposit[key]) for key, xyz in keys.items()),
                f'{arm} template has changed source atom coordinates or coverage')
        key_sets[arm] = set(keys)
        summaries[arm] = {'source_atom_count': len(keys), 'resolved_positions': len(set(remapped.res_id)),
                          'source_xyz_max_difference_A': 0.0}
    require(key_sets['split'] == key_sets['domain-iv'] < key_sets['ecd'], 'Shared target atom keys differ')
    historical_path = Path('results/experiment-20261114/input-specs/her2_target_refold.cif')
    digest = validate_file(historical_path, manifest['sources'][str(historical_path)]['sha256'])
    historical_cif = CIFFile.read(historical_path)
    historical = get_structure(historical_cif, model=1, use_author_fields=False)
    historical = canonical_target(historical, {'C': 'C', 'D': 'D'}, mappings['split'])
    require(np.array_equal(historical.res_id - 22,
                           historical_cif.block['atom_site']['auth_seq_id'].as_array(int)),
            'Historical template author and canonical positions differ')
    legacy = coordinate_keys(historical)
    common = sorted(key_sets['split'] & legacy.keys())
    moving, fixed = (np.array([source[key] for key in common], dtype=float) for source in (legacy, deposit))
    rotation, translation = target_transform(moving, fixed)
    aligned = rmsd(moving @ rotation + translation, fixed)
    require(len(common) == manifest['audit']['historical_crop_common_atoms']
            and np.isclose(aligned, manifest['audit']['historical_crop_aligned_rmsd_angstrom'], atol=1e-5),
            'Historical target template differs from coordinate audit')
    return {
        'arms': summaries, 'split_domain_iv_shared_atoms': len(key_sets['split']),
        'ecd_additional_atoms': len(key_sets['ecd'] - key_sets['split']),
        'historical_split_sha256': digest, 'historical_split_common_atoms': len(common),
        'historical_split_atom_rmsd_after_rigid_fit_A': aligned,
    }


def read_typed_cbor(stream, major):
    """Read one definite-length item header in the Foundry detailed score schema."""
    head = stream.read(1)
    require(bool(head) and head[0] >> 5 == major, f'Expected CBOR type {major}')
    count = head[0] & 31
    if count < 24:
        return count
    require(count <= 27, 'Unsupported indefinite/reserved CBOR item length')
    width = 1 << (count - 24)
    encoded = stream.read(width)
    require(len(encoded) == width, 'Truncated CBOR item length')
    return int.from_bytes(encoded, 'big')


def cbor_array(stream, tag, dtype, shape):
    require(read_typed_cbor(stream, 6) == tag, 'Unexpected detailed score typed-array tag')
    size = np.dtype(dtype).itemsize * int(np.prod(shape))
    require(read_typed_cbor(stream, 2) == size, 'Unexpected detailed score byte length')
    payload = stream.read(size)
    require(len(payload) == size, 'Truncated detailed score payload')
    return np.frombuffer(payload, dtype=dtype).reshape(shape)


def protenix_pae(path, atoms, n_tokens):
    # Foundry's CBOR uses RFC 8746 little-endian typed arrays and row-major matrices.
    fields = {
        'atom_predicted_local_distance_difference_test', 'token_token_predicted_aligned_error',
        'token_token_predicted_distance_error', 'atom_tokens',
    }
    with gzip.open(path, 'rb') as stream:
        require(read_typed_cbor(stream, 5) == len(fields), f'{path}: unexpected CBOR field count')
        matrices = {}
        for _ in range(len(fields)):
            name = stream.read(read_typed_cbor(stream, 3)).decode()
            require(name in fields, f'{path}: duplicate or unexpected detailed score field {name}')
            fields.remove(name)
            if name.startswith('token_token_'):
                require(read_typed_cbor(stream, 6) == 40 and read_typed_cbor(stream, 4) == 2
                        and read_typed_cbor(stream, 4) == 2
                        and tuple(read_typed_cbor(stream, 0) for _ in range(2)) == (n_tokens, n_tokens),
                        f'{path}: wrong matrix shape')
                matrices[name] = cbor_array(stream, 85, '<f4', (n_tokens, n_tokens))
            elif name == 'atom_tokens':
                indices = cbor_array(stream, 79, '<u8', (len(atoms),))
            else:
                lddt = cbor_array(stream, 85, '<f4', (len(atoms),))
        require(not stream.read(1), f'{path}: trailing data after CBOR map')
    require(np.allclose(lddt, atoms.b_factor, rtol=0, atol=0.0001),
            f'{path}: detailed atom confidences do not match CIF atom order')
    lookup = {key: i for i, key in enumerate(dict.fromkeys(zip(atoms.chain_id, atoms.res_id)))}
    expected = [lookup[key] for key in zip(atoms.chain_id, atoms.res_id)]
    require(np.array_equal(indices, expected), f'{path}: token-to-residue mapping differs from CIF')
    require(all(np.isfinite(array).all() and (array >= 0).all() for array in matrices.values()),
            f'{path}: invalid detailed confidence matrix')
    return matrices['token_token_predicted_aligned_error']


def rf3_summary_order(roles, token_chains, arrays, summary, path):
    # The sidecar documents atom/CIF token order, whereas the summary uses its own
    # order. Match every exposed pair mean and minimum back to sequence-derived roles.
    pairs = list(itertools.combinations(range(len(roles)), 2))
    checks = {}
    for name, matrix in arrays.items():
        for c1, c2 in itertools.permutations(roles, 2):
            block = matrix[np.ix_(token_chains == c1, token_chains == c2)]
            reverse = matrix[np.ix_(token_chains == c2, token_chains == c1)]
            checks[name, c1, c2] = float(block.min()), float((block.mean() + reverse.mean()) / 2)
        for key in (f'chain_pair_{name}', f'chain_pair_{name}_min'):
            rows = summary[key]
            require(len(rows) == len(roles) and all(len(row) == len(roles) for row in rows)
                    and all((value is None) == (i >= j)
                            for i, row in enumerate(rows) for j, value in enumerate(row)),
                    f'{path}: invalid RF3 summary layout for {key}')
    matches = []
    for order in itertools.permutations(roles):
        if all(abs(summary[f'chain_pair_{name}_min'][i][j] - checks[name, order[i], order[j]][0]) <= 0.00501
               and abs(summary[f'chain_pair_{name}'][i][j] - checks[name, order[i], order[j]][1]) <= 0.00501
               for name in arrays for i, j in pairs):
            matches.append(order)
    require(len(matches) == 1, f'{path}: RF3 indexed summaries have no unique chain order')
    return [roles[chain] for chain in matches[0]]


def rf3_confidence(name, atoms, roles, provider):
    folder = CAMPAIGN / RF3
    sidecar = read_json(folder / f'{name}_pae.json')
    summary = read_json(folder / f'{name}_summary_confidences.json')
    token_chains = atoms.chain_id[atoms.atom_name == 'CA']
    count = len(token_chains)
    require(sidecar['model'] == 'rf3' and sidecar['n_tokens'] == count
            and sidecar['token_order'] == 'tokens of the written structure, in file order',
            f'{name}: RF3 PAE token provenance mismatch')
    with np.load(folder / f'{name}_pae.npz', allow_pickle=False) as archive:
        exact_coverage(archive.files, sidecar['arrays'], f'{name} RF3 array fields')
        arrays = {}
        for key in ('pae', 'pde'):
            values, meta = archive[key], sidecar['arrays'][key]
            require(values.shape == (count, count) and list(values.shape) == meta['shape']
                    and str(values.dtype) == meta['dtype'] and meta['units'] == 'angstrom'
                    and np.isfinite(values).all() and (values >= 0).all(),
                    f'{name}: invalid RF3 {key} array')
            arrays[key] = values
    order = rf3_summary_order(roles, token_chains, arrays, summary, name)
    indices = list(map(str, range(len(order))))
    pairs = provider['pair_chains_iptm']
    exact_coverage(pairs, indices, f'{name} RF3 pair ipTM rows')
    for index in indices:
        exact_coverage(pairs[index], indices, f'{name} RF3 pair ipTM columns')
    matrix = np.array([[pairs[i][j] for j in indices] for i in indices])
    require(np.allclose(matrix, matrix.T)
            and np.allclose(matrix.diagonal(), [provider['chains_ptm'][i] for i in indices])
            and np.isclose(summary['iptm'], provider['iptm'], atol=0.000051)
            and np.isclose(summary['ptm'], provider['ptm'], atol=0.000051)
            and np.isclose(summary['overall_plddt'], provider['plddt'], atol=0.000051)
            and summary['has_clash'] == provider['has_clash']
            and (provider['n_tokens'], provider['n_chains'], provider['n_atoms']) == (count, len(order), len(atoms)),
            f'{name}: RF3 indexed confidence or receipt mismatch')
    pair_errors = {
        f'{order[i]}-{order[j]}': {f'{label}_A': summary[f'chain_pair_{label}'][i][j]
                                  for label in ('pae', 'pae_min', 'pde', 'pde_min')}
        for i, j in itertools.combinations(range(len(order)), 2)
    }
    return order, matrix, summary['iptm'], arrays['pae'], summary['has_clash'], pair_errors


def pose_geometry(atoms, roles, refs, sequence_mapping):
    _, ref_ca, ref_ab, ab_sites, footprint, shared, _ = refs
    resolved = {p: xyz for chain, role in roles.items() if role in ('C', 'D', 'T')
                for p, xyz in ca_coordinates(atoms[atoms.chain_id == chain]).items()}
    moving, fixed = (np.array([ca[p] for p in shared], dtype=float) for ca in (resolved, ref_ca))
    rotation, shift = target_transform(moving, fixed)
    fit = {'target_fit_ca_count': len(shared), 'target_fit_ca_rmsd_A': rmsd(moving @ rotation + shift, fixed)}
    if not sequence_mapping:
        return fit
    by_role = {role: ca_coordinates(atoms[atoms.chain_id == chain]) for chain, role in roles.items()
               if role in REFERENCE_CHAINS}
    ref_by_role = {role: ca_coordinates(ab) for role, ab in ref_ab.items()}
    keys = [(role, pos) for role, ca in ref_by_role.items() for pos in ca]
    interface_keys = [(role, pos) for role, sites in ab_sites.items() for pos in sites]
    antibody, reference = (np.array([ca[role][pos] for role, pos in keys], dtype=float)
                           for ca in (by_role, ref_by_role))
    index = [keys.index(key) for key in interface_keys]
    antibody = antibody @ rotation + shift
    antigen_sites = [shared.index(pos) for pos in sorted(footprint)]
    interface, ref_interface = (
        np.vstack([ab[index], ag[antigen_sites]]) for ab, ag in
        [(antibody, moving @ rotation + shift), (reference, fixed)]
    )
    return {
        **fit, 'antibody_mapped_ca_count': len(keys),
        'antibody_ca_rmsd_in_target_frame_A': rmsd(antibody, reference),
        'antibody_centroid_displacement_A': float(np.linalg.norm(antibody.mean(axis=0) - reference.mean(axis=0))),
        'antibody_ref_interface_ca_count': len(index),
        'antibody_ref_interface_ca_rmsd_in_target_frame_A': rmsd(antibody[index], reference[index]),
        'ref_interface_ca_count': len(interface),
        'ref_interface_ca_rmsd_in_target_frame_A': rmsd(interface, ref_interface),
    }


def footprint_pae(pae, canonical_atoms, roles, target_sites, ab_sites):
    # PAE is directional. Keep row/column directions explicit for both exporters.
    tokens = canonical_atoms[canonical_atoms.atom_name == 'CA']
    target = np.flatnonzero(np.isin(tokens.chain_id, [c for c, r in roles.items() if r not in ab_sites])
                            & np.isin(tokens.res_id, sorted(target_sites)))
    antibody = np.flatnonzero(np.logical_or.reduce([
        (tokens.chain_id == chain) & np.isin(tokens.res_id, ab_sites[role])
        for chain, role in roles.items() if role in ab_sites]))
    require(len(target) == len(target_sites) and len(antibody) == sum(map(len, ab_sites.values())),
            'Reference interface PAE residue mapping incomplete')
    return {'ref_interface_pae_ab_rows_target_cols_A': float(pae[np.ix_(antibody, target)].mean()),
            'ref_interface_pae_target_rows_ab_cols_A': float(pae[np.ix_(target, antibody)].mean())}


def antibody_mapping(name, panel):
    wt = panel['control_wt']
    row = panel[name]
    if 'cdrh3' not in row:
        return False
    wt_loop, loop = wt['cdrh3'], row['cdrh3']
    require(wt['vh_sequence'].count(wt_loop) == 1 and len(wt_loop) == len(loop)
            and wt['vh_sequence'].replace(wt_loop, loop) == row['vh_sequence']
            and wt['vl_sequence'] == row['vl_sequence'],
            f'{name} lacks a substitution-only WT antibody residue mapping')
    return True


def geometry_row(atoms, roles, refs, mapping, matched):
    _, _, _, ab_sites, footprint, shared, _ = refs
    mapped = canonical_target(atoms, roles, mapping)
    antibody = [chain for chain, role in roles.items() if role in ab_sites]
    target = [chain for chain, role in roles.items() if role in mapping]
    geometry = contacts(mapped, antibody, target)
    positions = {int(value.split(':')[1]) for value in geometry.pop('target_contact_residues')}
    geometry['antibody_contact_residues'] = sorted(
        (f'{roles[c]}:{p}' for c, p in (v.split(':') for v in geometry['antibody_contact_residues'])),
        key=lambda v: (v[0], int(v.split(':')[1])))
    common = positions & set(shared)
    recovered = common & footprint
    return mapped, {
        'target_contacts_p04626': sorted(positions), 'target_contacts_outside_common_p04626': sorted(positions - common),
        'wt_positions_recovered_p04626': sorted(recovered), 'wt_contact_count': len(recovered),
        'wt_footprint_recall': len(recovered) / len(footprint),
        'wt_footprint_jaccard': len(recovered) / len(common | footprint), **geometry,
        **pose_geometry(mapped, roles, refs, matched),
    }


def score_pairs(order, matrix):
    require(matrix.shape == (len(order), len(order)) and np.isfinite(matrix).all()
            and (matrix >= 0).all() and (matrix <= 1).all(),
            'Invalid chain pair ipTM')
    indices = {role: i for i, role in enumerate(order)}
    def pair(left, right):
        a, b = indices[left], indices[right]
        require(np.isclose(matrix[a, b], matrix[b, a], atol=1e-6), 'Asymmetric chain pair ipTM')
        return float(matrix[a, b])

    pairs = {f'iptm_{a}_{t}': pair(a, t) for t in order if t not in ('H', 'L') for a in ('H', 'L')}
    return {'vh_vl_iptm': pair('H', 'L'), 'max_antibody_target_pair_iptm': max(pairs.values()), **pairs}


def primary_criteria(rows):
    by_key = {(row['model'], row['target_context'], row['query_id']): row for row in rows}
    wt_split = by_key['protenix-v2', 'split', 'control_wt']
    wt, strong, weak, facs = [by_key['protenix-v2', 'domain-iv', name] for name in PANEL[:4]]
    checks = {
        'wt_contact_count_increases_from_split': wt['wt_contact_count'] > wt_split['wt_contact_count'],
        'wt_recovers_at_least_13_of_19': wt['wt_contact_count'] >= 13,
        'wt_has_no_sub_1p5A_heavy_atom_clashes': wt['heavy_atom_clash_pairs_below_1p5A'] == 0,
        'wt_recovers_more_contacts_than_facs_low': wt['wt_contact_count'] > facs['wt_contact_count'],
        'spr_strong_recovers_more_contacts_than_spr_weak': strong['wt_contact_count'] > weak['wt_contact_count'],
    }
    return {'criterion_met': all(checks.values()), 'checks': checks}


def audit(reference):
    archive = validate_archive()
    manifest, panel, canonical, queries, mappings, staged = validate_inputs()
    image = (CAMPAIGN / PROTENIX / 'image-id.txt').read_text().strip()
    require(image == IMAGE_ID, 'Unexpected Protenix image ID')
    protenix = {arm: validate_protenix_run(arm, queries[arm], manifest['arms'][arm]) for arm in ARMS}
    rf3_rows, rf3_run = validate_rf3_run(queries['domain-iv'])
    models = {}
    for arm in ARMS:
        for name in PANEL:
            path = protenix[arm][0][name]['path']
            models[arm, name] = validate_prediction(path, queries[arm][name], protenix=True)
    for name in PANEL:
        path = CAMPAIGN / RF3 / f'{name}.cif'
        models['rf3', name] = validate_prediction(path, queries['domain-iv'][name], protenix=False)

    refs = read_reference(reference, manifest, canonical, queries['split']['control_wt'], mappings)
    templates = audit_templates(manifest, mappings, queries, refs[0])
    rows = []
    for name in PANEL:
        for arm in [*ARMS, 'rf3']:
            atoms, roles = models[arm, name]
            model = 'rf3' if arm == 'rf3' else 'protenix-v2'
            context = 'domain-iv' if arm == 'rf3' else arm
            if arm == 'rf3':
                order, matrix, overall, pae, clash, pair_errors = rf3_confidence(name, atoms, roles, rf3_rows[name])
            else:
                files = protenix[arm][0][name]
                score = read_json(files['scores_path'])
                n_tokens = int(np.count_nonzero(atoms.atom_name == 'CA'))
                require((score['num_chains'], score['num_atoms'], score['num_tokens'], score['num_residues'])
                        == (len(roles), len(atoms), n_tokens, n_tokens),
                        f'{arm}/{name} Protenix structure and confidence counts differ')
                order = list(roles.values())
                matrix = np.asarray(score['chain_chain_interface_predicted_tm_score'], dtype=float)
                overall = float(score['interface_predicted_tm_score'])
                pae = protenix_pae(files['detailed_scores_path'], atoms, n_tokens)
                clash = score['clash_flag']
            matched = antibody_mapping(name, panel)
            mapped, geom = geometry_row(atoms, roles, refs, mappings[context], matched)
            require(0 <= overall <= 1 and isinstance(clash, bool), f'{arm}/{name}: invalid global confidence')
            rows.append({
                'query_id': name, 'model': model, 'target_context': context, 'seed': SEED, 'sample': SAMPLE,
                'cdrh3': panel[name].get('cdrh3'), 'spr_kd_nM': panel[name].get('kd_nM'),
                'facs_class': panel[name].get('facs_class'), 'cif_chain_roles': roles,
                'confidence_chain_order': order, 'overall_iptm': overall, 'provider_clash_flag': clash,
                'antibody_residue_mapping': '1-based variable positions, equal-length CDRH3 substitutions'
                if matched else None, **score_pairs(order, matrix), **geom,
                **({'rf3_chain_pair_errors': pair_errors} if arm == 'rf3' else {}),
                **(footprint_pae(pae, mapped, roles, refs[4], refs[3]) if name.startswith('control_') else {}),
            })
    provenance = {
        'coverage': {'protenix': 21, 'rf3': 7, 'inference_failures': []},
        'prespecified_protenix_domain_iv_rescue': primary_criteria(rows),
        'inputs_sha256': staged, 'archive': archive, 'protenix_image_id': image,
        'protenix': {arm: info for arm, (_, info) in protenix.items()}, 'rf3': rf3_run,
        'reference': {
            'source': manifest['sources']['pdb_cif']['url'], 'sha256': refs[-1],
            'target_chain': 'C', 'target_canonical_offset': 22, 'antibody_role_to_chain': REFERENCE_CHAINS,
            'target_ca_fit_positions_p04626': refs[5], 'wt_footprint_positions_p04626': sorted(refs[4]),
            'antibody_interface_positions': refs[3],
        },
        'source_templates': templates,
        'limits': ['Protenix logs give aggregate assignments, with target-only sequences and exact target coordinates '
                   'verified from the staged template files; per-query applied template masks were not exported.',
                   'RF3 model weights and submitted query payload were not archived by the provider.',
                   'RF3 indexed pair ipTM uses the inferred chain-summary order. The receipt has no explicit index map.',
                   'One sample per model/context and query, with no measurements made by these runs.',
                   'De novo antibody RMSD and interface PAE are omitted; no WT-equivalent antibody residue map.'],
    }
    return rows, provenance


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--reference', type=Path, required=True, help='Original 1N8Z mmCIF pinned in the input manifest')
    parser.add_argument('--output', type=Path, default=Path('results/context-panel'))
    args = parser.parse_args()
    rows, provenance = audit(args.reference)
    args.output.mkdir(parents=True, exist_ok=True)
    (args.output / 'audit.json').write_text(json.dumps({'provenance': provenance, 'predictions': rows},
                                                      indent=2, allow_nan=False) + '\n')
    # The JSON retains variable-length contact lists and all chain identities.
    fields = list(dict.fromkeys(k for row in rows for k, v in row.items() if not isinstance(v, (dict, list))))
    with (args.output / 'scores.csv').open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fields, extrasaction='ignore', lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)
    print(f'Audited {len(rows)} predictions with no inference failures.')


if __name__ == '__main__':
    main()
