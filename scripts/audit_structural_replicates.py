"""Audit the frozen, local OpenDDE HER2 domain IV structural replicates."""
from __future__ import annotations

import argparse
import csv
from itertools import permutations
import json
from pathlib import Path
from zipfile import BadZipFile, ZipFile

import numpy as np
from biotite.sequence import ProteinSequence

from scripts.audit_context_panel import (
    ARMS, CONFIG, PANEL, REFERENCE_CHAINS, antibody_mapping, exact_coverage,
    file_hash, footprint_pae, geometry_row, read_json, read_reference,
    rf3_summary_order, score_pairs, validate_file, validate_prediction,
)
from scripts.prepare_her2_context import canonical_sequence, require, sha256
from src.her2_structure_audit import contacts

MODEL = 'opendde-abag'
SEEDS = (42, 43, 44)
CAMPAIGN = Path('results/campaign/opendde-replicates')
RECEIPTS = Path('results/provider-jobs/opendde-replicates')
OUTPUT = Path('results/structural-replicates')
SUFFIXES = ('.cif', '_pae.npz', '_pae.json')
TOKEN_ORDER = 'tokens of the written structure, in file order'
WT_CUTOFF = 13
RANGE_METRICS = (
    'wt_contact_count', 'wt_footprint_recall', 'wt_footprint_jaccard',
    'target_fit_ca_rmsd_A', 'antibody_ca_rmsd_in_target_frame_A',
    'antibody_ref_interface_ca_rmsd_in_target_frame_A', 'ref_interface_ca_rmsd_in_target_frame_A',
    'overall_iptm', 'vh_vl_iptm', 'iptm_H_T', 'iptm_L_T',
    'heavy_atom_clash_pairs_below_1p5A', 'vh_vl_heavy_atom_clash_pairs_below_1p5A',
    'pae_mean_H_L_A', 'pae_mean_L_H_A', 'pae_mean_H_T_A', 'pae_mean_T_H_A',
    'pae_mean_L_T_A', 'pae_mean_T_L_A',
    'ref_interface_pae_ab_rows_target_cols_A', 'ref_interface_pae_target_rows_ab_cols_A',
)


def frozen_maps(manifest, canonical):
    mappings = {}
    for arm in ('split', 'domain-iv'):
        chains = manifest['arms'][arm]['target_chains']
        exact_coverage(chains, ARMS[arm], f'{arm} target mapping')
        mappings[arm] = {}
        for chain, details in chains.items():
            positions = details['query_position_to_canonical']
            require(positions == list(ARMS[arm][chain]) and len(positions) == details['length'],
                    f'{arm}/{chain}: unexpected canonical positions')
            sequence = ''.join(canonical[p - 1] for p in positions)
            require(sha256(sequence.encode()) == details['sequence_sha256'],
                    f'{arm}/{chain}: unexpected target sequence digest')
            mappings[arm][chain] = positions
    require(manifest['sifts']['mature_to_canonical_offset'] == 22, 'Unexpected 1N8Z SIFTS offset')
    return mappings


def frozen_queries(config=CONFIG):
    manifest_path = config / 'manifest.json'
    manifest = read_json(manifest_path)
    paths = ('P04626.fasta', 'domain-iv-queries.json')
    hashes = {name: validate_file(config / name, manifest['artifact_sha256'][name]) for name in paths}
    canonical = canonical_sequence((config / paths[0]).read_bytes())
    require(sha256(canonical.encode()) == manifest['canonical_sequence_sha256'],
            'Unexpected canonical sequence digest')
    mappings = frozen_maps(manifest, canonical)
    expected_target = ''.join(canonical[p - 1] for p in mappings['domain-iv']['T'])
    panel = {row['query_id']: row for row in manifest['panel']}
    items = read_json(config / paths[1])['queries']
    exact_coverage([row['query_id'] for row in manifest['panel']], PANEL, 'Manifest panel')
    exact_coverage([row['query_id'] for row in items], PANEL, 'Frozen query IDs')
    queries = {}
    for row in items:
        name = row['query_id']
        chains = row['chains']
        exact_coverage([chain['chain_id'] for chain in chains], ('T', 'H', 'L'), name, ordered=True)
        expected = {'T': expected_target}
        for role, kind in (('H', 'vh'), ('L', 'vl')):
            expected[role] = panel[name][f'{kind}_sequence']
            require(sha256(expected[role].encode()) == panel[name][f'{kind}_sha256'],
                    f'{name}: unexpected {kind} sequence digest')
        require(all(c['sequence'] == expected[c['chain_id']] and c['msa_path'] is None for c in chains),
                f'{name}: frozen sequence or MSA mismatch')
        queries[name] = expected
    return manifest, panel, canonical, queries, mappings, {**hashes, 'manifest.json': file_hash(manifest_path)}


def expected_params(seed):
    return {'seed': seed, 'diffusion_samples': 1, 'use_msa_server': False,
            'write_pae': True, 'output_format': 'cif'}


def validate_request(request, seed, queries):
    # ASVS V2.2: derive the admissible request from the frozen sequences and protocol.
    name = f'her2-domain-iv-opendde-seed{seed}'
    require(set(request) == {'model', 'name', 'params', 'targets'} and request['model'] == MODEL
            and request['name'] == name
            and json.dumps(request['params'], sort_keys=True) == json.dumps(expected_params(seed), sort_keys=True),
            f'Seed {seed}: request model, name or protocol mismatch')
    targets = request['targets']
    exact_coverage([item['name'] for item in targets], queries, f'Seed {seed} request targets')
    for item in targets:
        query = item['name']
        fasta = [f'>{role}\n{seq}' for role, seq in queries[query].items()]
        require(set(item) == {'name', 'content'} and item['content'].splitlines() == '\n'.join(fasta).splitlines(),
                f'Seed {seed}/{query}: expected exact T/H/L FASTA and no per-target parameters')
    return name


def validate_job(job, seed, name, queries):
    accepted = {**expected_params(seed), 'fast': False, 'use_tfg_guidance': False}
    require(job['model'] == MODEL and job['kind'] == 'predict' and job['name'] == name
            and isinstance(job['id'], str) and bool(job['id']) and job['total'] == len(queries)
            and json.dumps(job['params'], sort_keys=True) == json.dumps(accepted, sort_keys=True),
            f'Seed {seed}: accepted job settings mismatch')
    require(job['error'] is None and job['status'] in {'queued', 'running', 'succeeded'},
            f'Seed {seed}: job reports failure')
    if job['status'] == 'succeeded':
        require(job['results_ready'] is True and job['done'] == len(queries),
                f'Seed {seed}: incomplete finished job')
    if 'targets' in job:
        exact_coverage([target['id'] for target in job['targets']], queries, f'Seed {seed} job targets')
        if job['status'] == 'succeeded':
            require(all(target['state'] == 'done' for target in job['targets']),
                    f'Seed {seed}: incomplete target in finished job')


def submission(seed, queries, receipts):
    request_path = receipts / f'request-seed{seed}.json'
    job_path = receipts / f'job-seed{seed}.json'
    name = validate_request(read_json(request_path), seed, queries)
    job = read_json(job_path)
    validate_job(job, seed, name, queries)
    return job, {path.name: file_hash(path) for path in (request_path, job_path)}


def verify_result_receipt(path, job, rows, queries):
    receipt = read_json(path)
    require(receipt['job_id'] == job['id'] and receipt['ready'] is True
            and receipt['kind'] == 'predict' and receipt['rows'] == rows,
            f'{path}: results receipt differs from job or extracted rows')
    require(receipt['structures'] == {name: [f'{name}.cif'] for name in queries},
            f'{path}: unexpected structure paths')
    expected = {f'{name}{suffix}': kind for name in queries for suffix, kind in
                zip(SUFFIXES, ('structure', 'pae', 'pae_sidecar'))}
    artifacts = {}
    duplicate_sidecars = set()
    for item in receipt['artifacts']:
        member = item['path']
        require(member in expected and item['type'] == expected[member], f'{path}: unexpected artifact')
        if member in artifacts:
            require(item == artifacts[member] and item['type'] == 'pae_sidecar',
                    f'{path}: incompatible duplicate artifact')
            duplicate_sidecars.add(member)
        if item['type'] == 'structure':
            require(item['target'] == Path(member).stem, f'{path}: structure target mismatch')
        artifacts[member] = item
    exact_coverage(artifacts, expected, f'{path} artifact list')
    for kind in ('pae', 'pae_sidecar'):
        paths = receipt[kind]
        if kind == 'pae_sidecar':
            paths = list(dict.fromkeys(paths))
        exact_coverage(paths, [p for p, v in expected.items() if v == kind], f'{path} {kind} list')
    return sorted(duplicate_sidecars)


def validate_rows(rows, queries, seed):
    exact_coverage([row['id'] for row in rows], queries, f'Seed {seed} results')
    require(all(row['status'] == 'ok' and type(row['samples']) is int and row['samples'] == 1
                and row['msa'] is False and row['msa_depth'] == 1
                and row.get('seed', seed) == seed for row in rows),
            f'Seed {seed}: failed query or sample, seed or MSA mismatch')


def local_artifacts(folder, queries):
    required = {'results.json', *(f'{name}{suffix}' for name in queries for suffix in SUFFIXES)}
    require(folder.is_dir(), f'Awaiting extracted artifacts in {folder}')
    entries = {path.name: path for path in folder.iterdir()}
    archives = [entries.pop(name) for name in ('provider-results.zip', 'results.zip') if name in entries]
    require(len(archives) <= 1, f'{folder}: ambiguous provider archive')
    summaries = {f'{name}_summary_confidences.json' for name in queries}
    if entries.keys() & summaries:
        required |= summaries
    exact_coverage(entries, required, f'{folder} extracted artifacts')
    require(all(path.is_file() and not path.is_symlink() for path in entries.values()),
            f'{folder}: expected ordinary artifact files')
    # ASVS V11.4: hash read-only inputs with SHA-256 and compare the archive when available.
    hashes = {name: file_hash(path) for name, path in sorted(entries.items())}
    if archives:
        archive = archives[0]
        verify_archive(archive, entries, required)
        hashes[archive.name] = file_hash(archive)
    return hashes


def verify_archive(path, extracted, required):
    require(path.is_file() and not path.is_symlink(), f'{path}: expected ordinary archive file')
    with ZipFile(path) as zipped:
        exact_coverage(zipped.namelist(), required, f'{path} members')
        for member in zipped.infolist():
            local = extracted[member.filename]
            require(member.file_size == local.stat().st_size, f'{local}: archive member size differs')
            with zipped.open(member) as stream, local.open('rb') as contents:
                while data := stream.read(65536):
                    require(data == contents.read(len(data)), f'{local}: archive member differs')


def verify_download_hashes(path, seed, job_id, folder, hashes):
    receipt = read_json(path)
    require(receipt['seed'] == seed and receipt['job_id'] == job_id, f'{path}: seed or job mismatch')
    files = receipt['files']
    expected = {f'seed-{seed}/{name}': name for name in hashes}
    exact_coverage([item['path'] for item in files], expected, f'{path} hashed files')
    for item in files:
        name = expected[item['path']]
        require(item['sha256'] == hashes[name] and item['bytes'] == (folder / name).stat().st_size,
                f'{folder / name}: download receipt SHA-256 or size mismatch')


def completed_seed(seed, job, queries, folder, receipts):
    hashes = local_artifacts(folder, queries)
    rows = read_json(folder / 'results.json')
    validate_rows(rows, queries, seed)
    completed_path = receipts / f'completed-seed{seed}.json'
    complete = read_json(completed_path) if completed_path.is_file() else job
    validate_job(complete, seed, job['name'], queries)
    require(complete['status'] == 'succeeded' and complete['id'] == job['id'],
            f'Seed {seed}: awaiting matching completed job receipt')
    result_path = receipts / f'results-seed{seed}.json'
    file_path = receipts / f'files-seed{seed}.json'
    duplicate_sidecars = verify_result_receipt(result_path, job, rows, queries)
    verify_download_hashes(file_path, seed, job['id'], folder, hashes)
    receipt_paths = [path for path in (completed_path, result_path, file_path) if path.is_file()]
    return {row['id']: row for row in rows}, hashes, receipt_paths, duplicate_sidecars


def token_axes(atoms, roles, expected, sidecar, path):
    # These are the unique protein residue tokens in atom-file order. Check every
    # residue against its independently identified chain role and submitted sequence.
    ca = atoms[atoms.atom_name == 'CA']
    expected_count = sum(map(len, expected.values()))
    ordered_keys = list(dict.fromkeys(zip(atoms.chain_id, atoms.res_id)))
    require(sidecar['model'] == MODEL and sidecar['n_tokens'] == len(ca) == expected_count
            and sidecar['token_order'] == TOKEN_ORDER and list(zip(ca.chain_id, ca.res_id)) == ordered_keys
            and all(1 <= pos <= len(expected[roles[chain]])
                    and ProteinSequence.convert_letter_3to1(res) == expected[roles[chain]][pos - 1]
                    for chain, pos, res in zip(ca.chain_id, ca.res_id, ca.res_name)),
            f'{path}: PAE token order or full CIF sequence mismatch')
    return ca.chain_id


def error_arrays(path, sidecar, count):
    # ASVS V5.2: bound compressed member count/size before NumPy loads numeric data.
    fields = sidecar['arrays']
    require('pae' in fields and len(fields) <= 16, f'{path}: missing PAE or too many array fields')
    with ZipFile(path) as archive:
        exact_coverage(archive.namelist(), [f'{key}.npy' for key in fields], f'{path} members')
        require(all(item.file_size <= count ** 2 * 8 + 4096 for item in archive.infolist()),
                f'{path}: oversized confidence array member')
    matrices = {}
    with np.load(path, allow_pickle=False) as archive:
        for key in ('pae', 'pde'):
            if key not in fields:
                continue
            values, meta = archive[key], fields[key]
            require(values.shape == (count, count) and list(values.shape) == meta['shape']
                    and str(values.dtype) == meta['dtype'] and values.dtype.kind == 'f'
                    and meta['units'] == 'angstrom'
                    and np.isfinite(values).all() and (values >= 0).all(),
                    f'{path}: invalid {key} distances, dtype or token count')
            matrices[key] = values
    return matrices


def indexed_matrix(values, size, label):
    if isinstance(values, dict):
        indices = [str(i) for i in range(size)]
        exact_coverage(values, indices, f'{label} rows')
        for row in values.values():
            exact_coverage(row, indices, f'{label} columns')
        values = [[values[i][j] for j in indices] for i in indices]
    matrix = np.asarray(values, dtype=float)
    require(matrix.shape == (size, size), f'{label}: invalid pair ipTM shape')
    return matrix


def chain_mean_error_order(roles, tokens, arrays, provider, label):
    """Infer index roles by concordance of chain pTM and within-chain PAE/PDE."""
    size = len(roles)
    indices = [str(i) for i in range(size)]
    ptm = provider['chains_ptm']
    exact_coverage(ptm, indices, f'{label} chain pTM indices')
    scores = np.array([probability_score(ptm[i], f'{label} chain {i} pTM') for i in indices])
    masks = {chain: tokens == chain for chain in roles}
    require(set(arrays) == {'pae', 'pde'}, f'{label}: both PAE and PDE required for rank inference')
    means = {key: {role: float(matrix[np.ix_(masks[chain], masks[chain])].mean(dtype=float))
                   for chain, role in roles.items()} for key, matrix in arrays.items()}
    ranked_roles = [sorted(values, key=values.get) for values in means.values()]
    require(len(set(scores)) == size and all(len(set(values.values())) == size for values in means.values())
            and all(rank == ranked_roles[0] for rank in ranked_roles),
            f'{label}: chain PAE/PDE ranks do not resolve a unique confidence order')
    order = [''] * size
    for index, role in zip(np.argsort(-scores), ranked_roles[0]):
        order[index] = role
    evidence = {
        'method': 'descending indexed chain pTM matched to ascending within-chain mean PAE and PDE',
        'status': 'rank inference; no explicit provider index labels or numerical pair ipTM reconstruction',
        'within_chain_error_means_A': means,
        'indexed_chain_ptm': ptm,
    }
    return order, evidence


def probability_score(value, label):
    require(type(value) in (int, float) and np.isfinite(value) and 0 <= value <= 1,
            f'{label}: invalid probability-scale confidence')
    return float(value)


def indexed_iptm(provider, summary, order, label):
    sources = [(field, values[field]) for values in (provider, summary)
               for field in ('pair_chains_iptm', 'chain_pair_iptm') if field in values]
    if not sources:
        return {key: None for key in ('vh_vl_iptm', 'iptm_H_T', 'iptm_L_T')}, None
    size = len(order)
    field, values = sources[0]
    matrix = indexed_matrix(values, size, label)
    scores = score_pairs(order, matrix)
    for other_field, other in sources[1:]:
        require(np.allclose(matrix, indexed_matrix(other, size, label), rtol=0, atol=0.00501),
                f'{label}: incompatible {field} and {other_field}')
    if 'chains_ptm' in provider:
        entries = provider['chains_ptm']
        indices = [str(i) for i in range(size)]
        exact_coverage(entries, indices, f'{label} chain pTM')
        require(np.allclose(matrix.diagonal(), [entries[i] for i in indices], rtol=0, atol=1e-6),
                f'{label}: pair ipTM diagonal differs from indexed chain pTM')
    return scores, field


def check_confidence_receipt(provider, summary, n_atoms, n_tokens, n_chains, label):
    for key, count in (('n_atoms', n_atoms), ('n_tokens', n_tokens), ('n_chains', n_chains)):
        require(provider[key] == count, f'{label}: incorrect {key}')
    if 'n_residues' in provider:
        require(provider['n_residues'] == n_tokens, f'{label}: incorrect n_residues')
    for key, summary_key in (('iptm', 'iptm'), ('ptm', 'ptm'), ('plddt', 'overall_plddt')):
        available = [probability_score(record[k], f'{label}/{k}')
                     for record, k in ((provider, key), (summary, summary_key)) if k in record]
        require(len(available) < 2 or abs(available[0] - available[1]) <= 0.000051,
                f'{label}: result and confidence summary disagree on {key}')
    for record in (provider, summary):
        if 'has_clash' in record:
            require(type(record['has_clash']) is bool, f'{label}: invalid provider clash flag')
    if 'has_clash' in provider and 'has_clash' in summary:
        require(summary['has_clash'] == provider['has_clash'], f'{label}: provider clash flags disagree')


def pair_pae(pae, tokens, roles):
    return {
        f'pae_mean_{roles[a]}_{roles[b]}_A': float(pae[np.ix_(tokens == a, tokens == b)].mean(dtype=float))
        for a, b in permutations(roles, 2)
    }


def confidence(folder, name, atoms, roles, expected, provider):
    label = folder / name
    sidecar = read_json(folder / f'{name}_pae.json')
    summary_path = folder / f'{name}_summary_confidences.json'
    summary = read_json(summary_path) if summary_path.is_file() else {}
    tokens = token_axes(atoms, roles, expected, sidecar, label)
    arrays = error_arrays(folder / f'{name}_pae.npz', sidecar, len(tokens))
    check_confidence_receipt(provider, summary, len(atoms), len(tokens), len(roles), label)
    if summary_path.is_file():
        order = rf3_summary_order(roles, tokens, arrays, summary, label)
        evidence = {'method': 'match pair summary directional minima and bidirectional means',
                    'matched_arrays': list(arrays)}
    else:
        order, evidence = chain_mean_error_order(roles, tokens, arrays, provider, label)
    scores, iptm_field = indexed_iptm(provider, summary, order, label)
    return arrays['pae'], {
        'pae_token_chain_order': [roles[chain] for chain in dict.fromkeys(tokens)],
        'pae_token_chain_lengths': {roles[chain]: int((tokens == chain).sum()) for chain in roles},
        'pae_token_order_declaration': sidecar['token_order'],
        'indexed_confidence_chain_order_inferred': order,
        'confidence_axis_evidence': evidence,
        'pair_iptm_field': iptm_field,
        'pair_iptm_axis_basis': evidence['method'] if iptm_field else None,
        'provider_indexed_pair_iptm': provider.get('pair_chains_iptm'),
        'overall_iptm': summary.get('iptm', provider.get('iptm')),
        'provider_clash_flag': summary.get('has_clash', provider.get('has_clash')),
        **scores, **pair_pae(arrays['pae'], tokens, roles),
    }


def score_query(seed, name, folder, expected, provider, panel, refs, mapping):
    atoms, roles = validate_prediction(folder / f'{name}.cif', expected, protenix=False)
    pae, scores = confidence(folder, name, atoms, roles, expected, provider)
    matched = antibody_mapping(name, panel)
    mapped, geometry = geometry_row(atoms, roles, refs, mapping, matched)
    chains = {role: chain for chain, role in roles.items()}
    hl = contacts(atoms, [chains['H']], [chains['L']])
    return {
        'seed': seed, 'query_id': name, 'cif_chain_to_role': roles,
        'antibody_residue_mapping': '1-based variable positions with equal-length CDRH3 substitutions'
        if matched else None,
        'wt_contact_criterion_met': geometry['wt_contact_count'] >= WT_CUTOFF,
        'vh_vl_heavy_atom_clash_pairs_below_1p5A': hl['heavy_atom_clash_pairs_below_1p5A'],
        'vh_vl_minimum_heavy_atom_distance_A': hl['minimum_target_distance_A'],
        **scores, **geometry, **(footprint_pae(pae, mapped, roles, refs[4], refs[3]) if matched else {}),
    }


def query_ranges(rows):
    ranges = []
    for name in PANEL:
        group = [row for row in rows if row['query_id'] == name]
        exact_coverage([row['seed'] for row in group], SEEDS, name, ordered=True)
        metrics = {}
        for field in RANGE_METRICS:
            values = [row[field] for row in group if row.get(field) is not None]
            metrics[field] = {'min': min(values), 'max': max(values), 'n': len(values)} if values else None
        ranges.append({'query_id': name, 'seeds': [row['seed'] for row in group], 'metrics': metrics})
    return ranges


def compare_controls(rows):
    by_key = {(row['seed'], row['query_id']): row for row in rows}
    checks = []
    for seed in SEEDS:
        wt, strong, weak, facs = [by_key[seed, name] for name in PANEL[:4]]
        checks.append({
            'seed': seed, 'wt_at_least_13_of_19': wt['wt_contact_criterion_met'],
            'wt_recovered': wt['wt_contact_count'], 'wt_minus_facs_low_contacts':
                wt['wt_contact_count'] - facs['wt_contact_count'],
            'spr_strong_minus_weak_contacts': strong['wt_contact_count'] - weak['wt_contact_count'],
            'wt_has_no_antibody_target_heavy_atom_clashes': wt['heavy_atom_clash_pairs_below_1p5A'] == 0,
        })
    return checks


def seed_provenance(seed, job, hash_receipts, hashes, paths, duplicate_sidecars, receipts, provider):
    hash_receipts.update({path.name: file_hash(path) for path in paths})
    usage_path = receipts / f'usage-seed{seed}.json'
    usage = None
    if usage_path.is_file():
        record = read_json(usage_path)
        require(record['job'] == job['id'] and record['model'] == MODEL and record['status'] == 'succeeded',
                f'Seed {seed}: usage job mismatch')
        exact_coverage([item['id'] for item in record['targets']], provider, f'Seed {seed} usage targets')
        require(all(item['attempts'] == 1 for item in record['targets']), f'Seed {seed}: usage attempt mismatch')
        usage = {key: record[key] for key in ('currency', 'amount', 'billable_s')}
        hash_receipts[usage_path.name] = file_hash(usage_path)
    return {
        'seed': seed, 'job_id': job['id'], 'model': MODEL, 'params': job['params'],
        'status_in_submission_receipt': job['status'], 'completed_job_status': 'succeeded',
        'duplicate_result_receipt_sidecar_paths': duplicate_sidecars,
        'query_rows_with_explicit_seed': sum('seed' in row for row in provider.values()),
        'queries_completed': len(provider), 'receipts_sha256': hash_receipts, 'artifacts_sha256': hashes,
        'usage': usage,
    }


def audit(reference, campaign=CAMPAIGN, receipts=RECEIPTS):
    manifest, panel, canonical, queries, mappings, hashes = frozen_queries()
    jobs = {seed: submission(seed, queries, receipts) for seed in SEEDS}
    require(len({job['id'] for job, _ in jobs.values()}) == len(SEEDS), 'Seed job IDs must be distinct')
    runs, rows = [], []
    refs = read_reference(reference, manifest, canonical, queries['control_wt'], mappings)
    require(len(refs[5]) == 109 and len(refs[4]) == 19, 'Unexpected reference target or footprint size')
    for seed, (job, receipt_hashes) in jobs.items():
        folder = campaign / f'seed-{seed}'
        provider, artifacts, paths, duplicate_sidecars = completed_seed(seed, job, queries, folder, receipts)
        runs.append(seed_provenance(seed, job, receipt_hashes, artifacts, paths, duplicate_sidecars, receipts, provider))
        for name in PANEL:
            rows.append(score_query(seed, name, folder, queries[name], provider[name], panel,
                                    refs, mappings['domain-iv']))
    require(len(rows) == 21, 'Expected 21 validated structures')
    return {
        'predictions': rows, 'per_query_ranges': query_ranges(rows), 'controls_by_seed': compare_controls(rows),
        'provenance': {
            'coverage': {'seeds': list(SEEDS), 'query_ids': list(PANEL), 'predictions': len(rows)},
            'runs': runs, 'inputs_sha256': hashes,
            'request_policy': 'FASTA only; no supplied templates or restraints; accepted TFG guidance false',
            'seed_evidence': 'provider accepted job parameters; per-row seeds checked where present',
            'assay_evidence': {name: {key: panel[name][key] for key in ('kd_nM', 'all_wt_kd_nM', 'facs_class')
                                      if key in panel[name]} for name in PANEL[:4]},
            'model_weights_revision': None,
            'reference': {
                'sha256': refs[6], 'target_chain': 'C', 'target_canonical_offset': 22,
                'antibody_role_to_chain': REFERENCE_CHAINS, 'target_ca_fit_positions_p04626': refs[5],
                'wt_footprint_positions_p04626': sorted(refs[4]), 'antibody_interface_positions': refs[3],
            },
        },
    }


def write_results(output, result):
    rows = result['predictions']
    fields = list(dict.fromkeys(k for row in rows for k, v in row.items() if not isinstance(v, (dict, list))))
    output.mkdir(parents=True, exist_ok=True)
    payload = json.dumps(result, indent=2, allow_nan=False) + '\n'
    (output / 'audit.json').write_text(payload)
    with (output / 'scores.csv').open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fields, extrasaction='ignore', lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--reference', type=Path, required=True, help='Original 1N8Z mmCIF pinned in the manifest')
    parser.add_argument('--campaign', type=Path, default=CAMPAIGN)
    parser.add_argument('--receipts', type=Path, default=RECEIPTS)
    parser.add_argument('--output', type=Path, default=OUTPUT)
    args = parser.parse_args()
    try:
        result = audit(args.reference, args.campaign, args.receipts)
    except (ValueError, OSError, KeyError, TypeError, BadZipFile) as error:
        parser.exit(1, f'Audit stopped: {error}\n')
    write_results(args.output, result)
    print(f'Audited {len(result["predictions"])} OpenDDE predictions across seeds {SEEDS}.')


if __name__ == '__main__':
    main()
