"""Audit every completed pose in the matched HER2 template ablation on CPU."""

import argparse
import csv
import itertools
import json
import re
from decimal import Decimal
from pathlib import Path

import gemmi
import numpy as np
from Bio import SeqIO
from scipy.spatial import cKDTree

from scripts.audit_context_panel import (
    IMAGE_ID, REFERENCE_CHAINS, antibody_mapping, ca_coordinates, exact_coverage,
    file_hash, footprint_pae, geometry_row, protenix_pae, read_reference, rmsd, score_pairs,
    target_transform, validate_file, validate_prediction,
)
from scripts.prepare_her2_template_controls import (
    CONFIG, CONTEXT, MASK, OFFSET, PANEL, POSITIONS, QUERY_FILE, ROOT, SAMPLES, SEEDS,
    SOURCE_SHA256, SOURCE_URL, TEMPLATES, chain_sequences, frozen_context, make_templates,
    require, sha256, source_target,
)

OUT = ROOT / "results/template-controls"
MSA_CONFIG = ROOT / "configs/her2-target-msa"
SCHEMA = "foundry.cli.v1"
SUFFIXES = {
    "path": "_model.cif",
    "scores_path": "_scores.json",
    "detailed_scores_path": "_detailed_scores.cbor.gz",
}
ROLES = ("T", "H", "L")
MAPPING = {"T": list(POSITIONS)}
COMMON = [p for p in POSITIONS if p not in MASK]
EXPECTED = list(itertools.product(TEMPLATES, PANEL, SEEDS, SAMPLES))


def read_json(path):
    return json.loads(path.read_bytes())


def validate_staged(campaign, expected, alignment_path=None):
    staged = {}
    full_paths = {}
    for line in (campaign / "input.sha256").read_text().splitlines():
        digest, path = line.split(maxsplit=1)
        name = Path(path).name
        require(name not in staged, f"Duplicate staged input {name}")
        staged[name] = digest
        full_paths[name] = path
    require(staged == expected, "Staged input hashes differ from verified local inputs")
    require(alignment_path is None or full_paths.get("target.a3m") == alignment_path,
            "Query target MSA path differs from staged alignment")
    image = (campaign / "image-id.txt").read_text().strip()
    require(image == IMAGE_ID, "Unexpected local Protenix v2 image ID")
    return staged, image, full_paths


def validate_inputs(config, campaign, template_source, msa_config=None):
    canonical, panel, expected_queries, sources = frozen_context()
    manifest_path = config / "manifest.json"
    manifest = read_json(manifest_path)
    validated_6oge = validate_file(template_source, SOURCE_SHA256)
    require(manifest["sources"] == {SOURCE_URL: validated_6oge, **sources}
            and manifest["panel"] == panel
            and manifest["reference"]["common_canonical_positions"] == COMMON
            and manifest["run"] == {"model": "protenix-v2", "msa": False, "seeds": list(SEEDS),
                                    "samples": list(SAMPLES), "expected_poses": len(EXPECTED)},
            "Frozen experiment manifest mismatch")
    exact_coverage(manifest["artifact_sha256"], [QUERY_FILE, *TEMPLATES.values()], "Input artifact manifest")
    expected_staged = {name: validate_file(config / name, digest)
                       for name, digest in manifest["artifact_sha256"].items()}
    expected_staged["manifest.json"] = file_hash(manifest_path)

    require(read_json(config / QUERY_FILE) == expected_queries, "Query values differ from frozen controls")
    alignment_path = None
    if msa_config is not None:
        expected_staged, alignment_path = validate_target_msa(
            msa_config, expected_queries, expected_staged[TEMPLATES["complete"]])
    staged, image, staged_paths = validate_staged(campaign, expected_staged, alignment_path)

    deposit = gemmi.cif.read(str(template_source)).sole_block()
    templates, arms, bonds, crossing = make_templates(deposit, canonical)
    deposited_mapping = [
        {"query_index": i, "template_label_seq_id": i, "deposited_label_seq_id": p - OFFSET,
         "deposited_auth_seq_id": p, "canonical_position": p, "aa": canonical[p - 1]}
        for i, p in enumerate(POSITIONS, 1)
    ]
    expected_target = {
        "accession": "P04626", "query_chain": "T", "deposited_label_chain": "A",
        "deposited_author_chain": "A", "entity_id": "1", "positions": deposited_mapping,
        "internal_disulfides_p04626": bonds, "disulfides_crossing_crop_boundary_p04626": crossing,
    }
    require(manifest["arms"] == arms and manifest["target"] == expected_target,
            "Template coverage, residue mapping or disulfides disagree with 6OGE")
    for name, text in templates.items():
        actual = gemmi.cif.read(str(config / name)).sole_block()
        expected = gemmi.cif.read_string(text).sole_block()
        categories = expected.get_mmcif_category_names()
        require(actual.get_mmcif_category_names() == categories
                and all(actual.get_mmcif_category(cat) == expected.get_mmcif_category(cat)
                        for cat in categories), f"{name}: template differs from verified 6OGE crop")

    context = read_json(CONTEXT / "manifest.json")
    require(manifest["reference"]["sha256"] == context["sources"]["pdb_cif"]["sha256"],
            "Reference hash differs from frozen 1N8Z")
    queries = {q["query_id"]: chain_sequences(q) for q in expected_queries["queries"]}
    return manifest, context, canonical, queries, {p["query_id"]: p for p in panel}, deposit, staged, image, staged_paths


def validate_target_msa(config, expected_queries, template_hash):
    receipt = config / "manifest.json"
    manifest = read_json(receipt)
    alignment = config / "target.a3m"
    queries = config / QUERY_FILE
    require(manifest["template_sha256"] == template_hash, "MSA template differs from complete no-MSA crop")
    staged = {
        receipt.name: file_hash(receipt),
        queries.name: validate_file(queries, manifest["queries_sha256"]),
        alignment.name: validate_file(alignment, manifest["a3m_sha256"]),
        TEMPLATES["complete"]: template_hash,
    }
    target = expected_queries["queries"][0]["chains"][0]["sequence"]
    msa_queries = read_json(queries)
    alignment_path = msa_queries["queries"][0]["chains"][0]["msa_path"]
    require(isinstance(alignment_path, str) and Path(alignment_path).is_absolute()
            and Path(alignment_path).name == alignment.name,
            "Expected an explicit staged target.a3m path")
    expected = {"queries": [
        {**q, "chains": [{**c, "msa_path": alignment_path if c["chain_id"] == "T" else None}
                         for c in q["chains"]]} for q in expected_queries["queries"]
    ]}
    require(msa_queries == expected, "Target MSA queries change more than the target alignment path")

    with alignment.open() as stream:
        sequences = [str(row.seq) for row in SeqIO.parse(stream, "fasta")]
    require(len(sequences) == manifest["total_rows_including_query"] == 43
            and len(set(sequences)) == manifest["unique_aligned_sequences"] == 18
            and len(target) == manifest["query_length"]
            and sequences[0] == target and manifest["query_sha256"] == sha256(target.encode())
            and all(len(re.sub(r"[a-z.]", "", seq)) == len(target) for seq in sequences),
            "Target A3M alignment or query does not match the frozen manifest")
    return staged, alignment_path


def artifact_paths(query, seed, sample):
    # ASVS 5.3.2: Construct known paths from the frozen key. Receipt paths must match.
    stem = f"{query}/seed_{seed}/{query}_seed_{seed}_sample_{sample}"
    return {key: stem + suffix for key, suffix in SUFFIXES.items()}


def validate_run(folder, log_path, arm, directory_seed=None, staged_paths=None):
    result_path = folder / "result.json"
    config_path = folder / "config.json"
    run = read_json(result_path)
    provenance = run["provenance"]
    config = provenance["resolved_config"]
    seeds = config["seeds"]
    require(run["task"] == "predict" and run["schema_version"] == SCHEMA
            and provenance["module"] == {"name": "apheris-protenix", "version": "1.3.0"}
            and provenance["cli_contract_version"] == SCHEMA
            and read_json(config_path) == {"schema_version": SCHEMA, **config},
            f"{folder}: inconsistent runner metadata")
    require(bool(seeds) and len(seeds) == len(set(seeds)) and set(seeds) <= set(SEEDS)
            and (directory_seed is None or seeds == [directory_seed]),
            f"{folder}: unexpected seeds {seeds}")
    require(config["num_diffusion_samples"] == len(SAMPLES)
            and config["write_detailed_scores"] is True
            and config["query_failures"] == config["missing_msas"] == "fail"
            and config["weights"] == "/app/weights/protenix-v2.pt"
            and len(config["template_paths"]) == len(provenance["inputs"]) == 1
            and Path(config["template_paths"][0]).name == TEMPLATES[arm]
            and Path(provenance["inputs"][0]["path"]).name == QUERY_FILE
            and provenance["inputs"][0]["tier"] == "none",
            f"{folder}: unexpected model, inputs, template or sample count")
    if staged_paths is not None:
        require(config["template_paths"] == [staged_paths[TEMPLATES[arm]]]
                and provenance["inputs"][0]["path"] == staged_paths[QUERY_FILE],
                f"{folder}: runtime paths differ from the input hash receipt")

    predictions = run["result"]["predictions"]
    expected = list(itertools.product(PANEL, seeds, SAMPLES))
    keys = [(p["query_id"], p["seed"], p["sample"]) for p in predictions]
    exact_coverage(keys, expected, f"{folder} prediction receipt")
    model_paths = {}
    all_paths = {key: [] for key in SUFFIXES}
    for prediction, (query, seed, sample) in zip(predictions, keys, strict=True):
        require(prediction["query_index"] == PANEL.index(query), f"{folder}: incorrect query index")
        paths = artifact_paths(query, seed, sample)
        for key, path in paths.items():
            require(prediction[key] == path, f"{folder}: unexpected {key} for {(query, seed, sample)}")
            require((folder / path).resolve().is_relative_to(folder.resolve()),
                    f"{folder}: pose artifact escapes its run directory")
            all_paths[key].append(path)
        model_paths[arm, query, seed, sample] = {k: folder / v for k, v in paths.items()}
    for key, suffix in SUFFIXES.items():
        exact_coverage([str(p.relative_to(folder)) for p in folder.rglob(f"*{suffix}")],
                       all_paths[key], f"{folder} {key} artifacts")

    log = log_path.read_text()
    assigned = set(re.findall(r"Assigned templates to (\d+) entities across (\d+) queries", log))
    resolved = set(re.findall(r"Identified (\d+) resolved residue\(s\) from", log))
    polymer = set(re.findall(r"Read (\d+) chains? with a total of (\d+) residues from `_pdbx_poly_seq_scheme`", log))
    atom_counts = set(re.findall(r"Read (\d+) atoms from `_atom_site` table", log))
    initialization_seeds = {int(seed) for seed in re.findall(r"Seed set to (\d+)", log)}
    constraints = set(re.findall(r"Loaded constraint feature: #atom contact:(\d+) #contact:(\d+) #pocket:(\d+)", log))
    expected_coverage = len(POSITIONS) if arm == "complete" else len(COMMON)
    n = str(len(PANEL))
    # Lightning logs the first configured seed during global setup.
    # The pinned model's predict_step reseeds from batch["seed"] for each prediction.
    require(assigned == {(n, n)} and polymer == {("1", str(len(POSITIONS)))}
            and resolved == {str(expected_coverage)} and atom_counts == {"912" if arm == "complete" else "836"}
            and initialization_seeds == {seeds[0]} and constraints == {("0", "0", "0")},
            f"{folder}: template, seed or constraint log did not match the protocol")
    return model_paths, {
        "arm": arm, "seeds": seeds, "resolved_template_residues": expected_coverage,
        "template_assignments": len(PANEL), "pose_count": len(predictions),
        **({"lightning_initialization_seed": seeds[0]} if len(seeds) > 1 else {}),
        "receipt_sha256": {p.name: file_hash(p) for p in (result_path, config_path, log_path)},
        "config": config,
    }


def collect_runs(campaign, arms=TEMPLATES, staged_paths=None):
    layouts = {}
    for arm in arms:
        arm_folder = campaign / arm
        layouts[arm_folder / "result.json"] = (arm_folder, campaign / f"{arm}.log", arm, None)
        for seed in SEEDS:
            folder = arm_folder / f"seed-{seed}"
            layouts[folder / "result.json"] = (folder, folder.with_suffix(".log"), arm, seed)
    receipts = sorted(campaign.rglob("result.json"))
    require(set(receipts) <= layouts.keys(), "Unexpected Protenix run receipt location")
    models, runs = {}, {}
    for receipt in receipts:
        folder, *details = layouts[receipt]
        paths, info = validate_run(folder, *details, staged_paths=staged_paths)
        require(not models.keys() & paths.keys(), f"{folder}: duplicate arm/query/seed/sample across runs")
        models.update(paths)
        runs[str(folder.relative_to(campaign))] = info

    all_files = {p for paths in models.values() for p in paths.values()}
    for suffix in SUFFIXES.values():
        exact_coverage(campaign.rglob(f"*{suffix}"), [p for p in all_files if p.name.endswith(suffix)],
                       "Raw poses and receipts differ")
    common_configs = [
        {key: value for key, value in run["config"].items() if key not in {"seeds", "template_paths"}}
        for run in runs.values()
    ]
    require(not common_configs or all(c == common_configs[0] for c in common_configs),
            "Runner settings changed across template arms or seeds")
    return models, runs


def ecd_reference(deposit, canonical):
    atoms = source_target(deposit, canonical)
    ca = {int(a["auth_seq_id"]): [float(a[f"Cartn_{axis}"]) for axis in "xyz"]
          for a in atoms if a["label_atom_id"] == "CA"}
    require(set(COMMON) <= ca.keys(), "6OGE lacks a common C-alpha position")
    outside = [a for a in atoms if int(a["auth_seq_id"]) not in POSITIONS
               and a["type_symbol"] not in {"H", "D"}]
    coords = np.array([[float(a[f"Cartn_{axis}"]) for axis in "xyz"] for a in outside])
    ids = np.array([int(a["auth_seq_id"]) for a in outside])
    return np.array([ca[p] for p in COMMON]), cKDTree(coords), ids


def overlay_clashes(mapped, context):
    ref_ca, tree, ids = context
    target = ca_coordinates(mapped[mapped.chain_id == "T"])
    moving = np.array([target[p] for p in COMMON], dtype=float)
    rotation, shift = target_transform(moving, ref_ca)
    is_antibody = np.isin(mapped.chain_id, ["H", "L"]) & ~np.isin(mapped.element, ["H", "D"])
    coords = mapped.coord[is_antibody] @ rotation + shift
    hits = cKDTree(coords).sparse_distance_matrix(tree, 1.5, output_type="coo_matrix")
    clash = hits.data < 1.5
    return {
        "ecd_6oge_common_target_ca_rmsd_A": rmsd(moving @ rotation + shift, ref_ca),
        "ecd_6oge_outside_crop_heavy_atom_clash_pairs_below_1p5A": int(clash.sum()),
        "ecd_6oge_outside_crop_clash_positions_p04626": sorted(set(map(int, ids[hits.col[clash]]))),
        "ecd_6oge_outside_crop_minimum_distance_A": float(tree.query(coords)[0].min()),
    }


def score_pose(key, files, queries, panel, reference, ecd, campaign):
    arm, query, seed, sample = key
    atoms, roles = validate_prediction(files["path"], queries[query], protenix=True)
    require(tuple(roles) == tuple(roles.values()) == ROLES, "Unexpected Protenix confidence chain order")
    summary = read_json(files["scores_path"])
    tokens = int(np.count_nonzero(atoms.atom_name == "CA"))
    require((summary["num_chains"], summary["num_atoms"], summary["num_tokens"], summary["num_residues"])
            == (len(ROLES), len(atoms), tokens, tokens), "Protenix structure and confidence counts differ")
    overall = float(summary["interface_predicted_tm_score"])
    clash = summary["clash_flag"]
    require(0 <= overall <= 1 and isinstance(clash, bool), "Invalid global confidence or clash flag")
    pae = protenix_pae(files["detailed_scores_path"], atoms, tokens)
    matrix = np.asarray(summary["chain_chain_interface_predicted_tm_score"], dtype=float)
    require(antibody_mapping(query, panel), "Reference antibody position map is unavailable")
    mapped, geometry = geometry_row(atoms, roles, reference, MAPPING, matched=True)
    shared_count = len(set(geometry["target_contacts_p04626"]) & set(COMMON))
    pairs = score_pairs(ROLES, matrix)
    # Preserve ties in the exported decimal confidences.
    pair_mean = sum(Decimal(str(pairs[k])) for k in ("iptm_H_T", "iptm_L_T")) / 2
    return {
        "arm": arm, "query_id": query, "seed": seed, "sample": sample,
        "cdrh3": panel[query]["cdrh3"], "spr_kd_nM": panel[query].get("kd_nM"),
        "facs_class": panel[query].get("facs_class"),
        "cif_chain_roles": roles, "confidence_chain_order": list(ROLES),
        "overall_iptm": overall, "provider_clash_flag": clash,
        "row_artifact_sha256": {k: file_hash(path) for k, path in files.items()},
        "row_artifact_paths": {k: path.relative_to(campaign).as_posix() for k, path in files.items()},
        "wt_footprint_precision": geometry["wt_contact_count"] / shared_count if shared_count else None,
        "mean_antibody_target_pair_iptm": float(pair_mean),
        **pairs, **geometry,
        **footprint_pae(pae, mapped, roles, reference[4], reference[3]), **overlay_clashes(mapped, ecd),
    }


def rank_samples(samples):
    groups = []
    for arm, query, seed in itertools.product(TEMPLATES, PANEL, SEEDS):
        candidates = [row for row in samples
                      if (row["arm"], row["query_id"], row["seed"]) == (arm, query, seed)]
        if not candidates:
            continue
        by_confidence = sorted(candidates, key=lambda r: (-r["mean_antibody_target_pair_iptm"], r["sample"]))
        by_geometry = sorted(candidates, key=lambda r: (r["ref_interface_ca_rmsd_in_target_frame_A"], r["sample"]))
        for rank, row in enumerate(by_confidence, 1):
            row["confidence_rank_in_seed"] = rank
        for rank, row in enumerate(by_geometry, 1):
            row["reference_geometry_rank_in_seed_diagnostic"] = rank
        fields = ("sample", "confidence_rank_in_seed", "mean_antibody_target_pair_iptm",
                  "target_fit_ca_rmsd_A", "ref_interface_ca_rmsd_in_target_frame_A",
                  "wt_contact_count", "wt_footprint_recall", "wt_footprint_jaccard",
                  "heavy_atom_clash_pairs_below_1p5A", "ecd_6oge_outside_crop_heavy_atom_clash_pairs_below_1p5A",
                  "operational_conditions", "operational_pose_pass")
        groups.append({
            "arm": arm, "query_id": query, "seed": seed,
            "top1": {k: by_confidence[0][k] for k in fields},
            "best_reference_geometry_diagnostic": {k: by_geometry[0][k] for k in fields},
        })
    return groups


def operational_conditions(row):
    """Apply the frozen known-pose rule to each draw, regardless of selection."""
    return {
        "recovers_at_least_13_of_19": row["wt_contact_count"] >= 13,
        "target_ca_rmsd_at_most_2A": row["target_fit_ca_rmsd_A"] <= 2,
        "interface_ca_rmsd_at_most_4A": row["ref_interface_ca_rmsd_in_target_frame_A"] <= 4,
        "no_antibody_target_heavy_atom_clashes_below_1p5A": row["heavy_atom_clash_pairs_below_1p5A"] == 0,
    }


def operational_rule(groups, complete, eligible_arms=TEMPLATES):
    arms = {}
    for arm in eligible_arms:
        wt = [group for group in groups if group["arm"] == arm and group["query_id"] == "control_wt"]
        passes = [group["seed"] for group in wt if group["top1"]["operational_pose_pass"]]
        arms[arm] = {"wt_top1_seeds_observed": [group["seed"] for group in wt],
                     "passing_wt_top1_seeds": passes,
                     "qualifies": len(passes) >= 2 if complete else None}
    return {
        "source": "docs/research/structure-admission-rule.md",
        "status": "evaluated" if complete else "incomplete",
        "minimum_passing_seeds": 2, "arms": arms,
        "scope": "Operational heuristic for WT geometry reproduction.",
    }


def compare_ranked_pairs(groups):
    """Compare complete-minus-masked geometry at the same query and seed."""
    indexed = {(g["arm"], g["query_id"], g["seed"]): g["top1"] for g in groups}
    fields = ("wt_contact_count", "target_fit_ca_rmsd_A", "ref_interface_ca_rmsd_in_target_frame_A",
              "mean_antibody_target_pair_iptm", "heavy_atom_clash_pairs_below_1p5A")
    pairs = []
    for query, seed in itertools.product(PANEL, SEEDS):
        full, masked = (indexed.get((arm, query, seed)) for arm in TEMPLATES)
        if full is None or masked is None:
            continue
        pairs.append({"query_id": query, "seed": seed, "complete_sample": full["sample"],
                      "masked_sample": masked["sample"],
                      "complete_minus_masked": {field: full[field] - masked[field] for field in fields}})
    return pairs


def audit(campaign, reference_path, template_source, config=CONFIG, msa_config=None):
    manifest, context, canonical, queries, panel, deposit, staged, image, staged_paths = validate_inputs(
        config, campaign, template_source, msa_config)
    arms = ("complete",) if msa_config is not None else tuple(TEMPLATES)
    expected = [key for key in EXPECTED if key[0] in arms]
    models, runs = collect_runs(campaign, arms, staged_paths)
    refs = read_reference(reference_path, context, canonical, queries["control_wt"], {"split": {"T": COMMON}})
    require(refs[5] == COMMON and set(refs[1]) & set(POSITIONS) == set(COMMON)
            and len(refs[4]) == 19, "Unexpected reference target or footprint coverage")
    ecd = ecd_reference(deposit, canonical)
    missing = [key for key in expected if key not in models]
    samples = [score_pose(key, models[key], queries, panel, refs, ecd, campaign)
               for key in expected if key in models]
    for row in samples:
        row["msa_policy"] = "target_only" if msa_config is not None else "none"
        conditions = operational_conditions(row)
        row["operational_conditions"] = conditions
        row["operational_pose_pass"] = all(conditions.values())
    groups = rank_samples(samples)
    return {
        "status": "incomplete" if missing else "complete",
        "expected_pose_count": len(expected),
        "audited_pose_count": len(samples),
        "missing": [dict(zip(("arm", "query_id", "seed", "sample"), key)) for key in missing],
        "inputs_sha256": staged, "protenix_image_id": image,
        "data_origin": manifest["data_origin"],
        "msa_policy": "target_only" if msa_config is not None else "none",
        "selection_policy": {
            "group_by": ["arm", "query_id", "seed"],
            "score": "mean of VH–HER2 and VL–HER2 pair ipTM",
            "matrix_order": list(ROLES), "tie_break": "lower sample index",
            "reference_geometry_diagnostic": "minimum reference interface C-alpha RMSD after HER2 fit",
        },
        "template_assignment_evidence": "Aggregate runner logs and verified target-only input files. "
                                        "Per-query applied template masks are not in the export.",
        "reference": {
            "1n8z_sha256": refs[-1], "target_canonical_offset": 22,
            "antibody_role_to_chain": REFERENCE_CHAINS, "target_fit_positions_p04626": COMMON,
            "wt_footprint_positions_p04626": sorted(refs[4]), "antibody_interface_positions": refs[3],
            "6oge_sha256": SOURCE_SHA256, "6oge_outside_crop_heavy_atoms": len(ecd[2]),
        },
        "runs": runs,
        "selection": groups,
        "paired_top1_complete_minus_masked": compare_ranked_pairs(groups),
        "operational_rule": operational_rule(groups, complete=not missing, eligible_arms=arms),
        "predictions": samples,
    }


def write_audit(report, output):
    output.mkdir(parents=True, exist_ok=True)
    with (output / "audit.json").open("w") as stream:
        json.dump(report, stream, indent=2, allow_nan=False)
        stream.write("\n")
    records = report["predictions"]
    if records:
        fields = [k for k in records[0] if k != "operational_conditions"]
        fields.extend(records[0]["operational_conditions"])
    else:
        fields = ["arm", "query_id", "seed", "sample"]
    with (output / "scores.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fields, lineterminator="\n")
        writer.writeheader()
        for row in records:
            flat = {k: value for k, value in row.items() if k != "operational_conditions"}
            flat.update(row["operational_conditions"])
            writer.writerow({
                k: json.dumps(value, separators=(",", ":")) if isinstance(value, (list, dict)) else value
                for k, value in flat.items()
            })


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    campaign = parser.add_mutually_exclusive_group(required=True)
    campaign.add_argument("--campaign", type=Path, help="Root of extracted no-MSA Protenix outputs")
    campaign.add_argument("--msa-campaign", type=Path, help="Root of the complete-template, target-only MSA run")
    parser.add_argument("--reference", type=Path, required=True, help="Original 1N8Z mmCIF")
    parser.add_argument("--template-source", type=Path, required=True, help="Original 6OGE mmCIF")
    parser.add_argument("--config", type=Path, default=CONFIG)
    parser.add_argument("--msa-config", type=Path, default=MSA_CONFIG)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    report = audit(args.msa_campaign or args.campaign, args.reference, args.template_source, args.config,
                   args.msa_config if args.msa_campaign else None)
    default_output = OUT / "target-msa" if args.msa_campaign else OUT
    write_audit(report, args.output or default_output)
    print(f"{report['status']}: {report['audited_pose_count']}/{report['expected_pose_count']} poses audited")
    return 2 if report["status"] == "incomplete" else 0


if __name__ == "__main__":
    raise SystemExit(main())
