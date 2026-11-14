"""Inspect two pinned RFD4 metadata tables without loading structure shards."""

import argparse
import hashlib
import json
from pathlib import Path

REVISION = "aa86df1eb9e3529a0f39279f216dce3a95bb3dc2"
BASE = f"https://huggingface.co/datasets/nvidia/RFD4-Proteina/resolve/{REVISION}/"
FILES = {
    "immune_interfaces.parquet": {
        "path": "train_datasets/streaming-latest/pdb/2026_10_04_tmol_plip_cath/sampling/immune_interfaces.parquet",
        "bytes": 20926603,
        "sha256": "3ddd850105910a18bdf4273605d7d96c6a527468e533cc351c04722ea1774ca3",
    },
    "distilled_examples.parquet": {
        "path": "train_datasets/streaming-latest/synthetic/2026_09_18_antibody_antigen_esmfold2/examples.parquet",
        "bytes": 25488687,
        "sha256": "d4264ccca50832c8ed9edbe6c4c252de9a017069a71890551f573977b1b45f63",
    },
}


def audit(directory):
    import pyarrow.parquet as pq

    inputs, schemas = {}, {}
    for name, spec in FILES.items():
        path = directory / name
        if path.stat().st_size != spec["bytes"] or hashlib.sha256(path.read_bytes()).hexdigest() != spec["sha256"]:
            raise ValueError(f"Pinned input fingerprint mismatch: {name}")
        inputs[name] = {**spec, "url": BASE + spec["path"]}
        parquet = pq.ParquetFile(path)
        schemas[name] = {"rows": parquet.metadata.num_rows, "columns": parquet.schema_arrow.names}
    pdb = pq.read_table(directory / "immune_interfaces.parquet", columns=[
        "pdb_id", "ig_type", "assembly_id", "ig_pn_unit_iid", "antigen_pn_unit_iid",
        "deposition_date", "resolution", "passes_bond_filters", "passes_ligand_identity_filters",
    ]).to_pandas()
    distilled = pq.read_table(directory / "distilled_examples.parquet",
                              columns=["format", "source"]).to_pandas()
    hits = pdb[pdb.pdb_id.str.lower().isin(["1n8z", "6oge"])].sort_values(
        ["pdb_id", "assembly_id", "ig_pn_unit_iid", "antigen_pn_unit_iid"])
    return {
        "revision": REVISION, "inputs": inputs, "schemas": schemas,
        "downloaded_bytes": sum(row["bytes"] for row in inputs.values()),
        "pdb_immune_types": pdb.ig_type.value_counts().to_dict(),
        "distilled_formats": distilled["format"].value_counts().to_dict(),
        "distilled_sources": distilled.source.value_counts().to_dict(),
        "reference_entries_in_sampling_table": json.loads(hits.to_json(orient="records", date_format="iso")),
        "interpretation": [
            "Rows are sampling entries, not unique independent complexes or measured binders.",
            "Reference presence is in the published sampling table; exact checkpoint training exposure is not established.",
            "No release source filters were applied by this audit. Passing counts in the upstream card are different.",
            "Distillation metadata describe predicted structures and confidence, not measured KD or Tm.",
            "These training tables cannot be presented as independent tests of models trained on them.",
        ],
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--metadata-dir", type=Path, default=Path("data/external/rfd4-metadata"))
    parser.add_argument("--output", type=Path, default=Path("results/rfd4-training-metadata/audit.json"))
    args = parser.parse_args()
    result = audit(args.metadata_dir)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    print(json.dumps({k: result[k] for k in ("downloaded_bytes", "pdb_immune_types", "distilled_formats")}))


if __name__ == "__main__":
    main()
