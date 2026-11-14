"""Combine immutable neural runs without selecting an architecture or seed by SPR."""

import argparse
import json
from pathlib import Path

import numpy as np

from scripts import rank_spr_candidates as frozen
from scripts import train_her2_neural as neural
from src.her2_interaction_benchmark import load_frozen_features


def collect(runs, output, oof):
    if output.exists():
        raise FileExistsError(f"Output already exists: {output}")
    ids, oof_source = neural.load_primary_ids(oof)
    matrices, sources = {}, []
    for directory in runs:
        arrays, provenance = load_frozen_features(
            directory / "features.npz", directory / "metadata.json", ids,
        )
        metadata = provenance["extraction"]
        if metadata["status"] != "complete":
            raise ValueError(f"Incomplete neural run: {directory}")
        config = metadata["configuration"]
        overlap = metadata["audits"]["primary_spr_encoder_overlap"]
        counts = overlap["exact_cdrh3_match_count_by_split"]
        if counts["train"] or counts["val"]:
            raise ValueError(f"Encoder training overlap: {directory}")
        for split, expected in neural.FACS_SHA256.items():
            if metadata["inputs"][split]["sha256"] != expected:
                raise ValueError(f"FACS fingerprint mismatch: {directory}/{split}")
        if set(arrays) != {f"facs_{model}" for model in config["models"]}:
            raise ValueError(f"Feature/model mismatch: {directory}")
        augmentation = metadata.get("augmentation")
        prefix = "augmented" if augmentation else config["size"]
        renamed = []
        for key, values in arrays.items():
            model = key.removeprefix("facs_")
            name = f"{prefix}_{model}_seed{config['seed']}"
            if name in matrices:
                raise ValueError(f"Duplicate neural experiment: {name}")
            matrices[name] = values
            renamed.append(name)
        sources.append({"run": str(directory), "features": renamed, **provenance})
    if not matrices:
        raise ValueError("At least one completed neural run is required")
    output.mkdir(parents=True)
    path = output / "features.npz"
    np.savez_compressed(path, sequence_ids=np.asarray(ids, dtype=str), **matrices)
    metadata = {
        "ordered_sequence_ids": ids,
        "selection": "All supplied prespecified architectures and seeds; no SPR-based selection",
        "cohort": oof_source,
        "runs": sources,
        "code": frozen.source(Path(__file__)),
        "output": {"file": path.name, "sha256": frozen.source(path)["sha256"]},
    }
    (output / "metadata.json").write_text(json.dumps(metadata, indent=2, allow_nan=False) + "\n")
    print(json.dumps({"runs": len(runs), "features": list(matrices), "rows": len(ids)}))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--runs", nargs="+", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--oof", type=Path, default=neural.DEFAULT_OOF)
    args = parser.parse_args()
    collect(args.runs, args.output, args.oof)


if __name__ == "__main__":
    main()
