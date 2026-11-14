"""Generate experimental HER2 CDRH3 proposals and attach frozen FACS pseudolabels."""

import argparse
import csv
import hashlib
import io
import json
import os
import pickle
import platform
import random
import sys
import tempfile
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from zipfile import BadZipFile, ZipFile

import numpy as np

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from src import her2_benchmark as facs

ROOT = Path(__file__).resolve().parents[1]
SPLIT_ARCHIVE = Path("zenodo-10831512-trastuzumab-ngs/data.zip")
DEFAULT_OOF = ROOT / "results/spr-benchmark/oof_predictions.csv"
PRIMARY_COHORT = "Primary novelty cohort"
PRIMARY_COUNT = 158
# Same pinned input as scripts/extract_ablang_features.py.
OOF_SHA256 = "44dfc97781fb80652aa59bfb8226e3b71fc501053d12c1d6e66cf8ea1c180fbf"
CLASS_NAMES = ("low", "mid", "high")
MUTATION_COUNTS = (1, 2, 3)
POSITIONS = range(facs.LOOP_LENGTH)
ALTERNATIVES = {aa: facs.AMINO_ACIDS.replace(aa, "") for aa in facs.AMINO_ACIDS}
VALID_RESIDUES = frozenset(facs.AMINO_ACIDS)
BATCH_SIZE = 2048
TEACHER_WIDTHS = (2048, 2048, 256)
TEACHER_SELECTION = "highest validation high-class AP; ties use earliest epoch"


@dataclass(frozen=True, slots=True)
class Proposal:
    sequence_id: str
    parent_sequence: str
    mutation_count: int


def sha256(path):
    with path.open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def load_sequence_sources(data_dir, oof_path):
    """Read sequence columns only, while hashing the original bytes for provenance."""
    splits, inputs = {}, {}
    with ZipFile(data_dir / SPLIT_ARCHIVE) as archive:
        for name in facs.SPLIT_NAMES:
            frame, inputs[name] = facs.read_csv(
                f"{facs.SPLIT_PATH}{name}.csv", archive=archive, columns=("seq",)
            )
            values = facs.normalize_sequences(frame["seq"], name, length=facs.LOOP_LENGTH)
            splits[name] = frozenset(values)
            if not splits[name] or len(splits[name]) != len(values):
                raise ValueError(f"Expected nonempty, unique sequences in the FACS {name} split")

    oof, inputs["oof"] = facs.read_csv(oof_path, columns=("cohort", "cdrh3"))
    if inputs["oof"]["sha256"] != OOF_SHA256:
        raise ValueError("Frozen SPR OOF fingerprint mismatch")
    primary = facs.normalize_sequences(
        oof.loc[oof["cohort"].eq(PRIMARY_COHORT), "cdrh3"],
        PRIMARY_COHORT,
        length=facs.LOOP_LENGTH,
    )
    if len(primary) != PRIMARY_COUNT or not primary.is_unique:
        raise ValueError(f"Expected {PRIMARY_COUNT} unique primary SPR IDs")
    return splits, frozenset(primary), inputs


def exclusion_sets(splits, primary):
    return {
        "wild_type": {facs.WT_CDRH3},
        **{f"facs_{name}": splits[name] for name in facs.SPLIT_NAMES},
        "spr_primary": primary,
    }


def generate_proposals(splits, primary, *, count, seed, max_attempts=None):
    if count < 1:
        raise ValueError("Proposal count must be positive")
    if not splits["train"]:
        raise ValueError("The train split has no proposal parents")
    if max_attempts is None:
        max_attempts = count * 100

    parents = sorted(splits["train"])
    exclusions = exclusion_sets(splits, primary)
    rejections = dict.fromkeys(("cysteine", *exclusions, "duplicate"), 0)
    proposals, accepted_ids = [], set()
    attempts = 0
    rng = random.Random(seed)
    for start in range(0, count, len(MUTATION_COUNTS)):
        # Randomize the order within each block and retry the same count after rejection.
        # This preserves uniform mutation counts in the accepted data.
        block_size = min(len(MUTATION_COUNTS), count - start)
        for mutation_count in rng.sample(MUTATION_COUNTS, k=block_size):
            while True:
                if attempts >= max_attempts:
                    raise ValueError(
                        f"Proposal attempt limit reached ({attempts}); accepted "
                        f"{len(proposals)} of {count} with rejections {rejections}"
                    )
                attempts += 1
                parent = rng.choice(parents)
                residues = list(parent)
                for position in rng.sample(POSITIONS, k=mutation_count):
                    residues[position] = rng.choice(ALTERNATIVES[parent[position]])
                sequence = "".join(residues)
                exclusion = next((key for key, ids in exclusions.items() if sequence in ids), None)
                if "C" in sequence:
                    rejections["cysteine"] += 1
                elif exclusion is not None:
                    rejections[exclusion] += 1
                elif sequence in accepted_ids:
                    rejections["duplicate"] += 1
                else:
                    accepted_ids.add(sequence)
                    proposals.append(Proposal(sequence, parent, mutation_count))
                    break

    return proposals, {
        "attempts": attempts,
        "max_attempts": max_attempts,
        "rejections_first_match": rejections,
        "accepted_by_mutation_count": dict(sorted(Counter(p.mutation_count for p in proposals).items())),
    }


def audit_proposals(proposals, splits, primary, expected_count):
    ids = {p.sequence_id for p in proposals}
    audit = {
        "duplicate_rows": len(proposals) - len(ids),
        "cysteine_rows": sum("C" in p.sequence_id for p in proposals),
        "invalid_sequence_rows": sum(
            len(sequence) != facs.LOOP_LENGTH or not set(sequence) <= VALID_RESIDUES
            for sequence in ids
        ),
        "invalid_parent_rows": sum(p.parent_sequence not in splits["train"] for p in proposals),
        "incorrect_mutation_count_rows": sum(
            p.mutation_count not in MUTATION_COUNTS
            or p.mutation_count != sum(a != b for a, b in zip(p.sequence_id, p.parent_sequence))
            for p in proposals
        ),
        "exact_overlaps": {
            name: len(ids & excluded) for name, excluded in exclusion_sets(splits, primary).items()
        },
    }
    if len(proposals) != expected_count or any(audit["exact_overlaps"].values()) or any(
        value for key, value in audit.items() if key != "exact_overlaps"
    ):
        raise ValueError(f"Proposal exclusion audit failed for {len(proposals)} rows: {audit}")
    return audit


def audit_probabilities(probabilities, count):
    if (
        probabilities.shape != (count, len(CLASS_NAMES))
        or probabilities.dtype != np.float32
        or not np.isfinite(probabilities).all()
    ):
        raise ValueError("Expected finite float32 teacher probabilities of shape (proposal count, 3)")
    if (probabilities < 0).any() or (probabilities > 1).any() or not np.allclose(
        probabilities.sum(axis=1), 1.0, atol=1e-6, rtol=0
    ):
        raise ValueError("Teacher probabilities must be in [0, 1] and sum to one per row")
    values = probabilities.astype(np.float64)
    log_values = np.log(values, where=values > 0, out=np.zeros_like(values))
    entropy = -np.sum(values * log_values, axis=1)
    return {
        "mean_by_class": dict(zip(CLASS_NAMES, values.mean(axis=0).tolist())),
        "argmax_count_by_class": dict(
            zip(CLASS_NAMES, np.bincount(values.argmax(axis=1), minlength=len(CLASS_NAMES)).tolist())
        ),
        "entropy_nats": dict(zip(("p05", "median", "p95"), np.percentile(entropy, [5, 50, 95]).tolist())),
        "max_row_sum_error": float(np.max(np.abs(values.sum(axis=1) - 1.0))),
    }


def load_teacher(teacher_run, inputs, *, device):
    from scripts import train_her2_neural as neural

    metadata_path = teacher_run / "metadata.json"
    metadata_bytes = metadata_path.read_bytes()
    metadata = json.loads(metadata_bytes)
    if metadata["status"] != "complete" or metadata["schema_version"] != 1:
        raise ValueError("Expected a complete neural FACS teacher run with schema version 1")
    config = metadata["configuration"]
    mlp = metadata["models"]["mlp"]
    if (
        config["size"] != "large"
        or "mlp" not in config["models"]
        or tuple(config["architecture"]["mlp"]) != TEACHER_WIDTHS
        or mlp["latent_width"] != TEACHER_WIDTHS[-1]
        or config["class_to_index"] != dict(zip(CLASS_NAMES, range(len(CLASS_NAMES))))
        or config["class_to_index"] != neural.CLASS_TO_INDEX
        or config["residue_alphabet"] != facs.AMINO_ACIDS
        or config["checkpoint_selection"] != TEACHER_SELECTION
    ):
        raise ValueError("Teacher architecture, residue alphabet, class order or checkpoint selection mismatch")

    # Verify the producer code that fits the MLP only on the published FACS train split.
    for name, module in (("trainer", neural), ("facs_benchmark", facs)):
        if metadata["code"][name]["sha256"] != sha256(Path(module.__file__).resolve()):
            raise ValueError(f"Teacher {name} source fingerprint mismatch with the real-data-only training code")
    for name in facs.SPLIT_NAMES:
        actual = inputs[name]["sha256"]
        if actual != neural.FACS_SHA256[name] or metadata["inputs"][name]["sha256"] != actual:
            raise ValueError(f"Teacher FACS {name} split fingerprint mismatch")
    if metadata["inputs"]["facs"]["sha256"] != neural.FACS_SHA256["facs"]:
        raise ValueError("Teacher FACS mapping fingerprint mismatch")
    if metadata["inputs"]["frozen_oof"]["sha256"] != inputs["oof"]["sha256"]:
        raise ValueError("Teacher SPR OOF fingerprint mismatch")

    best_checkpoint = mlp["best_checkpoint"]
    history = mlp["history"]
    if not history or not all(np.isfinite(record["validation_high_ap"]) for record in history):
        raise ValueError("Teacher history must contain finite validation average precision")
    selected = neural.select_best_epoch(history)
    if (
        best_checkpoint["file"] != "best.pt"
        or best_checkpoint["sha256"] != selected["checkpoint"]["sha256"]
        or mlp["best_epoch"] != selected["epoch"]
    ):
        raise ValueError("Teacher checkpoint differs from the validation-selected MLP checkpoint")
    checkpoint_path = teacher_run / "mlp" / "best.pt"
    checkpoint_bytes = checkpoint_path.read_bytes()
    checkpoint_hash = hashlib.sha256(checkpoint_bytes).hexdigest()
    if best_checkpoint["sha256"] != checkpoint_hash:
        raise ValueError("Teacher checkpoint fingerprint mismatch")

    neural.configure_torch(config["seed"], device)
    torch = neural.torch
    model = neural.build_model("mlp", "large")
    state = torch.load(io.BytesIO(checkpoint_bytes), map_location="cpu", weights_only=True)
    if not isinstance(state, dict) or not all(
        isinstance(tensor, torch.Tensor) and tensor.dtype == torch.float32
        and torch.isfinite(tensor).all().item()
        for tensor in state.values()
    ):
        raise ValueError("Teacher checkpoint must contain finite float32 state_dict tensors")
    model.load_state_dict(state, strict=True)
    model.to(device).requires_grad_(False).eval()
    return model, {
        "metadata": {"source": str(metadata_path), "sha256": hashlib.sha256(metadata_bytes).hexdigest()},
        "checkpoint": {"source": str(checkpoint_path), "sha256": checkpoint_hash},
        "code": metadata["code"],
        "facs_mapping": metadata["inputs"]["facs"],
        "model": "large_mlp",
        "latent_width": TEACHER_WIDTHS[-1],
        "training_data": "Published FACS train split only",
        "checkpoint_selection": TEACHER_SELECTION,
        "selected_epoch": selected["epoch"],
        "seed": config["seed"],
        "torch_version": torch.__version__,
    }


def predict_probabilities(model, proposals, *, device):
    from scripts import train_her2_neural as neural

    torch = neural.torch
    model.eval()
    with torch.no_grad():
        tokens = torch.as_tensor(neural.encode_tokens(p.sequence_id for p in proposals), device=device)
        batches = [
            model(batch).softmax(dim=1, dtype=torch.float32).cpu()
            for batch in tokens.split(BATCH_SIZE)
        ]
        return torch.cat(batches).numpy()


def write_dataset(output, proposals, probabilities, metadata):
    output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".her2-synthetic-", dir=output.parent) as temporary:
        staging = Path(temporary) / "complete"
        staging.mkdir()
        npz = staging / "synthetic.npz"
        np.savez_compressed(
            npz,
            sequence_ids=np.array([p.sequence_id for p in proposals], dtype=f"<U{facs.LOOP_LENGTH}"),
            probabilities=probabilities,
        )
        csv_path = staging / "proposals.csv"
        with csv_path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.writer(handle, lineterminator="\n")
            writer.writerow([
                "sequence_id", "parent_sequence", "mutation_count",
                *(f"p_facs_{name}" for name in CLASS_NAMES),
            ])
            for proposal, row in zip(proposals, probabilities, strict=True):
                writer.writerow([
                    proposal.sequence_id, proposal.parent_sequence, proposal.mutation_count, *row,
                ])
        metadata["output"] = {
            "file": npz.name,
            "sha256": sha256(npz),
            "csv": {"file": csv_path.name, "sha256": sha256(csv_path)},
        }
        (staging / "metadata.json").write_text(
            json.dumps(metadata, indent=2, allow_nan=False) + "\n", encoding="utf-8"
        )
        # Reserve the exact path before publishing, including against a concurrent invocation.
        output.mkdir()
        try:
            staging.replace(output)
        except OSError:
            output.rmdir()
            raise


def run(args):
    if os.path.lexists(args.output):
        raise FileExistsError(f"Output already exists: {args.output}")
    splits, primary, inputs = load_sequence_sources(args.data_dir, args.oof)
    proposals, generation = generate_proposals(splits, primary, count=args.count, seed=args.seed)
    audit = audit_proposals(proposals, splits, primary, args.count)
    model, teacher = load_teacher(args.teacher_run, inputs, device=args.device)
    probabilities = predict_probabilities(model, proposals, device=args.device)
    metadata = {
        "schema_version": 1,
        "label_type": "pseudolabels",
        "label_description": "Teacher FACS class probabilities for unvalidated sequence proposals.",
        "measured_endpoints": [],
        "class_order": list(CLASS_NAMES),
        "seed": args.seed,
        "count": args.count,
        "device": args.device,
        "inputs": inputs,
        "teacher": teacher,
        "environment": {"python": platform.python_version(), "numpy": np.__version__},
        "generator_source": {"source": str(Path(__file__).resolve()), "sha256": sha256(Path(__file__))},
        "generation": {
            **generation,
            "parent_split": "train",
            "sampling": "Uniform parents, distinct positions and 19 alternate residues per position. "
                        "Balanced mutation counts in randomized blocks of up to three accepted rows.",
            "exclusion_set_sizes": {name: len(ids) for name, ids in exclusion_sets(splits, primary).items()},
            "rejection_priority": list(generation["rejections_first_match"]),
        },
        "exclusion_audit": audit,
        "probability_audit": audit_probabilities(probabilities, args.count),
    }
    write_dataset(args.output, proposals, probabilities, metadata)
    print(args.output / "metadata.json")


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", required=True, type=Path, help="Organizer hackathon-data directory")
    parser.add_argument("--oof", type=Path, default=DEFAULT_OOF, help="Pinned SPR OOF CSV")
    parser.add_argument("--teacher-run", required=True, type=Path, help="Completed, real-only large MLP run")
    parser.add_argument("--output", required=True, type=Path, help="New directory for this dataset")
    parser.add_argument("--count", type=int, default=100_000)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--device", choices=("cpu", "cuda"), default="cpu")
    args = parser.parse_args(argv)
    if args.count < 1:
        parser.error("--count must be positive")
    try:
        run(args)
    except (ValueError, OSError, KeyError, BadZipFile, RuntimeError, pickle.UnpicklingError) as exc:
        parser.exit(2, f"error: {exc}\n")


if __name__ == "__main__":
    main()
