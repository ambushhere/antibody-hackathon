"""Train FACS classifiers and export frozen CDRH3 features for the primary SPR cohort."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import platform
import random
import sys
from contextlib import contextmanager
from importlib.metadata import version
from pathlib import Path
from shutil import copyfileobj
from tempfile import TemporaryDirectory
from time import perf_counter
from zipfile import BadZipFile, ZipFile

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from src import her2_benchmark as facs

try:
    import torch
    from torch import nn
    from torch.nn import functional as F
except ModuleNotFoundError as exc:
    if exc.name != "torch":
        raise
    torch = None

ROOT = Path(__file__).resolve().parents[1]
SPLIT_ARCHIVE = Path("zenodo-10831512-trastuzumab-ngs/data.zip")
FACS_CSV = "tz_her2_facs_3classes_524346.csv"
FACS_ARCHIVE = Path("her2-trastuzumab") / f"{FACS_CSV}.zip"
# Uncompressed input hashes from results/spr-benchmark/metrics.json.
FACS_SHA256 = {
    "train": "b639077ff116daa8b8c3a2122b6ed9291361892469f6237963967ca0772cb262",
    "val": "ea5e2e0908f018f579e476cef4452f9de11ba9144513bffa783588d49123f5d1",
    "test": "b350c10c7df4eb5e1b966f76909a805fc7d8011486f1ddb5dd1d7a6e19556a06",
    "facs": "7ff5a95bcc1942bb6e34b1e4ed6c02a9c177ee30ca160440f8abe82ed787735e",
}
DEFAULT_OOF = ROOT / "results/spr-benchmark/oof_predictions.csv"
# Same frozen OOF source as scripts/extract_ablang_features.py.
OOF_SHA256 = "44dfc97781fb80652aa59bfb8226e3b71fc501053d12c1d6e66cf8ea1c180fbf"
PRIMARY_COHORT = "Primary novelty cohort"
PRIMARY_COUNT = 158
CLASS_TO_INDEX = {"low": 0, "mid": 1, "high": 2}
HIGH = CLASS_TO_INDEX["high"]
N_RESIDUES = len(facs.AMINO_ACIDS)
N_CLASSES = len(CLASS_TO_INDEX)
BATCH_SIZE = 2048
LEARNING_RATE = 3e-4
WEIGHT_DECAY = 0.01
DROPOUT = 0.1
CPU_THREADS = 4
MODELS = ("transformer", "mlp")
ARCHITECTURES = {
    "small": {"transformer": (64, 2, 4, 128), "mlp": (256, 64)},
    "large": {"transformer": (256, 6, 8, 1024), "mlp": (2048, 2048, 256)},
}
MAX_EPOCHS = {"small": 5, "large": 15}
HISTORICAL_FACS_LOGISTIC = {
    "source": "docs/BENCHMARK.md",
    "original_metrics_sha256": "229f4e90912aad73e4923f50f59970a53b922da8d4ce1ef0b8e413773fe36a78",
    "selected_model": "pairwise logistic regression",
    "selection": "validation high-class average precision",
    "split": "published FACS test",
    "roc_auc": 0.9858926153128575,
    "pr_auc_average_precision": 0.973718568352088,
}


def file_source(path):
    with path.open("rb") as handle:
        return {"source": str(path), "sha256": hashlib.file_digest(handle, "sha256").hexdigest()}


def require_hash(provenance, expected):
    if provenance["sha256"] != expected:
        raise ValueError(f"Frozen input SHA-256 mismatch for {provenance['source']}")


def load_facs_splits(data_dir):
    """Load the published remove_overlap splits and verify the organizer FACS mapping."""
    splits, inputs = {}, {}
    with ZipFile(data_dir / SPLIT_ARCHIVE) as archive:
        for name in facs.SPLIT_NAMES:
            splits[name], inputs[name] = facs.load_split(f"{facs.SPLIT_PATH}{name}.csv", archive=archive)
    split_audit = facs.audit_splits(splits)
    for name, source in inputs.items():
        require_hash(source, FACS_SHA256[name])

    mapping = data_dir / FACS_ARCHIVE
    # The existing mapping audit takes a CSV path. Copy the named member byte-for-byte.
    with ZipFile(mapping) as archive, TemporaryDirectory(prefix="her2-facs-") as temporary:
        members = [
            entry for entry in archive.infolist()
            if Path(entry.filename).name == FACS_CSV and "__MACOSX" not in Path(entry.filename).parts
        ]
        if len(members) != 1:
            raise ValueError(f"{mapping}: expected exactly one {FACS_CSV} member")
        extracted = Path(temporary) / FACS_CSV
        with archive.open(members[0]) as src, extracted.open("wb") as dst:
            copyfileobj(src, dst)
        inputs["facs"], mapping_audit = facs.audit_facs(extracted, splits)
        inputs["facs"]["source"] = f"{mapping}!{members[0].filename}"
        require_hash(inputs["facs"], FACS_SHA256["facs"])
    return splits, inputs, {"split_audit": split_audit, "facs_mapping": mapping_audit}


def load_primary_ids(path):
    """Read only the ordered primary CDRH3 IDs from the frozen SPR benchmark."""
    rows, source = facs.read_csv(path, columns=("cohort", "cdrh3"))
    require_hash(source, OOF_SHA256)
    ids = rows.loc[rows["cohort"].eq(PRIMARY_COHORT), "cdrh3"].tolist()
    if len(ids) != PRIMARY_COUNT or len(set(ids)) != PRIMARY_COUNT:
        raise ValueError(f"Expected {PRIMARY_COUNT} unique primary OOF IDs")
    validated = facs.normalize_sequences(pd.Series(ids), "Primary SPR IDs", length=facs.LOOP_LENGTH)
    if validated.tolist() != ids or facs.WT_CDRH3 in ids:
        raise ValueError("Expected canonical primary SPR IDs distinct from WT")
    return ids, source


def audit_primary_overlap(ids, splits):
    """Require zero exact encoder train/validation matches for the frozen cohort."""
    membership = {}
    for name, split in splits.items():
        sequences = set(split["seq"])
        membership[name] = [seq for seq in ids if seq in sequences]
    if membership["train"] or membership["val"]:
        raise ValueError(
            "Primary SPR IDs overlap FACS train/val: "
            f"train={membership['train']}, val={membership['val']}"
        )
    return {
        "rows": len(ids),
        "exact_cdrh3_match_count_by_split": {name: len(matches) for name, matches in membership.items()},
        "exact_cdrh3_match_ids_by_split": membership,
    }


def encode_tokens(sequences):
    normalized = facs.normalize_sequences(
        pd.Series(sequences, dtype=str), "CDRH3 tokens", length=facs.LOOP_LENGTH
    )
    ascii_codes = np.frombuffer("".join(normalized).encode("ascii"), dtype=np.uint8)
    return facs.RESIDUE_CODES[ascii_codes].astype(np.int64).reshape(-1, facs.LOOP_LENGTH)


def require_torch():
    if torch is None:
        raise ValueError("PyTorch is required for neural FACS training and inference")


def encode_split(frame, device):
    require_torch()
    labels = frame["class"].map(CLASS_TO_INDEX)
    if labels.isna().any():
        raise ValueError("Expected FACS classes low, mid and high")
    return (
        torch.as_tensor(encode_tokens(frame["seq"]), device=device),
        torch.as_tensor(labels.to_numpy(dtype=np.int64, copy=True), device=device),
    )


def position_one_hot(tokens):
    require_torch()
    return F.one_hot(tokens, num_classes=N_RESIDUES).flatten(1).to(torch.float32)


if torch is not None:
    class FACSMLP(nn.Module):
        def __init__(self, widths):
            super().__init__()
            layers = []
            input_width = facs.LOOP_LENGTH * N_RESIDUES
            for width in widths:
                layers.extend((nn.Linear(input_width, width), nn.GELU(), nn.Dropout(DROPOUT)))
                input_width = width
            self.latent_width = widths[-1]
            self.hidden = nn.Sequential(*layers)
            self.head = nn.Linear(self.latent_width, N_CLASSES)

        def encode(self, tokens):
            return self.hidden(position_one_hot(tokens))

        def forward(self, tokens):
            return self.head(self.encode(tokens))


    class FACSTransformer(nn.Module):
        def __init__(self, width, layers, heads, ff_width):
            super().__init__()
            self.latent_width = width
            self.residues = nn.Embedding(N_RESIDUES, width)
            self.cls = nn.Parameter(torch.empty(1, 1, width))
            self.positions = nn.Parameter(torch.empty(1, facs.LOOP_LENGTH + 1, width))
            nn.init.normal_(self.cls, std=0.02)
            nn.init.normal_(self.positions, std=0.02)
            self.encoder = nn.Sequential(*[
                nn.TransformerEncoderLayer(
                    d_model=width, nhead=heads, dim_feedforward=ff_width, dropout=DROPOUT,
                    activation="gelu", batch_first=True, norm_first=True,
                ) for _ in range(layers)
            ])
            self.norm = nn.LayerNorm(width)
            self.head = nn.Linear(width, N_CLASSES)

        def encode(self, tokens):
            cls = self.cls.expand(tokens.shape[0], -1, -1)
            inputs = torch.cat((cls, self.residues(tokens)), dim=1)
            return self.norm(self.encoder(inputs + self.positions)[:, 0])

        def forward(self, tokens):
            return self.head(self.encode(tokens))


def build_model(name, size="small"):
    require_torch()
    if size not in ARCHITECTURES or name not in MODELS:
        raise ValueError(f"Unknown FACS architecture: {name}/{size}")
    if name == "transformer":
        return FACSTransformer(*ARCHITECTURES[size][name])
    return FACSMLP(ARCHITECTURES[size][name])


def configure_torch(seed, device):
    require_torch()
    if not 0 <= seed < 2**32:
        raise ValueError("Seed must be between 0 and 2**32 - 1")
    if torch.device(device).type == "cuda":
        os.environ["CUBLAS_WORKSPACE_CONFIG"] = ":4096:8"
        if not torch.cuda.is_available():
            raise ValueError("CUDA was requested but is unavailable")
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    torch.set_num_threads(CPU_THREADS)
    torch.set_default_dtype(torch.float32)
    torch.set_float32_matmul_precision("highest")
    torch.backends.cudnn.allow_tf32 = False
    torch.backends.cudnn.benchmark = False
    torch.use_deterministic_algorithms(True)


def train_one_epoch(model, tokens, labels, optimizer, *, batch_size=BATCH_SIZE, generator=None):
    if not len(tokens) or len(tokens) != len(labels) or batch_size < 1:
        raise ValueError("Training requires a positive batch size and nonempty aligned tokens and labels")
    model.train()
    order = torch.randperm(len(tokens), device=tokens.device, generator=generator)
    total = torch.zeros((), dtype=torch.float32, device=tokens.device)
    for batch in order.split(batch_size):
        optimizer.zero_grad(set_to_none=True)
        loss = F.cross_entropy(model(tokens[batch]), labels[batch])
        loss.backward()
        optimizer.step()
        total += loss.detach() * len(batch)
    mean_loss = float(total) / len(tokens)
    if not math.isfinite(mean_loss):
        raise ValueError("Nonfinite training loss")
    return mean_loss


def predict(model, tokens, *, batch_size=BATCH_SIZE, include_features=False):
    if not len(tokens) or batch_size < 1:
        raise ValueError("Prediction requires a positive batch size and nonempty tokens")
    model.eval()
    probabilities, features = [], []
    with torch.no_grad():
        finite = torch.ones((), dtype=torch.bool, device=tokens.device)
        for batch in tokens.split(batch_size):
            latent = model.encode(batch)
            logits = model.head(latent)
            finite &= torch.isfinite(logits).all() & torch.isfinite(latent).all()
            probabilities.append(logits.softmax(dim=1)[:, HIGH])
            if include_features:
                features.append(latent)
        if not finite.item():
            raise ValueError("Nonfinite FACS logits or embeddings")
        scores = torch.cat(probabilities).cpu().numpy()
        values = torch.cat(features).cpu().numpy() if include_features else None
    if not np.isfinite(scores).all() or (values is not None and not np.isfinite(values).all()):
        raise ValueError("Nonfinite FACS predictions or embeddings")
    return scores, values


@contextmanager
def staged_file(path):
    temporary = path.with_name(f".{path.name}.tmp")
    try:
        yield temporary
        temporary.replace(path)
    finally:
        temporary.unlink(missing_ok=True)


def write_json(path, payload):
    with staged_file(path) as temporary:
        temporary.write_text(json.dumps(payload, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def save_checkpoint(path, model):
    # Store only CPU tensors in the state dictionary. Restore with weights_only=True.
    state = {key: tensor.detach().cpu() for key, tensor in model.state_dict().items()}
    if not all(torch.isfinite(tensor).all().item() for tensor in state.values()):
        raise ValueError("Nonfinite checkpoint tensor")
    with staged_file(path) as temporary:
        with temporary.open("wb") as handle:
            torch.save(state, handle)
    return {"file": path.name, "sha256": file_source(path)["sha256"]}


def select_best_epoch(history):
    """Select using validation high-class AP alone, keeping the earliest exact tie."""
    return max(history, key=lambda record: (record["validation_high_ap"], -record["epoch"]))


def fit_model(model, train, val, *, epochs, seed, output_dir, batch_size=BATCH_SIZE, log=print):
    """Train on FACS train, select on FACS validation and restore the selected weights."""
    if epochs < 1:
        raise ValueError("Expected at least one epoch")
    output_dir.mkdir()
    generator = torch.Generator(device=train[0].device).manual_seed(seed)
    optimizer = torch.optim.AdamW(model.parameters(), lr=LEARNING_RATE, weight_decay=WEIGHT_DECAY)
    val_high = val[1].cpu().numpy() == HIGH
    if not val_high.any() or val_high.all():
        raise ValueError("FACS validation requires high and non-high labels")
    history = []
    parameter_count = sum(parameter.numel() for parameter in model.parameters())
    started = perf_counter()
    for epoch in range(1, epochs + 1):
        epoch_started = perf_counter()
        loss = train_one_epoch(model, *train, optimizer, batch_size=batch_size, generator=generator)
        scores, _ = predict(model, val[0], batch_size=batch_size)
        ap = float(average_precision_score(val_high, scores))
        if not math.isfinite(ap):
            raise ValueError("Nonfinite validation high-class AP")
        checkpoint = save_checkpoint(output_dir / f"epoch-{epoch:02d}.pt", model)
        history.append({
            "epoch": epoch,
            "train_cross_entropy": loss,
            "validation_high_ap": ap,
            "checkpoint": checkpoint,
            "elapsed_seconds": perf_counter() - epoch_started,
        })
        best = select_best_epoch(history)
        result = {
            "parameter_count": parameter_count,
            "latent_width": model.latent_width,
            "history": history,
            "best_epoch": best["epoch"],
            "best_validation_high_ap": best["validation_high_ap"],
            "best_checkpoint": best["checkpoint"],
            "training_seconds": perf_counter() - started,
        }
        write_json(output_dir / "history.json", result)
        log(f"{output_dir.name} epoch {epoch}/{epochs} loss={loss:.6f} val_high_ap={ap:.6f} "
            f"best_epoch={best['epoch']} elapsed_seconds={history[-1]['elapsed_seconds']:.2f}")

    best_path = output_dir / result["best_checkpoint"]["file"]
    state = torch.load(best_path, map_location="cpu", weights_only=True)
    model.load_state_dict(state)
    model.eval()
    with staged_file(output_dir / "best.pt") as temporary:
        os.link(best_path, temporary)
    result["best_checkpoint"] = {**result["best_checkpoint"], "file": "best.pt"}
    write_json(output_dir / "history.json", result)
    return result


def run(args):
    if os.path.lexists(args.output):
        raise FileExistsError(f"Output already exists: {args.output}")
    configure_torch(args.seed, args.device)
    splits, inputs, audits = load_facs_splits(args.data_dir)
    ids, inputs["frozen_oof"] = load_primary_ids(args.oof)
    audits["primary_spr_encoder_overlap"] = audit_primary_overlap(ids, splits)
    code = {
        "trainer": file_source(Path(__file__).resolve()),
        "facs_benchmark": file_source(Path(facs.__file__).resolve()),
    }
    metadata = {
        "schema_version": 1,
        "status": "training",
        "cohort": PRIMARY_COHORT,
        "ordered_sequence_ids": ids,
        "inputs": inputs,
        "code": code,
        "audits": audits,
        "configuration": {
            **vars(args),
            "data_dir": str(args.data_dir),
            "oof": str(args.oof),
            "output": str(args.output),
            "batch_size": BATCH_SIZE,
            "class_to_index": CLASS_TO_INDEX,
            "residue_alphabet": facs.AMINO_ACIDS,
            "loss": "unweighted three-class cross entropy",
            "optimizer": "AdamW",
            "learning_rate": LEARNING_RATE,
            "weight_decay": WEIGHT_DECAY,
            "dropout": DROPOUT,
            "cpu_threads": CPU_THREADS,
            "dtype": "float32",
            "deterministic_algorithms": True,
            "checkpoint_selection": "highest validation high-class AP; ties use earliest epoch",
            "architecture": ARCHITECTURES[args.size],
        },
        "versions": {
            "python": platform.python_version(),
            "platform": platform.platform(),
            **{name: version(name) for name in ("torch", "numpy", "pandas", "scipy", "scikit-learn")},
            "torch_cuda": torch.version.cuda,
        },
        "historical_facs_logistic": HISTORICAL_FACS_LOGISTIC,
    }
    if args.device == "cuda":
        metadata["configuration"]["cuda_device_name"] = torch.cuda.get_device_name()
        metadata["configuration"]["cublas_workspace_config"] = os.environ["CUBLAS_WORKSPACE_CONFIG"]
        metadata["versions"]["cudnn"] = torch.backends.cudnn.version()

    args.output.mkdir(parents=True)
    write_json(args.output / "run.json", metadata)

    def log(message):
        with (args.output / "training.log").open("a", encoding="utf-8") as handle:
            handle.write(message + "\n")
        print(message, flush=True)

    log(f"Training {', '.join(args.models)} ({args.size}, seed={args.seed}, device={args.device})")
    tensors = {name: encode_split(split, args.device) for name, split in splits.items()}
    primary_tokens = torch.as_tensor(encode_tokens(ids), device=args.device)
    test_tokens, test_labels = tensors["test"]
    test_high = test_labels.cpu().numpy() == HIGH
    features, diagnostics, metrics = {}, {"sequence_id": ids}, {}
    for name in args.models:
        configure_torch(args.seed, args.device)
        model = build_model(name, args.size).to(args.device)
        result = fit_model(
            model, tensors["train"], tensors["val"], epochs=args.epochs, seed=args.seed,
            output_dir=args.output / name, log=log,
        )
        model.requires_grad_(False)
        test_scores, _ = predict(model, test_tokens)
        result["facs_test"] = facs.classification_metrics(test_high, test_scores)
        write_json(args.output / name / "metrics.json", result)
        primary_scores, features[f"facs_{name}"] = predict(model, primary_tokens, include_features=True)
        diagnostics[f"p_high_{name}"] = primary_scores
        metrics[name] = result
        log(f"{name} selected epoch {result['best_epoch']} "
            f"test_high_ap={result['facs_test']['pr_auc_average_precision']:.6f}")

    features_path = args.output / "features.npz"
    with staged_file(features_path) as temporary:
        with temporary.open("wb") as handle:
            np.savez_compressed(handle, sequence_ids=np.asarray(ids, dtype=str), **features)
    diagnostics_path = args.output / "primary_probabilities.csv"
    with staged_file(diagnostics_path) as temporary:
        pd.DataFrame(diagnostics).to_csv(temporary, index=False)
    metadata.update({
        "status": "complete",
        "models": metrics,
        "output": {
            "file": features_path.name,
            "sha256": file_source(features_path)["sha256"],
            "feature_shapes": {name: list(matrix.shape) for name, matrix in features.items()},
            "primary_probabilities": file_source(diagnostics_path),
        },
    })
    write_json(args.output / "metadata.json", metadata)
    write_json(args.output / "run.json", metadata)
    log(f"Complete: {args.output}")


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", required=True, type=Path, help="Organizer hackathon-data directory")
    parser.add_argument("--oof", type=Path, default=DEFAULT_OOF, help="Frozen SPR OOF CSV")
    parser.add_argument("--output", required=True, type=Path, help="New directory for this run")
    parser.add_argument("--size", choices=ARCHITECTURES, default="small")
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--epochs", type=int, default=5)
    parser.add_argument("--models", nargs="+", choices=MODELS, default=list(MODELS))
    parser.add_argument("--device", choices=("cpu", "cuda"), default="cpu")
    args = parser.parse_args(argv)
    if not 1 <= args.epochs <= MAX_EPOCHS[args.size]:
        parser.error(f"--epochs must be between 1 and {MAX_EPOCHS[args.size]} for --size {args.size}")
    if len(args.models) != len(set(args.models)):
        parser.error("--models must list each architecture at most once")
    try:
        run(args)
    except (ValueError, OSError, KeyError, BadZipFile) as exc:
        parser.exit(2, f"error: {exc}\n")


if __name__ == "__main__":
    main()
