"""Train the prespecified large Transformer with separate teacher pseudolabel loss."""

import argparse
import copy
import json
import math
import os
from pathlib import Path
from time import perf_counter

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score

from scripts import train_her2_neural as base
from src import her2_benchmark as facs

SYNTHETIC_WEIGHT = 0.2


def load_synthetic(directory, splits, primary_ids, device):
    metadata = json.loads((directory / "metadata.json").read_text())
    path = directory / "synthetic.npz"
    if base.file_source(path)["sha256"] != metadata["output"]["sha256"]:
        raise ValueError("Synthetic array fingerprint mismatch")
    if (metadata["label_type"] != "pseudolabels" or metadata["measured_endpoints"]
            or metadata["class_order"] != list(base.CLASS_TO_INDEX)
            or metadata["teacher"]["model"] != "large_mlp"
            or metadata["teacher"]["seed"] != 0):
        raise ValueError("Synthetic teacher or endpoint provenance differs from protocol")
    for name in facs.SPLIT_NAMES:
        if metadata["inputs"][name]["sha256"] != base.FACS_SHA256[name]:
            raise ValueError("Synthetic data used a different FACS split")
    with np.load(path, allow_pickle=False) as arrays:
        ids = arrays["sequence_ids"]
        probabilities = arrays["probabilities"]
        if ids.ndim != 1 or ids.dtype.kind != "U" or len(set(ids)) != len(ids):
            raise ValueError("Synthetic IDs must be a unique Unicode vector")
        tokens = base.encode_tokens(ids)
        forbidden = {facs.WT_CDRH3, *primary_ids}
        for frame in splits.values():
            forbidden.update(frame["seq"])
        if set(ids) & forbidden or any("C" in seq for seq in ids):
            raise ValueError("Synthetic sequence exclusion failed")
        if (len(ids) != metadata["count"] or probabilities.shape != (len(ids), 3)
                or not np.isfinite(probabilities).all()
                or (probabilities < 0).any() or (probabilities > 1).any()
                or not np.allclose(probabilities.sum(1), 1, atol=1e-6, rtol=0)):
            raise ValueError("Invalid synthetic class probabilities")
    tensors = (base.torch.as_tensor(tokens, device=device),
               base.torch.as_tensor(probabilities, dtype=base.torch.float32, device=device))
    return tensors, {
        "dataset_metadata": base.file_source(directory / "metadata.json"),
        "dataset": base.file_source(path),
        "teacher": metadata["teacher"],
        "rows": len(ids),
        "weight": SYNTHETIC_WEIGHT,
        "label_type": "pseudolabels; no measured endpoints",
        "exact_overlap_with_real_or_primary": 0,
    }


def train_epoch(model, real, synthetic, optimizer, real_generator, synthetic_generator,
                batch_size=base.BATCH_SIZE):
    torch = base.torch
    model.train()
    order = torch.randperm(len(real[0]), device=real[0].device, generator=real_generator)
    synthetic_order = torch.randperm(len(synthetic[0]), device=real[0].device,
                                     generator=synthetic_generator)
    cursor = 0
    totals = np.zeros(2, dtype=float)
    for indices in order.split(batch_size):
        # Preserve every real batch and cycle independently through synthetic rows.
        if cursor == len(synthetic_order):
            synthetic_order = torch.randperm(len(synthetic[0]), device=real[0].device,
                                             generator=synthetic_generator)
            cursor = 0
        pseudo_indices = synthetic_order[cursor:cursor + batch_size]
        cursor += len(pseudo_indices)
        optimizer.zero_grad(set_to_none=True)
        real_loss = torch.nn.functional.cross_entropy(model(real[0][indices]), real[1][indices])
        pseudo_logits = model(synthetic[0][pseudo_indices])
        pseudo_loss = -(synthetic[1][pseudo_indices] * pseudo_logits.log_softmax(-1)).sum(-1).mean()
        loss = real_loss + SYNTHETIC_WEIGHT * pseudo_loss
        if not torch.isfinite(loss):
            raise ValueError("Nonfinite augmented loss")
        loss.backward()
        optimizer.step()
        totals += np.array([real_loss.item(), pseudo_loss.item()]) * len(indices)
    return dict(zip(("real_cross_entropy", "synthetic_cross_entropy"), totals / len(real[0])))


def fit(model, train, val, synthetic, seed, directory, log):
    torch = base.torch
    directory.mkdir()
    generator = torch.Generator(device=train[0].device).manual_seed(seed)
    pseudo_generator = torch.Generator(device=train[0].device).manual_seed(seed + 10000)
    optimizer = torch.optim.AdamW(model.parameters(), lr=base.LEARNING_RATE,
                                 weight_decay=base.WEIGHT_DECAY)
    val_high = val[1].cpu().numpy() == base.HIGH
    history, started = [], perf_counter()
    for epoch in range(1, 16):
        tick = perf_counter()
        losses = train_epoch(model, train, synthetic, optimizer, generator, pseudo_generator)
        scores, _ = base.predict(model, val[0])
        ap = float(average_precision_score(val_high, scores))
        if not math.isfinite(ap):
            raise ValueError("Nonfinite validation AP")
        checkpoint = base.save_checkpoint(directory / f"epoch-{epoch:02d}.pt", model)
        history.append({"epoch": epoch, **losses, "validation_high_ap": ap,
                        "checkpoint": checkpoint, "elapsed_seconds": perf_counter() - tick})
        selected = base.select_best_epoch(history)
        result = {"history": history, "best_epoch": selected["epoch"],
                  "best_validation_high_ap": selected["validation_high_ap"],
                  "best_checkpoint": selected["checkpoint"],
                  "parameter_count": sum(p.numel() for p in model.parameters()),
                  "latent_width": model.latent_width, "training_seconds": perf_counter() - started}
        base.write_json(directory / "history.json", result)
        log(f"augmented transformer epoch {epoch}/15 val_high_ap={ap:.6f} "
            f"best_epoch={selected['epoch']} elapsed_seconds={history[-1]['elapsed_seconds']:.2f}")
    selected_path = directory / result["best_checkpoint"]["file"]
    model.load_state_dict(torch.load(selected_path, map_location="cpu", weights_only=True))
    os.link(selected_path, directory / "best.pt")
    result["best_checkpoint"] = {**result["best_checkpoint"], "file": "best.pt"}
    base.write_json(directory / "history.json", result)
    return result


def run(args):
    if args.output.exists():
        raise FileExistsError(f"Output already exists: {args.output}")
    metadata = json.loads((args.baseline / "metadata.json").read_text())
    config = metadata["configuration"]
    if (metadata["status"] != "complete" or config["size"] != "large"
            or config["epochs"] != 15 or config["seed"] not in (0, 1, 2)
            or "transformer" not in config["models"] or metadata.get("augmentation")):
        raise ValueError("Expected completed, prespecified real-only large baseline")
    if metadata["code"]["trainer"]["sha256"] != base.file_source(Path(base.__file__))["sha256"]:
        raise ValueError("Baseline trainer source differs")
    seed = config["seed"]
    base.configure_torch(seed, args.device)
    splits, inputs, audits = base.load_facs_splits(args.data_dir)
    ids, inputs["frozen_oof"] = base.load_primary_ids(args.oof)
    audits["primary_spr_encoder_overlap"] = base.audit_primary_overlap(ids, splits)
    synthetic, augmentation = load_synthetic(args.synthetic, splits, ids, args.device)
    metadata = copy.deepcopy(metadata)
    metadata.pop("models")
    metadata.pop("output")
    metadata.update({"status": "training", "inputs": inputs, "audits": audits,
                     "augmentation": augmentation, "paired_baseline": base.file_source(args.baseline / "metadata.json")})
    metadata["configuration"].update({"output": str(args.output), "device": args.device,
                                     "models": ["transformer"], "loss": "real CE + 0.2 synthetic soft-target CE",
                                     "training_data": "Real FACS train plus separate teacher pseudolabels"})
    metadata["code"]["augmentation_trainer"] = base.file_source(Path(__file__))
    args.output.mkdir(parents=True)
    base.write_json(args.output / "run.json", metadata)

    def log(message):
        with (args.output / "training.log").open("a") as handle:
            handle.write(message + "\n")
        print(message, flush=True)

    tensors = {name: base.encode_split(frame, args.device) for name, frame in splits.items()}
    base.configure_torch(seed, args.device)
    model = base.build_model("transformer", "large").to(args.device)
    result = fit(model, tensors["train"], tensors["val"], synthetic, seed,
                 args.output / "transformer", log)
    model.requires_grad_(False).eval()
    scores, _ = base.predict(model, tensors["test"][0])
    result["facs_test"] = facs.classification_metrics(tensors["test"][1].cpu().numpy() == base.HIGH, scores)
    base.write_json(args.output / "transformer" / "metrics.json", result)
    tokens = base.torch.as_tensor(base.encode_tokens(ids), device=args.device)
    scores, features = base.predict(model, tokens, include_features=True)
    path = args.output / "features.npz"
    np.savez_compressed(path, sequence_ids=np.asarray(ids, dtype=str), facs_transformer=features)
    pd.DataFrame({"sequence_id": ids, "p_high_transformer": scores}).to_csv(
        args.output / "primary_probabilities.csv", index=False)
    metadata.update({"status": "complete", "models": {"transformer": result},
                     "output": {"file": path.name, "sha256": base.file_source(path)["sha256"]}})
    base.write_json(args.output / "metadata.json", metadata)
    base.write_json(args.output / "run.json", metadata)
    log(f"Complete: {args.output}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, required=True)
    parser.add_argument("--baseline", type=Path, required=True)
    parser.add_argument("--synthetic", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--oof", type=Path, default=base.DEFAULT_OOF)
    parser.add_argument("--device", choices=("cpu", "cuda"), default="cuda")
    run(parser.parse_args())


if __name__ == "__main__":
    main()
