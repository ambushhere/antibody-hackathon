"""Compare frozen and final-block adaptation of the seed 0 FACS Transformer on SPR."""

import argparse
import hashlib
import io
import platform
import sys
from importlib.metadata import version
from pathlib import Path

import numpy as np

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from scripts import rank_spr_candidates as frozen
from scripts import train_her2_neural as neural
from src import her2_benchmark as facs
from src import her2_interaction_benchmark as interaction
from src import her2_spr_benchmark as spr

ROOT = Path(__file__).resolve().parents[1]
METRICS_PATH = ROOT / "results/spr-benchmark/metrics.json"
DISTANCE_PATH = ROOT / "results/spr-distance-diagnostic/distance_rows.csv"
BASELINE_PATH = ROOT / "results/neural-all-benchmark/oof_predictions.csv"
PROTOCOL_PATH = ROOT / "docs/research/spr-finetune-protocol.md"
METRICS_SHA256 = "9668a01ca8dcfdcc8e84b127c4d5f8de1580232090f5efd2925f608ae3c41edf"
BASELINE_SHA256 = "95e8dfbbfeb2a7d6daef8cc4582252ef3a341c29f1869eb4503977ec8140150e"
CHECKPOINT_SHA256 = "54c6dd9de3ded6f1112e7c1337b89b48be21c5348a97d6e96124e866bcafdd69"
SEED = 0
EPOCHS = (20, 50, 100)
HEAD_LR = 1e-3
BLOCK_LR = 1e-5
WEIGHT_DECAY = 0.01
FROZEN_HEAD = "frozen_head"
FINAL_BLOCK = "final_block"
ARMS = (FROZEN_HEAD, FINAL_BLOCK)
BASELINES = (interaction.POSITION, interaction.PAIRWISE)
COMPARISONS = [(FINAL_BLOCK, FROZEN_HEAD), *((arm, base) for arm in ARMS for base in BASELINES)]


def load_cohort(data_dir):
    metrics_source = neural.file_source(METRICS_PATH)
    neural.require_hash(metrics_source, METRICS_SHA256)
    benchmark, cohort = frozen.load_benchmark(METRICS_PATH)
    primary, inputs, _ = frozen.load_primary(data_dir, benchmark)
    oof, inputs["frozen_oof"] = frozen.load_oof(neural.DEFAULT_OOF, (interaction.POSITION, facs.SCORE))
    neural.require_hash(inputs["frozen_oof"], neural.OOF_SHA256)
    folds, aligned = frozen.original_folds(primary, oof, cohort)
    frozen.max_gap(primary[spr.TARGET], aligned[spr.TARGET], "Original OOF labels", frozen.OOF_TOLERANCE)
    distances, inputs["distance_diagnostic"] = interaction.load_distances(DISTANCE_PATH, primary, aligned)
    inputs["frozen_metrics"] = metrics_source

    columns = ("cdrh3", "fold", "group", spr.TARGET, *distances.columns, *BASELINES)
    baseline, inputs["ridge_oof"] = facs.read_csv(BASELINE_PATH, columns=columns)
    neural.require_hash(inputs["ridge_oof"], BASELINE_SHA256)
    baseline = baseline.set_index("cdrh3")
    if not baseline.index.is_unique or set(baseline.index) != set(aligned.index):
        raise ValueError("Ridge OOF sequence IDs differ from the frozen primary cohort")
    baseline = baseline.loc[aligned.index]
    if not np.array_equal(baseline[interaction.DISTANT].to_numpy(), distances[interaction.DISTANT].map(str)):
        raise ValueError("Ridge OOF distant subset differs from the frozen distance rows")
    for name in ("group", "fold", spr.TARGET, interaction.POSITION):
        frozen.max_gap(aligned[name], baseline[name], f"Ridge OOF {name}", frozen.OOF_TOLERANCE)
    for name in ("nearest_facs_trainval", "nearest_spr_outer_train"):
        frozen.max_gap(distances[name], baseline[name], f"Ridge OOF {name}", 0)
    scores = {name: frozen.finite_array(baseline[name], f"Ridge OOF {name}") for name in BASELINES}

    groups = primary["group"].to_numpy()
    inner_folds = []
    for (fold, train, _, _), reference in zip(folds, cohort["folds"]):
        inner = spr.group_splits(groups[train])
        membership = [
            {"train_n": len(fit), "val_n": len(val), "val_groups": np.unique(groups[train[val]]).tolist()}
            for fit, val in inner
        ]
        if membership != reference["inner_folds"]:
            raise ValueError(f"Inner groups differ from the frozen benchmark in outer fold {fold}")
        inner_folds.append(inner)
    return primary.join(distances).assign(fold=aligned["fold"].to_numpy()), folds, inner_folds, scores, inputs


def load_checkpoint(path):
    neural.require_torch()
    data = path.read_bytes()
    source = {"source": str(path), "sha256": hashlib.sha256(data).hexdigest()}
    # ASVS 1.5.2: verify the archived bytes before deserializing tensors on CPU.
    neural.require_hash(source, CHECKPOINT_SHA256)
    state = neural.torch.load(io.BytesIO(data), map_location="cpu", weights_only=True)
    return state, source


def build_regressor(state, arm, device):
    torch = neural.torch
    torch.manual_seed(SEED)
    model = neural.build_model("transformer", "large")
    model.load_state_dict(state, strict=True)
    model.requires_grad_(False)
    model.head = neural.nn.Linear(model.latent_width, 1)
    head = list(model.head.parameters())
    parameters = [{"params": head, "lr": HEAD_LR}]
    if arm == FINAL_BLOCK:
        model.encoder[-1].requires_grad_(True)
        model.norm.requires_grad_(True)
        parameters.append({
            "params": [*model.encoder[-1].parameters(), *model.norm.parameters()],
            "lr": BLOCK_LR,
        })
    # Eval mode disables dropout without disabling gradients in the unfrozen block.
    model.to(device).eval()
    optimizer = torch.optim.AdamW(parameters, weight_decay=WEIGHT_DECAY)
    return model, optimizer


def normalize_targets(train_targets):
    targets = frozen.finite_array(train_targets, "SPR training targets")
    mean, std = float(targets.mean()), float(targets.std(ddof=0))
    if not np.isfinite(std) or std <= 0:
        raise ValueError("SPR training targets must have positive finite standard deviation")
    return (targets - mean) / std, {"mean": mean, "std": std}


def fit_predict(state, arm, train_tokens, train_targets, predict_tokens, epochs):
    """Fit from scratch on training targets and predict at the requested epochs."""
    torch = neural.torch
    device = train_tokens.device
    normalized, transform = normalize_targets(train_targets)
    targets = torch.as_tensor(normalized, dtype=torch.float32, device=device)
    model, optimizer = build_regressor(state, arm, device)
    if arm == FROZEN_HEAD:
        with torch.no_grad():
            train_tokens, predict_tokens = model.encode(train_tokens), model.encode(predict_tokens)
        forward = model.head
    else:
        forward = model

    predictions = {}
    for epoch in range(1, max(epochs) + 1):
        optimizer.zero_grad(set_to_none=True)
        loss = neural.F.mse_loss(forward(train_tokens).squeeze(-1), targets)
        loss.backward()
        optimizer.step()
        if epoch in epochs:
            with torch.no_grad():
                values = forward(predict_tokens).squeeze(-1).cpu().numpy().astype(float)
            predictions[epoch] = frozen.finite_array(
                values * transform["std"] + transform["mean"], f"{arm} epoch {epoch} predictions"
            )
    return predictions, {
        "target_transform": transform,
        "trainable_parameters": [name for name, param in model.named_parameters() if param.requires_grad],
    }


def nested_oof(tokens, labels, groups, folds, inner_folds, state, log=print):
    scores = {arm: np.full(len(labels), np.nan) for arm in ARMS}
    records = []
    for (fold, train, test, _), inner in zip(folds, inner_folds):
        record = {
            "fold": fold,
            "train_groups": np.unique(groups[train]).tolist(),
            "test_groups": np.unique(groups[test]).tolist(),
            "inner_folds": [
                {
                    "train_groups": np.unique(groups[train[fit]]).tolist(),
                    "val_groups": np.unique(groups[train[val]]).tolist(),
                }
                for fit, val in inner
            ],
            "arms": {},
        }
        for arm in ARMS:
            validation = {epoch: np.full(len(train), np.nan) for epoch in EPOCHS}
            inner_fits = []
            for fit, val in inner:
                fit_rows = train[fit]
                predicted, fit_record = fit_predict(
                    state, arm, tokens[fit_rows], labels[fit_rows], tokens[train[val]], EPOCHS
                )
                for epoch in EPOCHS:
                    validation[epoch][val] = predicted[epoch]
                inner_fits.append(fit_record)
            errors = {
                epoch: spr.rmse(labels[train], frozen.finite_array(values, f"{arm} pooled validation"))
                for epoch, values in validation.items()
            }
            chosen = min(EPOCHS, key=lambda epoch: (errors[epoch], epoch))
            predicted, outer_fit = fit_predict(
                state, arm, tokens[train], labels[train], tokens[test], (chosen,)
            )
            scores[arm][test] = predicted[chosen]
            record["arms"][arm] = {
                "chosen_epoch": chosen, "inner_pooled_rmse": errors, "inner_fits": inner_fits,
                "outer_fit": outer_fit,
            }
            log(f"Outer fold {fold}, {arm}: selected epoch {chosen}")
        records.append(record)
    for name, values in scores.items():
        frozen.finite_array(values, f"{name} OOF predictions")
    return scores, records


def render_report(result):
    def display(metric):
        estimate = metric["estimate"]
        interval = metric["95_percent"]
        if estimate is None:
            return "undefined"
        if interval is None:
            return f"{estimate:.4f}"
        return f"{estimate:.4f} [{interval[0]:.4f}, {interval[1]:.4f}]"

    lines = [
        "# SPR adaptation benchmark", "",
        "Retrospective evaluation on the previously examined cohort with reportable KD. "
        "The target is −log10(KD in nM). Cohort reuse provides no independent replication "
        "or basis for candidate promotion.", "",
        "Brackets give 95% paired sequence group bootstrap intervals on fixed OOF predictions. "
        "A positive Spearman difference or negative RMSE difference favors the first model.", "",
    ]
    for name, subset in result["subsets"].items():
        label = "Pooled" if name == "pooled" else "Distant subset"
        lines.extend([
            f"## {label} ({subset['n']} variants, {subset['groups']} groups)", "",
            "| Model or paired difference | Spearman | RMSE |",
            "| --- | ---: | ---: |",
        ])
        for entries in (subset["metrics"], subset["paired_differences"]):
            lines.extend(
                f"| {model} | {display(metrics['spearman'])} | {display(metrics['rmse'])} |"
                for model, metrics in entries.items()
            )
        lines.append("")
    return "\n".join(lines)


def run(data_dir, checkpoint, output, device):
    if output.exists() or output.is_symlink():
        raise FileExistsError(f"Output already exists: {output}")
    state, weights = load_checkpoint(checkpoint)
    rows, folds, inner, scores, inputs = load_cohort(data_dir)
    inputs["checkpoint"] = weights
    code_paths = [
        Path(__file__).resolve(), PROTOCOL_PATH,
        *(Path(module.__file__).resolve() for module in (neural, frozen, interaction, facs, spr)),
    ]
    code = {str(path.relative_to(ROOT)): neural.file_source(path) for path in code_paths}
    neural.configure_torch(SEED, device)
    torch = neural.torch
    tokens = torch.as_tensor(neural.encode_tokens(rows["cdrh3"]), device=device)
    labels, groups = rows[spr.TARGET].to_numpy(), rows["group"].to_numpy()

    # Reserve a fresh directory atomically before any output is written.
    output.mkdir(parents=True, exist_ok=False)
    adapted, records = nested_oof(tokens, labels, groups, folds, inner, state,
                                  log=lambda message: print(message, flush=True))
    scores.update(adapted)
    mask = rows[interaction.DISTANT].to_numpy()
    result = {
        "status": "retrospective_exploratory_reused_cohort",
        "target": spr.TARGET,
        "cohort": frozen.COHORT,
        "configuration": {
            "seed": SEED, "architecture": neural.ARCHITECTURES["large"]["transformer"],
            "architecture_fields": ["width", "blocks", "heads", "feedforward_width"],
            "residue_alphabet": facs.AMINO_ACIDS, "device": device, "dtype": "float32",
            "device_name": torch.cuda.get_device_name() if device == "cuda" else platform.machine(),
            "deterministic_algorithms": torch.are_deterministic_algorithms_enabled(),
            "encoder_mode": "eval in both arms, including final block and norm",
            "dropout_enabled": False,
            "batch": "full training split", "optimizer": "AdamW",
            "head_lr": HEAD_LR, "final_block_and_norm_lr": BLOCK_LR, "weight_decay": WEIGHT_DECAY,
            "loss": "MSE on train mean and population std normalized target",
            "epochs": EPOCHS, "selection": "Pooled inner RMSE on original scale, ties prefer earlier epoch",
        },
        "bootstrap": {
            "draws": spr.BOOTSTRAP_DRAWS, "seed": spr.SEED,
            "method": "Paired percentile intervals from resampling original sequence groups within each subset",
        },
        "folds": records,
        "subsets": {
            "pooled": interaction.summarize(labels, scores, groups, COMPARISONS),
            interaction.DISTANT: interaction.summarize(
                labels[mask], {name: values[mask] for name, values in scores.items()}, groups[mask], COMPARISONS
            ),
        },
        "inputs": inputs, "code": code,
        "versions": {
            "python": platform.python_version(),
            **{name: version(name) for name in ("torch", "numpy", "pandas", "scipy", "scikit-learn")},
            "torch_cuda": torch.version.cuda,
        },
    }
    predictions_path = output / "oof_predictions.csv"
    # ASVS 1.2.10: output text is canonical residue IDs and fixed headers. All scores are numeric.
    with neural.staged_file(predictions_path) as temporary:
        rows.assign(**scores).to_csv(temporary, index=False)
    result["oof_predictions"] = neural.file_source(predictions_path)
    with neural.staged_file(output / "report.md") as temporary:
        temporary.write_text(render_report(result), encoding="utf-8")
    neural.write_json(output / "metrics.json", result)
    return result


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, required=True)
    parser.add_argument("--checkpoint", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--device", choices=("cpu", "cuda"), default="cuda")
    args = parser.parse_args(argv)
    try:
        run(args.data_dir, args.checkpoint, args.output, args.device)
    except (ValueError, OSError) as exc:
        parser.exit(2, f"error: {exc}\n")


if __name__ == "__main__":
    main()
