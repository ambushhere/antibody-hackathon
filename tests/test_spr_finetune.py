"""Check gradient boundaries and grouped SPR fit isolation on CPU."""

import json

import numpy as np
import pandas as pd
import pytest

from scripts import benchmark_spr_finetune as experiment
from scripts import train_her2_neural as neural
from src import her2_spr_benchmark as spr


@pytest.fixture
def torch_cpu():
    torch = pytest.importorskip("torch")
    original_threads = torch.get_num_threads()
    deterministic = torch.are_deterministic_algorithms_enabled()
    neural.configure_torch(0, "cpu")
    yield torch
    torch.set_num_threads(original_threads)
    torch.use_deterministic_algorithms(deterministic)


@pytest.fixture
def tiny_data(torch_cpu, monkeypatch):
    monkeypatch.setitem(neural.ARCHITECTURES["large"], "transformer", (16, 2, 4, 32))
    ids = [aa * 9 + tail for aa in "ACDEFGHIKL" for tail in "AC"]
    raw_tokens = neural.encode_tokens(ids)
    tokens = torch_cpu.as_tensor(raw_tokens)
    state = neural.build_model("transformer", "large").state_dict()
    groups = np.repeat(np.arange(10), 2)
    labels = np.linspace(-3.0, 0.0, num=len(ids))
    folds = [
        (fold, train, test, spr.ALPHAS[0])
        for fold, (train, test) in enumerate(spr.group_splits(groups), start=1)
    ]
    inner = [spr.group_splits(groups[train]) for _, train, _, _ in folds]
    return ids, tokens, labels, groups, folds, inner, state


def test_both_arms_start_equal_and_only_intended_parameters_can_update(torch_cpu, tiny_data):
    _, tokens, _, _, _, _, state = tiny_data
    frozen, _ = experiment.build_regressor(state, experiment.FROZEN_HEAD, "cpu")
    adapted, optimizer = experiment.build_regressor(state, experiment.FINAL_BLOCK, "cpu")
    assert not any(module.training for module in frozen.modules())
    assert not any(module.training for module in adapted.modules())
    torch_cpu.testing.assert_close(frozen.head.weight, adapted.head.weight, rtol=0, atol=0)
    torch_cpu.testing.assert_close(frozen.head.bias, adapted.head.bias, rtol=0, atol=0)
    assert frozen.head.out_features == adapted.head.out_features == 1
    torch_cpu.testing.assert_close(frozen(tokens), adapted(tokens), rtol=1e-5, atol=1e-6)
    assert [group["lr"] for group in optimizer.param_groups] == [1e-3, 1e-5]
    assert [group["weight_decay"] for group in optimizer.param_groups] == [0.01, 0.01]

    assert {name for name, param in frozen.named_parameters() if param.requires_grad} == {
        "head.weight", "head.bias"
    }
    initial = {name: param.detach().clone() for name, param in adapted.named_parameters()}
    adapted(tokens).square().mean().backward()
    optimizer.step()
    for name, param in adapted.named_parameters():
        trainable = name.startswith(("head.", "encoder.1.", "norm."))
        assert param.requires_grad is trainable
        assert (param.grad is not None) is trainable
        if not trainable:
            assert torch_cpu.equal(initial[name], param)
    assert not torch_cpu.equal(initial["encoder.1.linear1.weight"], adapted.encoder[1].linear1.weight)
    assert not torch_cpu.equal(initial["norm.weight"], adapted.norm.weight)


def test_outer_labels_never_enter_selection_or_fit(tiny_data, monkeypatch):
    _, tokens, labels, groups, folds, inner, state = tiny_data
    monkeypatch.setattr(experiment, "EPOCHS", (1, 2))
    fold, train, test, alpha = folds[0]
    # Put this fold first in both runs. Perturbing its held-out labels may affect other folds.
    original, records = experiment.nested_oof(tokens, labels, groups, folds, inner, state)
    changed_labels = labels.copy()
    changed_labels[test] += 10_000
    changed, changed_records = experiment.nested_oof(tokens, changed_labels, groups, folds, inner, state)

    record = records[0]
    assert record["fold"] == fold
    assert not set(record["train_groups"]) & set(record["test_groups"])
    for split in record["inner_folds"]:
        assert set(split["train_groups"]) | set(split["val_groups"]) == set(record["train_groups"])
        assert not set(split["train_groups"]) & set(split["val_groups"])
    for arm in experiment.ARMS:
        assert record["arms"][arm] == changed_records[0]["arms"][arm]
        np.testing.assert_array_equal(original[arm][test], changed[arm][test])
        outer_transform = record["arms"][arm]["outer_fit"]["target_transform"]
        assert outer_transform == {"mean": labels[train].mean(), "std": labels[train].std(ddof=0)}
        for (fit, _), fit_record in zip(inner[0], record["arms"][arm]["inner_fits"]):
            train_targets = labels[train[fit]]
            assert fit_record["target_transform"] == {
                "mean": train_targets.mean(), "std": train_targets.std(ddof=0),
            }


def test_normalization_rejects_constant_training_targets():
    normalized, transform = experiment.normalize_targets([1.0, 5.0])
    np.testing.assert_array_equal(normalized, [-1.0, 1.0])
    assert transform == {"mean": 3.0, "std": 2.0}
    with pytest.raises(ValueError, match="standard deviation"):
        experiment.normalize_targets([1.0, 1.0])


def test_cpu_smoke_cli(tiny_data, torch_cpu, tmp_path, monkeypatch):
    ids, _, labels, groups, folds, inner, state = tiny_data
    checkpoint = tmp_path / "initial.pt"
    torch_cpu.save(state, checkpoint)
    output = tmp_path / "experiment"
    cli = ["--data-dir", str(tmp_path), "--checkpoint", str(checkpoint),
           "--output", str(output), "--device", "cpu"]
    with pytest.raises(SystemExit, match="2"):
        experiment.main(cli)
    assert not output.exists()

    fingerprint = neural.file_source(checkpoint)["sha256"]
    monkeypatch.setattr(experiment, "CHECKPOINT_SHA256", fingerprint)
    monkeypatch.setattr(experiment, "EPOCHS", (1, 2, 3))
    monkeypatch.setattr(spr, "BOOTSTRAP_DRAWS", 5)
    rows = pd.DataFrame({
        "cdrh3": ids, spr.TARGET: labels, "group": groups, experiment.interaction.DISTANT: groups < 5,
    })
    for fold, _, test, _ in folds:
        rows.loc[test, "fold"] = fold
    baselines = {experiment.BASELINES[0]: labels[::-1], experiment.BASELINES[1]: labels * 0.5}
    monkeypatch.setattr(experiment, "load_cohort", lambda _: (rows, folds, inner, baselines, {}))
    experiment.main(cli)

    metrics_path = output / "metrics.json"
    metrics = json.loads(metrics_path.read_text())
    assert metrics["inputs"]["checkpoint"]["sha256"] == fingerprint
    assert len(metrics["folds"]) == 5
    assert metrics["subsets"]["pooled"]["n"] == 20
    assert metrics["subsets"]["distant75"]["n"] == 10
    assert "final_block minus frozen_head" in metrics["subsets"]["pooled"]["paired_differences"]
    predictions = pd.read_csv(output / "oof_predictions.csv")
    assert predictions["cdrh3"].tolist() == ids
    assert np.isfinite(predictions[[*experiment.ARMS, *experiment.BASELINES]].to_numpy()).all()
    assert {path.name for path in output.iterdir()} == {"metrics.json", "report.md", "oof_predictions.csv"}
    with pytest.raises(SystemExit, match="2"):
        experiment.main(cli)
    assert json.loads(metrics_path.read_text()) == metrics
