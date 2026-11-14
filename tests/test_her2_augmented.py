"""Test the augmented HER2 trainer's provenance, loss and selection contracts."""

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from scripts import train_her2_augmented as augmented
from scripts import train_her2_neural as base
from src import her2_benchmark as facs


@pytest.fixture
def tensor_api(monkeypatch):
    """Keep provenance tests runnable when optional PyTorch is absent."""
    if base.torch is None:
        class NumpyTensorAPI:
            float32 = np.float32

            @staticmethod
            def as_tensor(values, *, dtype=None, device=None):
                return np.asarray(values, dtype=dtype)

        monkeypatch.setattr(base, "torch", NumpyTensorAPI)
    return base.torch


def _digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_synthetic(directory, *, ids=None, probabilities=None):
    directory.mkdir(parents=True, exist_ok=True)
    ids = ["FFFFFFFFFF", "GGGGGGGGGG"] if ids is None else ids
    probabilities = (
        np.asarray([[0.2, 0.3, 0.5]] * len(ids), dtype=np.float32)
        if probabilities is None else np.asarray(probabilities)
    )
    path = directory / "synthetic.npz"
    np.savez_compressed(path, sequence_ids=np.asarray(ids), probabilities=probabilities)
    metadata = {
        "label_type": "pseudolabels",
        "measured_endpoints": [],
        "class_order": list(base.CLASS_TO_INDEX),
        "teacher": {"model": "large_mlp", "seed": 0},
        "inputs": {
            name: {"sha256": base.FACS_SHA256[name]} for name in facs.SPLIT_NAMES
        },
        "count": len(ids),
        "output": {"sha256": _digest(path)},
    }
    (directory / "metadata.json").write_text(json.dumps(metadata), encoding="utf-8")
    return directory


@pytest.fixture
def real_context():
    splits = {
        "train": pd.DataFrame({"seq": ["AAAAAAAAAA"]}),
        "val": pd.DataFrame({"seq": ["DDDDDDDDDD"]}),
        "test": pd.DataFrame({"seq": ["EEEEEEEEEE"]}),
    }
    primary_ids = ["IIIIIIIIII"]
    return splits, primary_ids


def load(directory, context):
    splits, primary_ids = context
    return augmented.load_synthetic(directory, splits, primary_ids, "cpu")


def test_synthetic_source_round_trip_is_unicode_and_records_protocol(
    tmp_path, tensor_api, real_context
):
    directory = write_synthetic(tmp_path / "synthetic")

    (tokens, probabilities), audit = load(directory, real_context)

    assert len(tokens) == len(probabilities) == 2
    assert audit["rows"] == 2
    assert audit["weight"] == augmented.SYNTHETIC_WEIGHT == 0.2
    assert audit["teacher"] == {"model": "large_mlp", "seed": 0}
    assert audit["label_type"] == "pseudolabels; no measured endpoints"
    assert audit["exact_overlap_with_real_or_primary"] == 0


def test_changed_synthetic_file_is_rejected(tmp_path, tensor_api, real_context):
    directory = write_synthetic(tmp_path / "synthetic")
    path = directory / "synthetic.npz"
    path.write_bytes(path.read_bytes() + b"tampered")

    with pytest.raises(ValueError, match="Synthetic array fingerprint mismatch"):
        load(directory, real_context)


def test_synthetic_metadata_must_pin_all_facs_splits(tmp_path, tensor_api, real_context):
    directory = write_synthetic(tmp_path / "synthetic")
    metadata_path = directory / "metadata.json"
    metadata = json.loads(metadata_path.read_text())
    metadata["inputs"]["val"]["sha256"] = "0" * 64
    metadata_path.write_text(json.dumps(metadata), encoding="utf-8")

    with pytest.raises(ValueError, match="different FACS split"):
        load(directory, real_context)


@pytest.mark.parametrize("field,value", [
    ("class_order", ["high", "mid", "low"]),
    ("teacher", {"model": "other_teacher", "seed": 0}),
    ("teacher", {"model": "large_mlp", "seed": 17}),
    ("label_type", "measured_labels"),
    ("measured_endpoints", ["SPR"]),
])
def test_protocol_mismatch_is_rejected(tmp_path, tensor_api, real_context, field, value):
    directory = write_synthetic(tmp_path / "synthetic")
    metadata_path = directory / "metadata.json"
    metadata = json.loads(metadata_path.read_text())
    metadata[field] = value
    metadata_path.write_text(json.dumps(metadata), encoding="utf-8")

    with pytest.raises(ValueError, match="Synthetic teacher or endpoint provenance"):
        load(directory, real_context)


@pytest.mark.parametrize("sequence", [
    "AAAAAAAAAA",  # FACS train
    "DDDDDDDDDD",  # FACS validation
    "EEEEEEEEEE",  # FACS test
    "IIIIIIIIII",  # frozen primary cohort
    facs.WT_CDRH3,
])
def test_exact_overlap_with_each_protected_source_is_rejected(
    tmp_path, tensor_api, real_context, sequence
):
    directory = write_synthetic(tmp_path / "synthetic", ids=[sequence])

    with pytest.raises(ValueError, match="Synthetic sequence exclusion failed"):
        load(directory, real_context)


def test_cysteine_sequence_is_rejected(tmp_path, tensor_api, real_context):
    directory = write_synthetic(tmp_path / "synthetic", ids=["CCCCCCCCCC"])

    with pytest.raises(ValueError, match="Synthetic sequence exclusion failed"):
        load(directory, real_context)


@pytest.mark.parametrize("ids,probabilities", [
    (["FFFFFFFFFF", "FFFFFFFFFF"], [[0.2, 0.3, 0.5]] * 2),
    ([b"FFFFFFFFFF"], [[0.2, 0.3, 0.5]]),
    (["FFFFFFFFFF"], [[0.2, 0.3]]),
    (["FFFFFFFFFF"], [[0.2, float("nan"), 0.8]]),
    (["FFFFFFFFFF"], [[-0.1, 0.3, 0.8]]),
    (["FFFFFFFFFF"], [[0.2, 0.3, 0.4]]),
])
def test_malformed_synthetic_arrays_are_rejected(
    tmp_path, tensor_api, real_context, ids, probabilities
):
    directory = write_synthetic(
        tmp_path / "synthetic", ids=ids, probabilities=probabilities
    )
    metadata_path = directory / "metadata.json"
    metadata = json.loads(metadata_path.read_text())
    metadata["output"]["sha256"] = _digest(directory / "synthetic.npz")
    metadata["count"] = len(ids)
    metadata_path.write_text(json.dumps(metadata), encoding="utf-8")

    with pytest.raises(ValueError):
        load(directory, real_context)


def _torch_or_skip():
    return pytest.importorskip("torch")


def test_augmented_real_batches_match_baseline_order_and_cover_each_row_once():
    torch = _torch_or_skip()

    class RecordingModel(torch.nn.Module):
        def __init__(self):
            super().__init__()
            self.linear = torch.nn.Linear(1, 3)
            self.seen = []

        def forward(self, values):
            self.seen.append(values.detach().clone())
            return self.linear(values.to(torch.float32))

    class NoOpOptimizer:
        def __init__(self, parameters):
            self.parameters = list(parameters)

        def zero_grad(self, set_to_none=True):
            for parameter in self.parameters:
                parameter.grad = None

        def step(self):
            pass

    real = (
        torch.arange(7, dtype=torch.int64).reshape(-1, 1),
        torch.tensor([0, 1, 2, 0, 1, 2, 0]),
    )
    synthetic = (
        torch.arange(100, 104, dtype=torch.int64).reshape(-1, 1),
        torch.tensor([[0.2, 0.3, 0.5]] * 4, dtype=torch.float32),
    )
    seed = 23

    baseline_model = RecordingModel()
    baseline_model.load_state_dict(RecordingModel().state_dict())
    base.train_one_epoch(
        baseline_model, *real, NoOpOptimizer(baseline_model.parameters()), batch_size=3,
        generator=torch.Generator(device="cpu").manual_seed(seed),
    )
    baseline_order = torch.cat(baseline_model.seen).flatten()

    augmented_model = RecordingModel()
    augmented_model.load_state_dict(baseline_model.state_dict())
    augmented.train_epoch(
        augmented_model, real, synthetic, NoOpOptimizer(augmented_model.parameters()),
        torch.Generator(device="cpu").manual_seed(seed),
        torch.Generator(device="cpu").manual_seed(seed + 10000), batch_size=3,
    )
    augmented_real_order = torch.cat([
        batch for batch in augmented_model.seen if bool((batch < 50).all())
    ]).flatten()

    assert torch.equal(augmented_real_order, baseline_order)
    assert sorted(augmented_real_order.tolist()) == list(range(len(real[0])))


def test_soft_target_loss_is_finite_and_contributes_weighted_gradients():
    torch = _torch_or_skip()

    class TinyClassifier(torch.nn.Module):
        def __init__(self):
            super().__init__()
            self.linear = torch.nn.Linear(1, 3)

        def forward(self, values):
            return self.linear(values.to(torch.float32))

    class CaptureOptimizer:
        def __init__(self, parameters):
            self.parameters = list(parameters)
            self.gradients = None

        def zero_grad(self, set_to_none=True):
            for parameter in self.parameters:
                parameter.grad = None

        def step(self):
            self.gradients = [parameter.grad.detach().clone() for parameter in self.parameters]

    model = TinyClassifier()
    with torch.no_grad():
        model.linear.weight.copy_(torch.tensor([[0.2], [-0.4], [0.7]]))
        model.linear.bias.copy_(torch.tensor([0.1, 0.3, -0.2]))
    real = (torch.tensor([[1.2]]), torch.tensor([0]))
    synthetic = (torch.tensor([[2.4]]), torch.tensor([[0.1, 0.7, 0.2]]))

    reference = TinyClassifier()
    reference.load_state_dict(model.state_dict())
    real_loss = torch.nn.functional.cross_entropy(reference(real[0]), real[1])
    pseudo_loss = -(
        synthetic[1] * reference(synthetic[0]).log_softmax(dim=-1)
    ).sum(dim=-1).mean()
    (real_loss + augmented.SYNTHETIC_WEIGHT * pseudo_loss).backward()
    expected_gradients = [parameter.grad.detach().clone() for parameter in reference.parameters()]

    real_only = TinyClassifier()
    real_only.load_state_dict(model.state_dict())
    torch.nn.functional.cross_entropy(real_only(real[0]), real[1]).backward()
    real_only_gradients = [parameter.grad.detach().clone() for parameter in real_only.parameters()]

    optimizer = CaptureOptimizer(model.parameters())
    totals = augmented.train_epoch(
        model, real, synthetic, optimizer,
        torch.Generator(device="cpu").manual_seed(0),
        torch.Generator(device="cpu").manual_seed(10000), batch_size=1,
    )

    assert np.isfinite(totals["real_cross_entropy"])
    assert np.isfinite(totals["synthetic_cross_entropy"])
    for actual, expected in zip(optimizer.gradients, expected_gradients, strict=True):
        torch.testing.assert_close(actual, expected)
    assert any(
        not torch.allclose(actual, real_grad)
        for actual, real_grad in zip(optimizer.gradients, real_only_gradients, strict=True)
    )


def test_fit_runs_fixed_epochs_and_restores_validation_selected_checkpoint(monkeypatch, tmp_path):
    torch = _torch_or_skip()

    class TinyClassifier(torch.nn.Module):
        def __init__(self):
            super().__init__()
            self.latent_width = 1
            self.linear = torch.nn.Linear(1, 3)

        def forward(self, values):
            return self.linear(values.to(torch.float32))

    train = (torch.tensor([[0.2], [0.8]]), torch.tensor([2, 0]))
    val = (torch.tensor([[0.3], [0.9]]), torch.tensor([2, 0]))
    synthetic = (
        torch.tensor([[1.3], [1.7]]),
        torch.tensor([[0.1, 0.2, 0.7], [0.6, 0.3, 0.1]]),
    )
    calls = 0

    def controlled_predict(_model, _tokens):
        nonlocal calls
        calls += 1
        # Epoch 2 has the best validation AP; test rows and training loss do not select it.
        return (np.array([0.9, 0.1]) if calls == 2 else np.array([0.1, 0.9])), None

    monkeypatch.setattr(base, "predict", controlled_predict)
    torch.manual_seed(7)
    model = TinyClassifier()
    output = tmp_path / "transformer"
    result = augmented.fit(model, train, val, synthetic, 5, output, lambda _: None)

    selected = torch.load(output / "epoch-02.pt", map_location="cpu", weights_only=True)
    last = torch.load(output / "epoch-15.pt", map_location="cpu", weights_only=True)
    best_link = torch.load(output / "best.pt", map_location="cpu", weights_only=True)
    assert len(result["history"]) == calls == 15
    assert result["best_epoch"] == 2
    assert result["best_validation_high_ap"] == pytest.approx(1.0)
    assert result["best_checkpoint"]["file"] == "best.pt"
    assert all(torch.equal(model.state_dict()[key], selected[key]) for key in selected)
    assert all(torch.equal(best_link[key], selected[key]) for key in selected)
    assert any(not torch.equal(selected[key], last[key]) for key in selected)
    assert result["history"][1]["validation_high_ap"] > result["history"][0]["validation_high_ap"]
