"""Exercise FACS-only checkpoint selection and the frozen SPR feature contract."""

import hashlib
import json
from zipfile import ZipFile

import numpy as np
import pandas as pd
import pytest

from scripts import train_her2_neural as neural
from src import her2_benchmark as facs


@pytest.fixture
def torch_cpu():
    torch = pytest.importorskip("torch")
    original_threads = torch.get_num_threads()
    original_deterministic = torch.are_deterministic_algorithms_enabled()
    neural.configure_torch(0, "cpu")
    yield torch
    torch.set_num_threads(original_threads)
    torch.use_deterministic_algorithms(original_deterministic)


def toy_frames():
    splits = {}
    for name, residues in zip(facs.SPLIT_NAMES, ("ACDEFG", "HIK", "LMN")):
        sequences = [aa * facs.LOOP_LENGTH for aa in residues]
        classes = (["low", "mid", "high"] * 2)[:len(sequences)]
        splits[name] = pd.DataFrame({
            "seq": sequences,
            "label": [int(label == "high") for label in classes],
            "class": classes,
            "cluster_centre": sequences,
        })
    return splits


def write_toy_data(directory, splits):
    archive = directory / neural.SPLIT_ARCHIVE
    archive.parent.mkdir(parents=True)
    hashes = {}
    with ZipFile(archive, "w") as output:
        for name, split in splits.items():
            data = split.to_csv(index=False).encode()
            output.writestr(f"{facs.SPLIT_PATH}{name}.csv", data)
            hashes[name] = hashlib.sha256(data).hexdigest()
    mapping = directory / neural.FACS_ARCHIVE
    mapping.parent.mkdir(parents=True)
    unique = pd.concat(splits.values()).drop_duplicates("seq")
    csv = unique[["seq", "class"]].rename(columns={"seq": "cdrh3", "class": "binding_class"})
    data = csv.to_csv(index=False).encode()
    with ZipFile(mapping, "w") as output:
        output.writestr(neural.FACS_CSV, data)
    hashes["facs"] = hashlib.sha256(data).hexdigest()
    return hashes


def test_position_tokens_and_one_hot_match_facs(torch_cpu):
    seqs = [" ACDEFGHIKL ", "CAAAAAAAAA"]
    tokens = neural.encode_tokens(seqs)
    np.testing.assert_array_equal(tokens, [list(range(10)), [1] + [0] * 9])
    assert tokens.dtype == np.int64
    one_hot = neural.position_one_hot(torch_cpu.as_tensor(tokens))
    assert one_hot.dtype == torch_cpu.float32
    np.testing.assert_array_equal(one_hot.numpy(), facs.sequence_features(pd.Series(seqs), "position").toarray())
    with pytest.raises(ValueError, match="unknown residues"):
        neural.encode_tokens(["ACDEFGHIKX"])


@pytest.mark.parametrize("name", neural.MODELS)
@pytest.mark.parametrize("size,width", [("small", 64), ("large", 256)])
def test_model_shapes_and_eval_repeatability(torch_cpu, name, size, width):
    tokens = torch_cpu.as_tensor(neural.encode_tokens(["ACDEFGHIKL", "ACDEFGHIKY"]))
    model = neural.build_model(name, size).eval()
    with torch_cpu.no_grad():
        assert model(tokens).shape == (2, 3)
    scores, features = neural.predict(model, tokens, include_features=True)
    repeated_scores, repeated_features = neural.predict(model, tokens, include_features=True)
    assert features.shape == (2, width) and features.dtype == np.float32
    np.testing.assert_array_equal(features, repeated_features)
    np.testing.assert_array_equal(scores, repeated_scores)
    assert np.isfinite(features).all()
    assert np.all((0 <= scores) & (scores <= 1))
    if name == "transformer":
        # CLS occupies index 0. Causal masking would hide even the final residue from CLS.
        assert not np.allclose(features[0], features[1])
        assert model.positions.shape == (1, 11, width)
        assert all(layer.norm_first for layer in model.encoder)


@pytest.mark.parametrize("name", neural.MODELS)
def test_train_step_updates_weights(torch_cpu, name):
    tokens, labels = neural.encode_split(toy_frames()["train"], "cpu")
    model = neural.build_model(name)
    optimizer = torch_cpu.optim.AdamW(model.parameters(), lr=neural.LEARNING_RATE,
                                     weight_decay=neural.WEIGHT_DECAY)
    before = model.head.weight.detach().clone()
    loss = neural.train_one_epoch(
        model, tokens, labels, optimizer, batch_size=3,
        generator=torch_cpu.Generator().manual_seed(0),
    )
    assert np.isfinite(loss)
    assert not torch_cpu.equal(model.head.weight, before)


def test_checkpoint_rule_uses_validation_and_earliest_tie():
    history = [
        {"epoch": 4, "validation_high_ap": 0.8, "test_high_ap": 1.0},
        {"epoch": 2, "validation_high_ap": 0.8},
        {"epoch": 1, "validation_high_ap": 0.7, "test_high_ap": 1.0},
        {"epoch": 3, "validation_high_ap": 0.6, "test_high_ap": 1.0},
    ]
    assert neural.select_best_epoch(history) is history[1]


def test_interrupted_fit_preserves_completed_epoch(torch_cpu, tmp_path):
    splits = toy_frames()
    train, val = (neural.encode_split(splits[name], "cpu") for name in ("train", "val"))
    output = tmp_path / "mlp"

    def interrupt(_):
        raise InterruptedError("stop after first persisted epoch")

    with pytest.raises(InterruptedError):
        neural.fit_model(neural.build_model("mlp"), train, val, epochs=2, seed=0,
                         output_dir=output, batch_size=3, log=interrupt)
    result = json.loads((output / "history.json").read_text())
    assert result["best_epoch"] == 1
    assert len(result["history"]) == 1
    checkpoint = output / result["best_checkpoint"]["file"]
    assert neural.file_source(checkpoint)["sha256"] == result["best_checkpoint"]["sha256"]
    assert torch_cpu.load(checkpoint, map_location="cpu", weights_only=True)


@pytest.mark.parametrize("conflict,match", [
    ("duplicate", "duplicate CDRH3"),
    ("sequence", "split overlap"),
    ("cluster", "split overlap.*cluster_centre"),
])
def test_loader_rejects_published_split_overlap(tmp_path, monkeypatch, conflict, match):
    splits = toy_frames()
    if conflict == "duplicate":
        splits["train"].loc[1, "seq"] = splits["train"].loc[0, "seq"]
    elif conflict == "sequence":
        splits["val"].loc[0, "seq"] = splits["train"].loc[0, "seq"]
    else:
        splits["test"].loc[0, "cluster_centre"] = splits["val"].loc[0, "cluster_centre"]
    monkeypatch.setattr(neural, "FACS_SHA256", write_toy_data(tmp_path, splits))
    with pytest.raises(ValueError, match=match):
        neural.load_facs_splits(tmp_path)


def test_loader_rejects_unpinned_split_bytes(tmp_path):
    write_toy_data(tmp_path, toy_frames())
    with pytest.raises(ValueError, match="Frozen input SHA-256 mismatch"):
        neural.load_facs_splits(tmp_path)


def test_primary_order_hash_and_overlap(tmp_path):
    ids, source = neural.load_primary_ids(neural.DEFAULT_OOF)
    assert len(ids) == len(set(ids)) == 158
    assert hashlib.sha256(("\n".join(ids) + "\n").encode()).hexdigest() == (
        "637ff639540aa806ad4940a1facf297934abce954520d47a275c752eccabfd00"
    )
    assert source["sha256"] == neural.OOF_SHA256
    changed = tmp_path / "oof.csv"
    changed.write_bytes(neural.DEFAULT_OOF.read_bytes() + b"\n")
    with pytest.raises(ValueError, match="Frozen input SHA-256 mismatch"):
        neural.load_primary_ids(changed)
    splits = toy_frames()
    splits["test"].loc[0, "seq"] = ids[0]
    audit = neural.audit_primary_overlap(ids, splits)
    assert audit["exact_cdrh3_match_count_by_split"] == {"train": 0, "val": 0, "test": 1}
    splits["val"].loc[0, "seq"] = ids[0]
    with pytest.raises(ValueError, match="Primary SPR IDs overlap FACS train/val"):
        neural.audit_primary_overlap(ids, splits)


def test_cpu_smoke_cli_writes_safe_checkpoints_and_loader_compatible_features(
    torch_cpu, tmp_path, monkeypatch
):
    from src.her2_interaction_benchmark import load_frozen_features

    source = tmp_path / "toy-organizer-data"
    output = tmp_path / "run"
    splits = toy_frames()
    monkeypatch.setattr(neural, "FACS_SHA256", write_toy_data(source, splits))
    args = ["--data-dir", str(source), "--output", str(output), "--epochs", "2", "--device", "cpu"]
    test_tokens = neural.encode_tokens(splits["test"]["seq"])
    test_calls = []
    original_predict = neural.predict

    def track_test(model, tokens, **kwargs):
        if np.array_equal(tokens.cpu().numpy(), test_tokens):
            name = "transformer" if isinstance(model, neural.FACSTransformer) else "mlp"
            assert (output / name / "best.pt").is_file()
            test_calls.append(name)
        return original_predict(model, tokens, **kwargs)

    monkeypatch.setattr(neural, "predict", track_test)
    neural.main(args)
    assert test_calls == ["transformer", "mlp"]
    metadata_path = output / "metadata.json"
    metadata = json.loads(metadata_path.read_text())
    ids, _ = neural.load_primary_ids(neural.DEFAULT_OOF)
    features, _ = load_frozen_features(output / "features.npz", metadata_path, ids)
    assert set(features) == {"facs_transformer", "facs_mlp"}
    assert all(matrix.shape == (158, 64) for matrix in features.values())
    probabilities = pd.read_csv(output / "primary_probabilities.csv")
    assert probabilities["sequence_id"].tolist() == ids
    assert metadata["configuration"]["seed"] == 0
    assert metadata["historical_facs_logistic"]["pr_auc_average_precision"] == 0.973718568352088
    assert metadata["audits"]["primary_spr_encoder_overlap"]["exact_cdrh3_match_count_by_split"] == {
        "train": 0, "val": 0, "test": 0
    }

    for name in neural.MODELS:
        metrics = metadata["models"][name]
        history = json.loads((output / name / "history.json").read_text())
        assert len(history["history"]) == 2
        best = neural.select_best_epoch(history["history"])
        assert metrics["best_epoch"] == best["epoch"]
        assert metrics["best_checkpoint"]["sha256"] == best["checkpoint"]["sha256"]
        for checkpoint in [metrics["best_checkpoint"], *(row["checkpoint"] for row in history["history"])]:
            path = output / name / checkpoint["file"]
            assert hashlib.sha256(path.read_bytes()).hexdigest() == checkpoint["sha256"]
            state = torch_cpu.load(path, map_location="cpu", weights_only=True)
            assert state and all(isinstance(value, torch_cpu.Tensor) for value in state.values())
        selected = torch_cpu.load(output / name / "best.pt", map_location="cpu", weights_only=True)
        model = neural.build_model(name)
        model.load_state_dict(selected)
        primary_tokens = torch_cpu.as_tensor(neural.encode_tokens(ids))
        replay_scores, replay_features = neural.predict(model, primary_tokens, include_features=True)
        np.testing.assert_array_equal(replay_features, features[f"facs_{name}"])
        np.testing.assert_allclose(replay_scores, probabilities[f"p_high_{name}"], rtol=1e-6)
    with pytest.raises(SystemExit, match="2"):
        neural.main(args)
    assert json.loads(metadata_path.read_text()) == metadata

    replay = tmp_path / "replay"
    neural.main(["--data-dir", str(source), "--output", str(replay), "--epochs", "2", "--models", "mlp"])
    replay_features, _ = load_frozen_features(replay / "features.npz", replay / "metadata.json", ids)
    assert set(replay_features) == {"facs_mlp"}
    np.testing.assert_array_equal(replay_features["facs_mlp"], features["facs_mlp"])
