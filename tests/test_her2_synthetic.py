"""Check the HER2 synthetic proposal and frozen pseudolabel contracts."""

import csv
import hashlib
import json
from pathlib import Path
from zipfile import ZipFile

import numpy as np
import pytest

from scripts import generate_her2_synthetic as synthetic


@pytest.fixture
def sources(tmp_path, monkeypatch):
    splits = {
        "train": {"AGGDGFYAMD", "AAGDGFYAMD", "CGGDGFYAMD"},
        "val": {"SAGDGFYAMD"},
        "test": {"GAGDGFYAMD"},
    }
    archive_path = tmp_path / synthetic.SPLIT_ARCHIVE
    archive_path.parent.mkdir()
    with ZipFile(archive_path, "w") as archive:
        for name, ids in splits.items():
            # Deliberately invalid values verify the generator never needs FACS labels.
            archive.writestr(
                f"{synthetic.facs.SPLIT_PATH}{name}.csv",
                "seq,label,class\n" + "".join(f"{seq},unused,unused\n" for seq in sorted(ids)),
            )
    oof = tmp_path / "oof.csv"
    oof.write_text(
        f"cdrh3,cohort,neg_log10_kd_nm\nAAGAGFYAMD,{synthetic.PRIMARY_COHORT},unused\n"
        "NAGDGFYAMD,Other cohort,unused\n",
        encoding="utf-8",
    )
    monkeypatch.setattr(synthetic, "OOF_SHA256", synthetic.sha256(oof))
    monkeypatch.setattr(synthetic, "PRIMARY_COUNT", 1)
    return tmp_path, oof


def test_generation_is_reproducible_audited_and_balanced(sources):
    splits, primary, inputs = synthetic.load_sequence_sources(*sources)
    assert primary == {"AAGAGFYAMD"}
    proposals, stats = synthetic.generate_proposals(splits, primary, count=301, seed=0)
    reordered = {name: tuple(reversed(sorted(ids))) for name, ids in splits.items()}
    assert (proposals, stats) == synthetic.generate_proposals(reordered, primary, count=301, seed=0)
    assert proposals != synthetic.generate_proposals(splits, primary, count=301, seed=1)[0]
    assert sorted(stats["accepted_by_mutation_count"].values()) == [100, 100, 101]
    assert stats["attempts"] == len(proposals) + sum(stats["rejections_first_match"].values())
    assert stats["rejections_first_match"]["cysteine"] > 0
    assert any("C" in p.parent_sequence for p in proposals)
    assert len({p.sequence_id for p in proposals}) == 301
    for proposal in proposals:
        assert proposal.parent_sequence in splits["train"]
        assert len(proposal.sequence_id) == 10
        assert "C" not in proposal.sequence_id
        assert proposal.sequence_id != synthetic.facs.WT_CDRH3
        assert proposal.sequence_id not in primary
        assert all(proposal.sequence_id not in ids for ids in splits.values())
        assert proposal.mutation_count == sum(
            a != b for a, b in zip(proposal.sequence_id, proposal.parent_sequence, strict=True)
        )
    assert synthetic.audit_proposals(proposals, splits, primary, 301)["exact_overlaps"] == {
        "wild_type": 0, "facs_train": 0, "facs_val": 0, "facs_test": 0, "spr_primary": 0,
    }
    with ZipFile(sources[0] / synthetic.SPLIT_ARCHIVE) as archive:
        assert inputs["val"]["sha256"] == hashlib.sha256(
            archive.read(f"{synthetic.facs.SPLIT_PATH}val.csv")
        ).hexdigest()


def test_oof_fingerprint_mismatch_fails(sources):
    _, oof = sources
    oof.write_bytes(oof.read_bytes().replace(b"unused", b"modified"))
    with pytest.raises(ValueError, match="SPR OOF fingerprint mismatch"):
        synthetic.load_sequence_sources(*sources)


def test_generation_cannot_loop_indefinitely():
    splits = {"train": {"CCCCCCCCCC"}, "val": set(), "test": set()}
    with pytest.raises(ValueError, match=r"attempt limit reached \(7\); accepted 0"):
        synthetic.generate_proposals(splits, set(), count=5, seed=0, max_attempts=7)


@pytest.mark.parametrize("sequence_id", [
    "AAGDGFYAMD", "SAGDGFYAMD", "GAGDGFYAMD", "AAGAGFYAMD", synthetic.facs.WT_CDRH3,
])
def test_audit_rejects_each_excluded_cohort(sources, sequence_id):
    splits, primary, _ = synthetic.load_sequence_sources(*sources)
    parent = "AGGDGFYAMD"
    mutations = sum(a != b for a, b in zip(sequence_id, parent, strict=True))
    with pytest.raises(ValueError, match="exclusion audit failed"):
        synthetic.audit_proposals([synthetic.Proposal(sequence_id, parent, mutations)], splits, primary, 1)


@pytest.mark.parametrize("probabilities", [
    [[0.4, 0.4, 0.4]], [[float("nan"), 0.2, 0.8]], [[-0.1, 0.2, 0.9]], [[0.2, 0.8]],
])
def test_probability_audit_rejects_invalid_rows(probabilities):
    with pytest.raises(ValueError, match="[Tt]eacher probabilities"):
        synthetic.audit_probabilities(np.asarray(probabilities, dtype=np.float32), 1)


def test_dataset_round_trip_and_refuses_overwrite(sources):
    splits, primary, _ = synthetic.load_sequence_sources(*sources)
    proposals, _ = synthetic.generate_proposals(splits, primary, count=6, seed=0)
    probabilities = np.array([[0.2, 0.0, 0.8]] * 6, dtype=np.float32)
    audit = synthetic.audit_probabilities(probabilities, 6)
    assert audit["argmax_count_by_class"] == {"low": 0, "mid": 0, "high": 6}
    assert audit["entropy_nats"]["median"] == pytest.approx(0.5004024)
    output = sources[0] / "synthetic"
    synthetic.write_dataset(output, proposals, probabilities, {"label_type": "pseudolabels"})
    with np.load(output / "synthetic.npz", allow_pickle=False) as dataset:
        assert dataset["sequence_ids"].dtype.kind == "U"
        assert dataset["sequence_ids"].tolist() == [p.sequence_id for p in proposals]
        assert dataset["probabilities"].dtype == np.float32
        np.testing.assert_array_equal(dataset["probabilities"], probabilities)
    csv_path = output / "proposals.csv"
    assert b"\r" not in csv_path.read_bytes()
    with csv_path.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    assert [r["parent_sequence"] for r in rows] == [p.parent_sequence for p in proposals]
    assert [int(r["mutation_count"]) for r in rows] == [p.mutation_count for p in proposals]
    before = (output / "metadata.json").read_bytes()
    manifest = json.loads(before)["output"]
    assert manifest["sha256"] == synthetic.sha256(output / "synthetic.npz")
    assert manifest["csv"]["sha256"] == synthetic.sha256(csv_path)
    with pytest.raises(FileExistsError):
        synthetic.write_dataset(output, proposals, probabilities, {"label_type": "changed"})
    assert (output / "metadata.json").read_bytes() == before


@pytest.fixture
def teacher_fixture(sources, monkeypatch):
    torch = pytest.importorskip("torch")
    from scripts import train_her2_neural as neural

    data_dir, oof = sources
    _, _, inputs = synthetic.load_sequence_sources(data_dir, oof)
    monkeypatch.setattr(neural, "FACS_SHA256", {
        **{name: inputs[name]["sha256"] for name in synthetic.facs.SPLIT_NAMES},
        "facs": neural.FACS_SHA256["facs"],
    })
    teacher_run = data_dir / "teacher"
    model_dir = teacher_run / "mlp"
    model_dir.mkdir(parents=True)
    model = neural.build_model("mlp", "large")
    with torch.no_grad():
        for parameter in model.parameters():
            parameter.zero_()
        model.head.bias.copy_(torch.tensor([0.0, 1.0, 2.0]))
    checkpoint = neural.save_checkpoint(model_dir / "best.pt", model)
    metadata = {
        "schema_version": 1,
        "status": "complete",
        "code": {
            "trainer": {"sha256": synthetic.sha256(Path(neural.__file__))},
            "facs_benchmark": {"sha256": synthetic.sha256(Path(synthetic.facs.__file__))},
        },
        "configuration": {
            "size": "large",
            "models": ["mlp"],
            "architecture": neural.ARCHITECTURES["large"],
            "class_to_index": neural.CLASS_TO_INDEX,
            "residue_alphabet": synthetic.facs.AMINO_ACIDS,
            "checkpoint_selection": synthetic.TEACHER_SELECTION,
            "seed": 0,
        },
        "inputs": {
            **inputs,
            "facs": {"sha256": neural.FACS_SHA256["facs"]},
            "frozen_oof": inputs["oof"],
        },
        "models": {"mlp": {
            "latent_width": model.latent_width,
            "best_epoch": 1,
            "best_checkpoint": checkpoint,
            "history": [
                {"epoch": 1, "validation_high_ap": 0.75, "checkpoint": {**checkpoint, "file": "epoch-01.pt"}},
                {"epoch": 2, "validation_high_ap": 0.50, "checkpoint": {"sha256": "unused"}},
            ],
        }},
    }
    (teacher_run / "metadata.json").write_text(json.dumps(metadata), encoding="utf-8")
    return teacher_run, metadata, inputs


def test_cli_cpu_smoke_with_frozen_teacher(sources, teacher_fixture):
    data_dir, oof = sources
    teacher_run, _, _ = teacher_fixture
    output = data_dir / "smoke"
    synthetic.main([
        "--data-dir", str(data_dir), "--oof", str(oof), "--teacher-run", str(teacher_run),
        "--output", str(output), "--count", "9", "--device", "cpu",
    ])
    expected = np.exp(np.arange(3, dtype=np.float32))
    expected /= expected.sum()
    with np.load(output / "synthetic.npz", allow_pickle=False) as result:
        assert result["probabilities"].dtype == np.float32
        np.testing.assert_allclose(result["probabilities"], np.tile(expected, (9, 1)), rtol=1e-6)
    metadata = json.loads((output / "metadata.json").read_text())
    assert metadata["seed"] == 0
    assert metadata["label_type"] == "pseudolabels"
    assert metadata["class_order"] == ["low", "mid", "high"]
    assert metadata["teacher"]["selected_epoch"] == 1
    assert metadata["teacher"]["checkpoint"]["sha256"] == synthetic.sha256(teacher_run / "mlp/best.pt")


@pytest.mark.parametrize("target, expected", [
    ("train", "FACS train split fingerprint mismatch"),
    ("val", "FACS val split fingerprint mismatch"),
    ("test", "FACS test split fingerprint mismatch"),
    ("source", "trainer source fingerprint mismatch"),
    ("checkpoint", "checkpoint fingerprint mismatch"),
    ("class_order", "class order"),
])
def test_teacher_rejects_inconsistent_provenance(teacher_fixture, target, expected):
    teacher_run, metadata, inputs = teacher_fixture
    if target in synthetic.facs.SPLIT_NAMES:
        metadata["inputs"][target]["sha256"] = "changed"
    elif target == "source":
        metadata["code"]["trainer"]["sha256"] = "changed"
    elif target == "checkpoint":
        with (teacher_run / "mlp/best.pt").open("ab") as handle:
            handle.write(b"changed")
    else:
        metadata["configuration"]["class_to_index"] = {"high": 0, "mid": 1, "low": 2}
    (teacher_run / "metadata.json").write_text(json.dumps(metadata), encoding="utf-8")
    with pytest.raises(ValueError, match=expected):
        synthetic.load_teacher(teacher_run, inputs, device="cpu")
