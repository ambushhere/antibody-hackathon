"""Extract frozen paired AbLang2 features for the primary HER2 SPR cohort."""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import platform
import tempfile
from collections import defaultdict
from importlib.metadata import version
from pathlib import Path
from time import perf_counter

import gemmi
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OOF = ROOT / "results/spr-benchmark/oof_predictions.csv"
CANONICAL = ROOT / "configs/her2-context/P04626.fasta"
PROTOCOL = ROOT / "docs/research/interaction-protocol.md"
REVISION = "586af3083c32b5fb2f0a1c855f2ad4f1cad15ec3"
PACKAGE_SHA = "9fbb5d5cc62e7d9dfb8f99dfa8df81e8bb04779a15c82496e00c570a93cb11c2"
WEIGHTS_URL = "https://zenodo.org/records/10185169/files/ablang2-weights.tar.gz"
CIF_URL = "https://files.rcsb.org/download/1N8Z.cif"
FROZEN_SHA = {
    "oof": "44dfc97781fb80652aa59bfb8226e3b71fc501053d12c1d6e66cf8ea1c180fbf",
    "spr": "4d6f53a6183f377670e9577b3233a7c4c14888965ac63824863219d13424fa9b",
    "cif": "9f22f02b6e514de56dab8a48bc31a16302a51397779e64c398ad3a7f05cb19a2",
    "canonical": "544da5f05e6961a6a51f6fb7540d48b61f57cdfe5c8cfe4da455bb3118574cdb",
    "model.pt": "56d6f07862a6f824f88c8707bbc03e4026c9db762be2d3041e9767e2e6f86386",
    "hparams.json": "33474b76d78e01beb24b3621cd44d42b0813333d21e5b806168580bd836c5c8d",
}
PRIMARY = "Primary novelty cohort"
N_PRIMARY = 158
AMINO_ACIDS = frozenset("ACDEFGHIKLMNPQRSTVWY")
WT_LOOP = "WGGDGFYAMD"
# Frozen SPR scaffold, independently checked against deposited 1N8Z polymer sequences.
VH_PREFIX = (
    "EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKGRFTI"
    "SADTSKNTAYLQMNSLRAEDTAVYYCSR"
)
VH_SUFFIX = "YWGQGTLVTVSS"
WT_VH = VH_PREFIX + WT_LOOP + VH_SUFFIX
WT_VL = (
    "DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRFSGSRSG"
    "TDFTLTISSLQPEDFATYYCQQHYTTPPTFGQGTKVEIK"
)
VH_RESIDUE_INDICES = tuple(range(len(VH_PREFIX), len(VH_PREFIX) + len(WT_LOOP)))
# One explicit start token precedes VH. The hidden states retain all token positions.
H3_TOKEN_INDICES = tuple(index + 1 for index in VH_RESIDUE_INDICES)
N_TOKENS = len(WT_VH) + len(WT_VL) + 5
SPECIAL_IDS = {"<": 0, "-": 21, ">": 22, "*": 23, "X": 24, "|": 25}
WIDTH = 480
BATCH_SIZE = 4
CONTACT_CUTOFF = 4.5
HER2_OFFSET = 22


def require(condition, message):
    if not condition:
        raise ValueError(message)


def digest(data):
    return hashlib.sha256(data).hexdigest()


def check_sha(path, key):
    with path.open("rb") as handle:
        require(hashlib.file_digest(handle, "sha256").hexdigest() == FROZEN_SHA[key],
                f"Frozen {key} SHA-256 mismatch for {path}")


def relative_path(path):
    path = path.resolve()
    return str(path.relative_to(ROOT)) if path.is_relative_to(ROOT) else str(path)


def csv_rows(data):
    return csv.DictReader(io.StringIO(data.decode("utf-8-sig")))


def validate_pairs(oof_rows, spr_rows, expected_count=N_PRIMARY):
    ids = [row["cdrh3"] for row in oof_rows if row["cohort"] == PRIMARY]
    require(len(ids) == expected_count and len(set(ids)) == expected_count,
            f"Expected {expected_count} unique primary OOF IDs")
    pairs = defaultdict(set)
    for row in spr_rows:
        pairs[row["cdrh3"]].add((row["vh_sequence"], row["vl_sequence"]))
    require(pairs[WT_LOOP] == {(WT_VH, WT_VL)}, "SPR WT VH/VL pairing mismatch")
    for sequence_id in ids:
        require(len(sequence_id) == len(WT_LOOP) and set(sequence_id) <= AMINO_ACIDS
                and sequence_id != WT_LOOP, f"Invalid CDRH3 ID {sequence_id!r}")
        expected_vh = VH_PREFIX + sequence_id + VH_SUFFIX
        require(pairs[sequence_id] == {(expected_vh, WT_VL)},
                f"Missing or ambiguous SPR VH/VL context for {sequence_id}")
    return ids


def cif_rows(block, category):
    columns = block.get_mmcif_category(category)
    return [dict(zip(columns, values)) for values in zip(*columns.values())]


def coordinates(atoms):
    xyz = np.array([[float(atom[f"Cartn_{axis}"]) for axis in "xyz"] for atom in atoms])
    require(xyz.ndim == 2 and xyz.shape[1] == 3 and np.isfinite(xyz).all(),
            "Missing or malformed coordinate array")
    return xyz


def residue_letter(monomer):
    return gemmi.find_tabulated_residue(monomer).one_letter_code


def structural_mapping(cif, canonical_fasta):
    block = gemmi.cif.read_string(cif.decode()).sole_block()
    require(block.name == "1N8Z", "Expected the deposited 1N8Z structure")
    polymers = {row["entity_id"]: row for row in cif_rows(block, "_entity_poly.")}
    for entity, chain, sequence in (("2", "B", WT_VH), ("1", "A", WT_VL)):
        polymer = polymers[entity]
        deposited = "".join(polymer["pdbx_seq_one_letter_code_can"].split())
        require(polymer["pdbx_strand_id"] == chain and deposited.startswith(sequence),
                f"SPR variable domain does not match 1N8Z chain {chain}")

    header, *canonical_lines = canonical_fasta.decode().splitlines()
    require(header.startswith(">sp|P04626|") and " SV=1" in header,
            "Expected UniProt HER2 P04626 sequence version 1")
    canonical = "".join(canonical_lines)
    target_polymer = polymers["3"]
    require(target_polymer["pdbx_strand_id"] == "C"
            and "".join(target_polymer["pdbx_seq_one_letter_code_can"].split()) == canonical[22:629],
            "1N8Z chain C does not match canonical HER2 P04626 23–629")
    reference, = (row for row in cif_rows(block, "_struct_ref_seq.")
                  if row["pdbx_strand_id"] == "C")
    require(tuple(reference[key] for key in ("pdbx_db_accession", "pdbx_auth_seq_align_beg",
                                             "pdbx_auth_seq_align_end", "db_align_beg",
                                             "db_align_end")) == ("P04626", "1", "607", "23", "629"),
            "Unexpected 1N8Z to P04626 reference alignment")

    scheme = cif_rows(block, "_pdbx_poly_seq_scheme.")
    protein = [
        atom for atom in cif_rows(block, "_atom_site.")
        if atom["group_PDB"] == "ATOM" and atom["pdbx_PDB_model_num"] == "1"
        and atom["type_symbol"] not in ("H", "D") and float(atom["occupancy"]) > 0
    ]
    target = [atom for atom in protein if atom["label_asym_id"] == "C"]
    require(all(atom["auth_asym_id"] == "C" and atom["label_entity_id"] == "3"
                and atom["label_seq_id"] == atom["auth_seq_id"] and not atom["pdbx_PDB_ins_code"]
                and residue_letter(atom["label_comp_id"])
                == canonical[int(atom["auth_seq_id"]) + HER2_OFFSET - 1] for atom in target),
            "Ambiguous HER2 coordinate numbering")
    target_xyz = coordinates(target)

    mapping = []
    for index, (residue_index, token_index, aa) in enumerate(
        zip(VH_RESIDUE_INDICES, H3_TOKEN_INDICES, WT_LOOP), start=1
    ):
        position = str(residue_index + 1)
        entries = [row for row in scheme if row["asym_id"] == "B" and row["seq_id"] == position]
        require(len(entries) == 1, f"Ambiguous 1N8Z heavy polymer position {position}")
        entry = entries[0]
        require(entry["entity_id"] == "2" and entry["pdb_strand_id"] == "B"
                and entry["auth_seq_num"] == position and not entry["pdb_ins_code"]
                and residue_letter(entry["mon_id"]) == aa,
                f"WT loop position {index} does not map to 1N8Z heavy author position {position}")
        heavy = [atom for atom in protein if atom["label_asym_id"] == "B"
                 and atom["label_seq_id"] == position]
        require(any(atom["label_atom_id"] == "CA" for atom in heavy)
                and all(atom["auth_asym_id"] == "B" and atom["label_entity_id"] == "2"
                        and atom["auth_seq_id"] == position and not atom["pdbx_PDB_ins_code"]
                        and residue_letter(atom["auth_comp_id"]) == aa for atom in heavy),
                f"Missing or ambiguous 1N8Z heavy coordinate residue B{position}")
        distance_sq = ((coordinates(heavy)[:, None, :] - target_xyz[None, :, :]) ** 2).sum(axis=2)
        target_contacts = sorted({
            int(target[j]["auth_seq_id"])
            for j in np.flatnonzero(distance_sq.min(axis=0) < CONTACT_CUTOFF ** 2)
        })
        mapping.append({
            "window_position_1based": index,
            "wt_aa": aa,
            "vh_residue_index_0based": residue_index,
            "model_output_token_index_0based": token_index,
            "1n8z_heavy_author_chain": "B",
            "1n8z_heavy_author_residue": int(position),
            "1n8z_heavy_label_seq_id": int(position),
            "minimum_her2_heavy_atom_distance_angstrom": float(np.sqrt(distance_sq.min())),
            "1n8z_her2_author_contacts": target_contacts,
            "p04626_her2_contacts": [number + HER2_OFFSET for number in target_contacts],
        })
    contact_indices = [
        row["window_position_1based"] - 1 for row in mapping if row["p04626_her2_contacts"]
    ]
    require(contact_indices, "No physical HER2 contacts in the frozen window10")
    return mapping, contact_indices


def format_pairs(ids):
    return [f"<{VH_PREFIX}{sequence_id}{VH_SUFFIX}>|<{WT_VL}>" for sequence_id in ids]


def tokenize(tokenizer, ids):
    require({key: tokenizer.aa_to_token[key] for key in SPECIAL_IDS} == SPECIAL_IDS,
            "Unexpected special token IDs")
    pairs = format_pairs(ids)
    tokens = tokenizer(pairs, pad=True, w_extra_tkns=False, device="cpu")
    require(tuple(tokens.shape) == (len(ids), N_TOKENS), "Unexpected paired token array shape")
    for sequence_id, row, pair in zip(ids, tokens, pairs):
        require(tokenizer.decode(row) == pair, f"Pair token round-trip failed for {sequence_id}")
        require(tokenizer.decode(row[list(H3_TOKEN_INDICES)]) == sequence_id,
                f"Window10 token indexing failed for {sequence_id}")
        require(not set(row.tolist()).intersection({21, 23, 24}),
                f"Padding, mask or unknown token found for {sequence_id}")
    return tokens


def pool_embeddings(hidden, batch_size, contact_indices):
    require(hidden.shape == (batch_size, N_TOKENS, WIDTH), "Unexpected AbRep hidden state shape")
    require(hidden.dtype == np.float32 and np.isfinite(hidden).all(),
            "AbRep hidden states must be finite float32")
    require(len(contact_indices) > 0 and sorted(set(contact_indices)) == list(contact_indices)
            and all(0 <= index < len(WT_LOOP) for index in contact_indices),
            "Invalid or empty WT contact mask")
    h3 = hidden[:, H3_TOKEN_INDICES, :]
    return h3.mean(axis=1), h3[:, contact_indices, :].mean(axis=1)


def load_model(weights):
    import ablang2
    import torch
    from ablang2.models.ablang2.ablang import AbLang
    from ablang2.models.ablang2.tokenizers import ABtokenizer

    package = Path(ablang2.__file__).resolve().parent
    code_hash = hashlib.sha256()
    for path in sorted(package.rglob("*.py")):
        code_hash.update(path.relative_to(package).as_posix().encode() + b"\0")
        code_hash.update(path.read_bytes() + b"\0")
    require(code_hash.hexdigest() == PACKAGE_SHA, "Installed AbLang2 source differs from pinned revision")
    for name in ("model.pt", "hparams.json"):
        check_sha(weights / name, name)
    params = json.loads((weights / "hparams.json").read_text())
    # Match the constructor arguments in upstream fetch_ablang2 with an explicit safe local load.
    model = AbLang(
        vocab_size=params["vocab_size"], hidden_embed_size=params["hidden_embed_size"],
        n_attn_heads=params["n_attn_heads"], n_encoder_blocks=params["n_encoder_blocks"],
        padding_tkn=params["pad_tkn"], mask_tkn=params["mask_tkn"],
        layer_norm_eps=params["layer_norm_eps"], a_fn=params["a_fn"],
    )
    state = torch.load(weights / "model.pt", map_location="cpu", weights_only=True)
    model.load_state_dict(state, strict=True)
    model.eval().requires_grad_(False)
    require(all(not module.training for module in model.modules()), "Model is not entirely in eval mode")
    require(all(parameter.device.type == "cpu" and parameter.dtype == torch.float32
                for parameter in model.parameters()), "Expected CPU float32 model parameters")
    return model, ABtokenizer(), params


def write_output(output, ids, primary, ablation, metadata):
    output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".ablang-", dir=output.parent) as temporary:
        staging = Path(temporary) / "complete"
        staging.mkdir()
        npz = staging / "features.npz"
        np.savez_compressed(npz, sequence_ids=np.asarray(ids, dtype="<U10"),
                            cdrh3_mean=primary, wt_contact_mean=ablation)
        metadata["output"] = {"file": npz.name, "sha256": digest(npz.read_bytes()),
                              "feature_shape_each": list(primary.shape)}
        (staging / "metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")
        staging.rename(output)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--spr", required=True, type=Path, help="Original frozen SPR CSV")
    parser.add_argument("--cif", required=True, type=Path, help="Deposited 1N8Z mmCIF")
    parser.add_argument("--weights-dir", type=Path, help="Directory with verified model.pt and hparams.json")
    parser.add_argument("--output-dir", type=Path, default=ROOT / "results/ablang-features")
    modes = parser.add_mutually_exclusive_group()
    modes.add_argument("--validate-only", action="store_true", help="Verify cohort and geometry without Torch")
    modes.add_argument("--smoke", action="store_true", help="Encode only the first batch without writing features")
    args = parser.parse_args()
    if not args.validate_only and args.weights_dir is None:
        parser.error("--weights-dir is required for inference")
    if not args.validate_only and not args.smoke and args.output_dir.exists():
        parser.error("--output-dir already exists; use a fresh directory")

    paths = {"oof": OOF, "spr": args.spr, "cif": args.cif, "canonical": CANONICAL}
    payloads = {key: path.read_bytes() for key, path in paths.items()}
    for key, data in payloads.items():
        require(digest(data) == FROZEN_SHA[key], f"Frozen {key} SHA-256 mismatch for {paths[key]}")
    ids = validate_pairs(csv_rows(payloads["oof"]), csv_rows(payloads["spr"]))
    mapping, contact_indices = structural_mapping(payloads["cif"], payloads["canonical"])
    metadata = {
        "schema_version": 1,
        "inputs": {key: {"path": relative_path(path), "sha256": FROZEN_SHA[key]}
                   for key, path in paths.items()},
        "extractor": {"path": relative_path(Path(__file__)), "sha256": digest(Path(__file__).read_bytes())},
        "interaction_protocol_sha256": digest(PROTOCOL.read_bytes()),
        "cohort": PRIMARY,
        "ordered_sequence_ids": ids,
        "pairing": {"heavy_prefix": VH_PREFIX, "heavy_suffix": VH_SUFFIX, "wt_window10": WT_LOOP,
                    "light_sequence": WT_VL, "token_format": "<VH>|<VL>", "token_count": N_TOKENS,
                    "special_token_ids": SPECIAL_IDS, "heavy_length": len(WT_VH), "light_length": len(WT_VL),
                    "model_input_sha256": digest(("\n".join(format_pairs(ids)) + "\n").encode("ascii"))},
        "structural_mask": {"pdb": "1N8Z", "url": CIF_URL, "model": 1,
                            "definition": "Observed positive occupancy polymer non-H/D atom pairs to HER2 chain C",
                            "distance_strictly_less_than_angstrom": CONTACT_CUTOFF,
                            "her2_p04626_minus_author": HER2_OFFSET,
                            "window_indices_0based": contact_indices, "window_residues": mapping},
        "pooling": {"cdrh3_mean": "Mean of final AbRep hidden states at the ten window positions",
                    "wt_contact_mean": "Mean of the same states restricted to the WT contact mask"},
    }
    if args.validate_only:
        print(json.dumps(metadata, indent=2))
        return

    import torch

    torch.set_num_threads(1)
    torch.set_num_interop_threads(1)
    torch.use_deterministic_algorithms(True)
    torch.manual_seed(0)
    torch.set_default_dtype(torch.float32)
    model, tokenizer, params = load_model(args.weights_dir)
    metadata["model"] = {
        "name": "ablang2-paired", "implementation": f"https://github.com/oxpig/AbLang2/tree/{REVISION}",
        "revision": REVISION, "installed_python_tree_sha256": PACKAGE_SHA, "code_license": "BSD-3-Clause",
        "weights_record": "10.5281/zenodo.10185169", "weights_url": WEIGHTS_URL,
        "weights_license": "BSD-3-Clause-Clear",
        "archive_sha256": "7acb23b880e3b69b3fea8b7d6cacac0fda147daafd99edc411e0ff5c61252ac5",
        "archive_zenodo_md5": "6425b35e9b83750fde67a6a6d240d996",
        "weights_sha256": {name: FROZEN_SHA[name] for name in ("model.pt", "hparams.json")},
        "hparams": params,
    }
    metadata["runtime"] = {
        "python": platform.python_version(), "platform": platform.platform(), "device": "cpu",
        "dtype": "float32", "batch_size": BATCH_SIZE, "threads": 1, "interop_threads": 1,
        "deterministic_algorithms": True, "eval": True, "masking": False, "seed": 0,
        "versions": {name: version(name) for name in
                     ("ablang2", "torch", "numpy", "gemmi", "einops", "rotary-embedding-torch", "requests")},
    }
    selected_ids = ids[:BATCH_SIZE] if args.smoke else ids
    primary, ablation = [], []
    started = perf_counter()
    with torch.inference_mode():
        for start in range(0, len(selected_ids), BATCH_SIZE):
            batch_ids = selected_ids[start:start + BATCH_SIZE]
            tokens = tokenize(tokenizer, batch_ids)
            hidden = model.AbRep(tokens).last_hidden_states.cpu().numpy()
            full_mean, contact_mean = pool_embeddings(hidden, len(batch_ids), contact_indices)
            primary.append(full_mean)
            ablation.append(contact_mean)
    elapsed = perf_counter() - started
    print(f"Encoded {len(selected_ids)} paired VH/VL contexts in {elapsed:.2f} s on one CPU thread")
    if not args.smoke:
        metadata["runtime"]["inference_seconds"] = round(elapsed, 3)
        write_output(args.output_dir, ids, np.concatenate(primary), np.concatenate(ablation), metadata)


if __name__ == "__main__":
    main()
