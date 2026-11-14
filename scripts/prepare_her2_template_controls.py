"""Build the matched, target-only 6OGE template ablation on CPU."""

import argparse
import json
from pathlib import Path

import gemmi

from scripts.prepare_her2_context import (
    CONFIG as CONTEXT, ROOT, canonical_sequence, chain_sequences, json_text,
    require, residue_aa, rows, sha256, write_rows,
)

CONFIG = ROOT / "configs/her2-template-controls"
SOURCE_URL = "https://files.rcsb.org/download/6OGE.cif"
SOURCE_SHA256 = "33a32b66e64aac934a7cfebfd78903027d3b4bf2a7f510e37c2ace4e4b52f18a"
POSITIONS = range(511, 630)
MASK = range(603, 613)
OFFSET = 22  # 6OGE chain A deposited label index + 22 = author position = P04626.
PANEL = ("control_wt", "control_facs_low", "control_spr_strong", "control_spr_weak")
TEMPLATES = {"complete": "complete-target.cif", "masked": "masked-target.cif"}
QUERY_FILE = "queries.json"
SEEDS = (42, 43, 44)
SAMPLES = (1, 2)


def frozen_context():
    manifest_path = CONTEXT / "manifest.json"
    manifest = json.loads(manifest_path.read_bytes())
    fasta = (CONTEXT / "P04626.fasta").read_bytes()
    require(sha256(fasta) == manifest["sources"]["p04626_fasta"]["sha256"],
            "Frozen P04626 FASTA hash mismatch")
    canonical = canonical_sequence(fasta)
    require(sha256(canonical.encode()) == manifest["canonical_sequence_sha256"],
            "Frozen P04626 sequence hash mismatch")
    query_path = CONTEXT / "domain-iv-queries.json"
    data = query_path.read_bytes()
    require(sha256(data) == manifest["artifact_sha256"][query_path.name],
            "Frozen domain IV queries hash mismatch")
    original = json.loads(data)["queries"]
    require(tuple(q["query_id"] for q in original if q["query_id"] in PANEL) == PANEL,
            "Frozen control coverage or order changed")
    panel = [row for row in manifest["panel"] if row["query_id"] in PANEL]
    require(tuple(row["query_id"] for row in panel) == PANEL, "Frozen manifest panel changed")
    queries = []
    for query, record in zip(original[:len(PANEL)], panel, strict=True):
        require(query["query_id"] == record["query_id"], "Frozen query order changed")
        chains = chain_sequences(query)
        require(list(chains) == ["T", "H", "L"], "Frozen query chain order changed")
        for role, name in (("H", "vh"), ("L", "vl")):
            require(chains[role] == record[f"{name}_sequence"]
                    and sha256(chains[role].encode()) == record[f"{name}_sha256"],
                    f"{query['query_id']}: frozen sequence or hash mismatch")
        chains["T"] = canonical[POSITIONS.start - 1:POSITIONS.stop - 1]
        queries.append({"query_id": query["query_id"], "chains": [
            {"chain_id": role, "sequence": sequence, "msa_path": None}
            for role, sequence in chains.items()
        ]})
    sources = {
        str(path.relative_to(ROOT)): sha256(path.read_bytes())
        for path in (manifest_path, query_path, CONTEXT / "P04626.fasta")
    }
    return canonical, panel, {"queries": queries}, sources


def source_target(block, canonical):
    reference, = [r for r in rows(block, "_struct_ref.")
                  if r["entity_id"] == "1" and r["pdbx_db_accession"] == "P04626"]
    alignment, = [r for r in rows(block, "_struct_ref_seq.") if r["ref_id"] == reference["id"]]
    require(block.name.upper() == "6OGE" and alignment["pdbx_strand_id"] == "A"
            and alignment["pdbx_db_accession"] == "P04626"
            and [int(alignment[k]) for k in ("seq_align_beg", "seq_align_end",
                                             "db_align_beg", "db_align_end",
                                             "pdbx_auth_seq_align_beg", "pdbx_auth_seq_align_end")]
            == [1, 622, 23, 644, 23, 644], "6OGE chain A deposited alignment changed")
    polymer, = [r for r in rows(block, "_entity_poly.") if r["entity_id"] == "1"]
    require(polymer["pdbx_strand_id"] == "A"
            and "".join(polymer["pdbx_seq_one_letter_code_can"].split()) == canonical[22:644],
            "6OGE chain A polymer sequence differs from P04626")
    atoms = [a for a in rows(block, "_atom_site.")
             if a["group_PDB"] == "ATOM" and a["label_asym_id"] == "A"]
    require(all(a["auth_asym_id"] == "A" and a["label_entity_id"] == "1"
                and a["pdbx_PDB_model_num"] == "1" and not a["pdbx_PDB_ins_code"]
                and not a["label_alt_id"] and int(a["auth_seq_id"]) == int(a["label_seq_id"]) + OFFSET
                and a["auth_comp_id"] == a["label_comp_id"]
                and residue_aa(a["label_comp_id"]) == canonical[int(a["auth_seq_id"]) - 1]
                for a in atoms), "6OGE chain A coordinate mapping changed")
    return atoms


def make_templates(block, canonical):
    sequence = canonical[POSITIONS.start - 1:POSITIONS.stop - 1]
    deposited = [a for a in source_target(block, canonical) if int(a["auth_seq_id"]) in POSITIONS]
    require([int(a["auth_seq_id"]) for a in deposited if a["label_atom_id"] == "CA"] == list(POSITIONS),
            "6OGE crop lacks exact C-alpha coverage")
    keys = {(a["auth_seq_id"], a["label_atom_id"]) for a in deposited}
    require(len(keys) == len(deposited), "Duplicate deposited atoms in 6OGE crop")

    scheme = [r for r in rows(block, "_pdbx_poly_seq_scheme.")
              if r["asym_id"] == "A" and int(r["seq_id"]) + OFFSET in POSITIONS]
    require(len(scheme) == len(POSITIONS), "Incorrect 6OGE scheme coverage")
    seq = []
    for i, (pos, record) in enumerate(zip(POSITIONS, scheme, strict=True), 1):
        require(int(record["seq_id"]) + OFFSET == int(record["pdb_seq_num"])
                == int(record["auth_seq_num"]) == pos
                and record["pdb_strand_id"] == "A" and record["entity_id"] == "1"
                and residue_aa(record["mon_id"]) == sequence[i - 1], "6OGE scheme mapping changed")
        record.update(seq_id=str(i), ndb_seq_num=str(i))
        seq.append(dict(entity_id="1", num=str(i), mon_id=record["mon_id"], hetero="n"))

    internal, crossing = [], []
    for record in rows(block, "_struct_conn."):
        if (record["conn_type_id"] != "disulf"
                or record["ptnr1_label_asym_id"] != "A" or record["ptnr2_label_asym_id"] != "A"):
            continue
        p, q = (int(record[f"ptnr{i}_auth_seq_id"]) for i in (1, 2))
        if p in POSITIONS and q in POSITIONS:
            require(p not in MASK and q not in MASK, "Coordinate mask intersects an internal disulfide")
            for i, pos in enumerate((p, q), 1):
                require(int(record[f"ptnr{i}_label_seq_id"]) + OFFSET == pos
                        and (str(pos), "SG") in keys, "Incorrect 6OGE disulfide endpoint")
                record[f"ptnr{i}_label_seq_id"] = str(pos - POSITIONS.start + 1)
            internal.append(record)
        elif p in POSITIONS or q in POSITIONS:
            crossing.append([p, q])

    components = {r["mon_id"] for r in seq}
    shared = {
        "_entity.": [dict(id="1", type="polymer", pdbx_description="HER2 target")],
        "_entity_poly.": [dict(entity_id="1", type="polypeptide(L)", nstd_linkage="no",
                               nstd_monomer="no", pdbx_seq_one_letter_code=sequence,
                               pdbx_seq_one_letter_code_can=sequence, pdbx_strand_id="A")],
        "_struct_asym.": [dict(id="A", entity_id="1")],
        "_entity_poly_seq.": seq,
        "_pdbx_poly_seq_scheme.": scheme,
        "_chem_comp.": [r for r in rows(block, "_chem_comp.") if r["id"] in components],
        "_struct_conn.": internal,
    }
    templates, arms = {}, {}
    for arm, filename in TEMPLATES.items():
        doc = gemmi.cif.Document()
        crop = doc.add_new_block("6OGE")
        for category in ("_entry.", "_exptl.", "_pdbx_audit_revision_history.", "_refine."):
            data = block.get_mmcif_category(category)
            if data:
                crop.set_mmcif_category(category, data)
        for category, records in shared.items():
            write_rows(crop, category, records)
        mask = tuple(MASK) if arm == "masked" else ()
        coords = [dict(a, label_seq_id=str(int(a["auth_seq_id"]) - POSITIONS.start + 1))
                  for a in deposited if int(a["auth_seq_id"]) not in mask]
        write_rows(crop, "_atom_site.", coords)
        templates[filename] = doc.as_string()
        arms[arm] = {"template": filename, "query": QUERY_FILE,
                     "masked_canonical_positions": list(mask),
                     "coordinate_residues": len(POSITIONS) - len(mask), "atom_count": len(coords)}
    disulfides = [[int(r[f"ptnr{i}_auth_seq_id"]) for i in (1, 2)] for r in internal]
    return templates, arms, disulfides, crossing


def prepare(source_path, output):
    data = source_path.read_bytes()
    # ASVS 2.2.1: allow only the audited RCSB deposit before parsing mmCIF.
    require(sha256(data) == SOURCE_SHA256, "Original RCSB 6OGE mmCIF SHA-256 mismatch")
    block = gemmi.cif.read_string(data.decode()).sole_block()
    canonical, panel, queries, sources = frozen_context()
    templates, arms, disulfides, crossing = make_templates(block, canonical)
    artifacts = {QUERY_FILE: json_text(queries), **templates}
    manifest = {
        "sources": {SOURCE_URL: SOURCE_SHA256, **sources},
        "reference": {
            "url": "https://files.rcsb.org/download/1N8Z.cif",
            "sha256": json.loads((CONTEXT / "manifest.json").read_bytes())["sources"]["pdb_cif"]["sha256"],
            "common_canonical_positions": [p for p in POSITIONS if p not in MASK],
        },
        "data_origin": (
            "Both 1N8Z and 6OGE appear in upstream RFD4 sampling metadata "
            "(results/rfd4-training-metadata/audit.json). This is a reproduction control. "
            "Exact checkpoint exposure is unestablished."
        ),
        "target": {
            "accession": "P04626", "query_chain": "T", "deposited_label_chain": "A",
            "deposited_author_chain": "A", "entity_id": "1",
            "positions": [
                {"query_index": i, "template_label_seq_id": i, "deposited_label_seq_id": p - OFFSET,
                 "deposited_auth_seq_id": p, "canonical_position": p, "aa": canonical[p - 1]}
                for i, p in enumerate(POSITIONS, 1)
            ],
            "internal_disulfides_p04626": disulfides,
            "disulfides_crossing_crop_boundary_p04626": crossing,
        },
        "panel": panel,
        "arms": arms,
        "run": {"model": "protenix-v2", "msa": False, "seeds": list(SEEDS),
                "samples": list(SAMPLES), "expected_poses": len(PANEL) * len(TEMPLATES) * len(SEEDS) * len(SAMPLES)},
        "artifact_sha256": {name: sha256(payload.encode()) for name, payload in artifacts.items()},
    }
    output.mkdir(parents=True, exist_ok=True)
    for name, payload in {**artifacts, "manifest.json": json_text(manifest)}.items():
        (output / name).write_bytes(payload.encode())
    return manifest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-cif", type=Path, required=True, help=SOURCE_URL)
    parser.add_argument("--output", type=Path, default=CONFIG)
    args = parser.parse_args()
    prepare(args.source_cif, args.output)


if __name__ == "__main__":
    main()
