"""Prepare matched HER2 queries and target-only 1N8Z templates without running inference."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import textwrap
from pathlib import Path
from urllib.request import urlopen

import gemmi
import numpy as np
from scipy.spatial.distance import cdist

ROOT = Path(__file__).resolve().parents[1]
CONFIG = ROOT / "configs/her2-context"
CAMPAIGN = ROOT / "results/campaign/campaign-20261114"
RESULTS = ROOT / "results/experiment-20261114"
P04626 = "https://rest.uniprot.org/uniprotkb/P04626.fasta"
PDB_CIF = "https://files.rcsb.org/download/1N8Z.cif"
RCSB = "https://data.rcsb.org/rest/v1/core/polymer_entity/1N8Z/3"
PDBE = "https://www.ebi.ac.uk/pdbe/api/mappings/uniprot/1n8z"
OFFSET = 22
DEPOSITED_END = 607
CROPS = {"C": (489, 580), "D": (591, 607)}
ARMS = {"split": CROPS, "domain-iv": {"T": (489, 630)}, "ecd": {"T": (1, 630)}}
WT_CDRH3 = "WGGDGFYAMD"
DE_NOVO = "denovo_fullseq_k1_1"
DESIGN_SPEC = RESULTS / "input-specs/denovo-full-sequence.json"
HOTSPOT_SPEC_KEY = "her2_denovo_full_sequence"


def require(condition, message):
    if not condition:
        raise ValueError(message)


def sha256(data):
    return hashlib.sha256(data).hexdigest()


def read_source(path, url):
    if path is not None:
        return path.read_bytes()
    with urlopen(url, timeout=30) as response:
        return response.read()


def rows(block, category):
    columns = block.get_mmcif_category(category)
    return [dict(zip(columns, values)) for values in zip(*columns.values())]


def write_rows(block, category, records):
    require(bool(records), f"Empty required category {category}")
    block.set_mmcif_category(category, {key: [row[key] for row in records] for key in records[0]})


def canonical_sequence(fasta):
    header, *sequence_lines = fasta.decode().splitlines()
    sequence = "".join(sequence_lines)
    require(header.startswith(">sp|P04626|") and " SV=1" in header, "Expected P04626 sequence version 1")
    require(len(sequence) == 1255, "Unexpected P04626 length")
    require(set(sequence) <= set("ACDEFGHIKLMNPQRSTVWY"), "Invalid P04626 sequence")
    return sequence


def validate_sifts(rcsb, pdbe, canonical):
    identification = rcsb["rcsb_polymer_entity_container_identifiers"]
    require(identification["entity_id"] == "3" and identification["asym_ids"] == ["C"]
            and identification["auth_asym_ids"] == ["C"], "Unexpected RCSB entity or chain")
    alignments = [
        item for item in rcsb["rcsb_polymer_entity_align"]
        if item["provenance_source"] == "SIFTS" and item["reference_database_accession"] == "P04626"
    ]
    expected = {"entity_beg_seq_id": 1, "length": DEPOSITED_END, "ref_beg_seq_id": OFFSET + 1}
    require(len(alignments) == 1 and alignments[0]["aligned_regions"] == [expected], "RCSB SIFTS changed")
    mapping, = pdbe["1n8z"]["UniProt"]["P04626"]["mappings"]
    require((mapping["entity_id"], mapping["chain_id"], mapping["struct_asym_id"]) == (3, "C", "C")
            and (mapping["unp_start"], mapping["unp_end"]) == (23, 629)
            and mapping["identity"] == 1, "PDBe SIFTS disagrees with RCSB")
    for key, number in (("start", 1), ("end", DEPOSITED_END)):
        endpoint = mapping[key]
        require(endpoint == {"author_residue_number": number, "author_insertion_code": "",
                             "residue_number": number}, "Unexpected PDBe SIFTS endpoint")
    require(rcsb["entity_poly"]["pdbx_seq_one_letter_code_can"] == canonical[OFFSET:629],
            "RCSB polymer sequence differs from canonical")
    return expected, mapping


def target_sequence(canonical, bounds):
    low, high = bounds
    return canonical[low + OFFSET - 1:high + OFFSET]


def protein_atoms(block, chain, entity=None):
    return [
        row for row in rows(block, "_atom_site.")
        if row["group_PDB"] == "ATOM" and row["auth_asym_id"] == chain
        and row["pdbx_PDB_model_num"] == "1" and (entity is None or row["label_entity_id"] == entity)
    ]


def residue_aa(name):
    return gemmi.find_tabulated_residue(name).one_letter_code


def xyz(atom):
    return tuple(float(atom[f"Cartn_{axis}"]) for axis in "xyz")


def atom_coordinates(atoms):
    result = {
        (int(atom["auth_seq_id"]), atom["label_comp_id"], atom["label_atom_id"], atom["label_alt_id"]): xyz(atom)
        for atom in atoms
    }
    require(len(result) == len(atoms), "Duplicate atom in coordinate audit")
    return result


def residue_numbers(atoms):
    return {int(atom["auth_seq_id"]) for atom in atoms}


def footprint(target, antibody):
    selected = [atom for atom in target if atom["type_symbol"] not in ("H", "D")]
    distances = cdist([xyz(atom) for atom in selected], [
        xyz(atom) for atom in antibody if atom["type_symbol"] not in ("H", "D")
    ])
    return sorted({int(atom["auth_seq_id"]) for atom, nearest in zip(selected, distances.min(axis=1))
                   if nearest < 4.5})


def audit_coordinates(full, organizer, legacy, canonical, local_pdb):
    atoms = protein_atoms(full, "C", "3")
    resolved = residue_numbers(atoms)
    missing = {int(row["auth_seq_id"]) for row in rows(full, "_pdbx_unobs_or_zero_occ_residues.")
               if row["auth_asym_id"] == "C"}
    require(set(range(1, DEPOSITED_END + 1)) - resolved == missing, "Deposited gap audit failed")
    require(all(atom["auth_seq_id"] == atom["label_seq_id"]
                and atom["label_asym_id"] == "C" and not atom["pdbx_PDB_ins_code"]
                and residue_aa(atom["auth_comp_id"]) == canonical[int(atom["auth_seq_id"]) + OFFSET - 1]
                for atom in atoms), "Deposited coordinate mapping differs from SIFTS")
    if local_pdb is not None:
        chain = gemmi.read_structure(str(local_pdb))[0]["C"]
        local = {(res.seqid.num, res.name, atom.name, atom.altloc.strip("\x00 ") or False): tuple(atom.pos)
                 for res in chain if res.het_flag == "A" for atom in res}
        require(local == atom_coordinates(atoms), "Local PDB polymer coordinates differ from RCSB mmCIF")

    crop_atoms = []
    for chain, (low, high) in CROPS.items():
        cropped = protein_atoms(organizer, chain)
        deposited = [atom for atom in atoms if low <= int(atom["auth_seq_id"]) <= high]
        require(residue_numbers(cropped) == set(range(low, high + 1)), f"Crop {chain} numbering changed")
        require(atom_coordinates(cropped) == atom_coordinates(deposited), f"Crop {chain} coordinates changed")
        crop_atoms.extend(cropped)

    full_contacts = footprint(atoms, protein_atoms(full, "A") + protein_atoms(full, "B"))
    crop_contacts = footprint(crop_atoms, protein_atoms(organizer, "A") + protein_atoms(organizer, "B"))
    require(full_contacts == crop_contacts, "Full Fab and crop Fv contact footprints differ")

    historical = atom_coordinates(protein_atoms(legacy, "C") + protein_atoms(legacy, "D"))
    crop = atom_coordinates(crop_atoms)
    common = sorted(crop.keys() & historical.keys())
    x, y = (np.array([coordinates[key] for key in common]) for coordinates in (crop, historical))
    x -= x.mean(axis=0)
    y -= y.mean(axis=0)
    u, _, v = np.linalg.svd(x.T @ y)
    rotation = u @ np.diag([1, 1, np.linalg.det(u @ v)]) @ v
    displacement = np.linalg.norm(x @ rotation - y, axis=1)
    return atoms, {
        "rcsb_target_polymer_atoms": len(atoms),
        "rcsb_resolved_target_residues": len(resolved),
        "rcsb_unobserved_deposited_mature_positions": sorted(missing),
        "rcsb_unobserved_deposited_canonical_positions": [r + OFFSET for r in sorted(missing)],
        "organizer_crop_atom_count_including_oxt": len(crop),
        "historical_crop_common_atoms": len(common),
        "historical_crop_aligned_rmsd_angstrom": float(np.sqrt(np.mean(displacement**2))),
        "local_pdb_atom_coordinates_verified": local_pdb is not None,
        "footprint_cutoff_angstrom_exclusive": 4.5,
        "footprint_mature_positions": full_contacts,
        "footprint_canonical_positions": [r + OFFSET for r in full_contacts],
    }


def make_template(full, atoms, chains, canonical):
    """Emit query-shaped polymer tables and keep only deposited, non-OXT protein atoms."""
    doc = gemmi.cif.Document()
    block = doc.add_new_block("1N8Z")
    for category in ("_entry.", "_exptl.", "_pdbx_audit_revision_history.", "_refine."):
        block.set_mmcif_category(category, full.get_mmcif_category(category))
    entity, polymer, asym, seq, scheme, coordinates, unobserved = ([] for _ in range(7))
    query_location = {}
    for entity_number, (chain, bounds) in enumerate(chains.items(), 1):
        entity_id = str(entity_number)
        first, last = bounds
        sequence = target_sequence(canonical, bounds)
        entity.append(dict(id=entity_id, type="polymer", pdbx_description="HER2 target"))
        polymer.append(dict(entity_id=entity_id, type="polypeptide(L)", nstd_linkage="no",
                            nstd_monomer="no", pdbx_seq_one_letter_code=sequence,
                            pdbx_seq_one_letter_code_can=sequence, pdbx_strand_id=chain))
        asym.append(dict(id=chain, entity_id=entity_id))
        selected = [atom for atom in atoms if first <= int(atom["auth_seq_id"]) <= last
                    and atom["label_atom_id"] != "OXT"]
        resolved = residue_numbers(selected)
        for index, aa in enumerate(sequence, 1):
            mature = first + index - 1
            query_location[mature] = (chain, str(index))
            monomer = gemmi.expand_one_letter_sequence(aa, gemmi.ResidueKind.AA)[0]
            seq.append(dict(entity_id=entity_id, num=str(index), mon_id=monomer, hetero="n"))
            auth_num = str(mature) if mature in resolved else None
            auth_name = monomer if mature in resolved else None
            scheme.append(dict(asym_id=chain, entity_id=entity_id, seq_id=str(index), mon_id=monomer,
                               ndb_seq_num=str(index), pdb_seq_num=str(mature), auth_seq_num=auth_num,
                               pdb_mon_id=auth_name, auth_mon_id=auth_name, pdb_strand_id=chain,
                               pdb_ins_code=False, hetero="n"))
            if mature not in resolved and mature <= DEPOSITED_END:
                unobserved.append(dict(id=str(len(unobserved) + 1), PDB_model_num="1", polymer_flag="Y",
                                       occupancy_flag="1", auth_asym_id=chain, auth_comp_id=monomer,
                                       auth_seq_id=str(mature), PDB_ins_code=None, label_asym_id=chain,
                                       label_comp_id=monomer, label_seq_id=str(index)))
        for atom in selected:
            coordinates.append(dict(atom, label_asym_id=chain, auth_asym_id=chain, label_entity_id=entity_id,
                                    label_seq_id=str(int(atom["auth_seq_id"]) - first + 1)))
    components = {row["mon_id"] for row in seq}
    for category, records in (
        ("_entity.", entity), ("_entity_poly.", polymer), ("_struct_asym.", asym),
        ("_entity_poly_seq.", seq), ("_pdbx_poly_seq_scheme.", scheme),
        ("_chem_comp.", [row for row in rows(full, "_chem_comp.") if row["id"] in components]),
        ("_atom_site.", coordinates),
    ):
        write_rows(block, category, records)
    if unobserved:
        write_rows(block, "_pdbx_unobs_or_zero_occ_residues.", unobserved)
    disulfides = []
    for row in rows(full, "_struct_conn."):
        if row["conn_type_id"] != "disulf" or row["ptnr1_label_asym_id"] != row["ptnr2_label_asym_id"]:
            continue
        if row["ptnr1_label_asym_id"] != "C":
            continue
        left, right = (query_location.get(int(row[f"ptnr{i}_auth_seq_id"])) for i in (1, 2))
        if left is None or right is None or left[0] != right[0]:
            continue
        for i, (chain, label_id) in enumerate((left, right), 1):
            row.update({f"ptnr{i}_label_asym_id": chain, f"ptnr{i}_auth_asym_id": chain,
                        f"ptnr{i}_label_seq_id": label_id})
        disulfides.append(row)
    if disulfides:
        write_rows(block, "_struct_conn.", disulfides)
    return doc.as_string()


def chain_sequences(query):
    return {chain["chain_id"]: chain["sequence"] for chain in query["chains"]}


def antibody_panel(historical_queries, facs, spr, identity, ranking, scores, optimizer_ids):
    queries = {query["query_id"]: chain_sequences(query) for query in historical_queries}
    wt = queries["control_wt"]

    with facs.open(newline="") as handle:
        low = next(row for row in csv.DictReader(handle)
                   if row["binding_class"] == "low" and "C" not in row["cdrh3"])
    require(wt["H"].count(WT_CDRH3) == 1 and wt["L"] == low["vl_sequence"], "Unexpected WT background")
    require((low["vh_sequence"], low["vl_sequence"]) ==
            (queries["control_facs_low"]["H"], queries["control_facs_low"]["L"]), "FACS-low changed")
    with spr.open(newline="") as handle:
        measured = []
        for row in csv.DictReader(handle):
            if (row["binding_result"] != "quantified_binding" or not row["kd_nM"]
                    or row["vl_sequence"] != wt["L"]
                    or row["vh_sequence"] != wt["H"].replace(WT_CDRH3, row["cdrh3"])):
                continue
            kd = float(row["kd_nM"])
            require(math.isfinite(kd) and kd > 0, "Invalid measured KD")
            measured.append((kd, row))
    wt_kds = [kd for kd, row in measured if row["cdrh3"] == WT_CDRH3]
    require(bool(wt_kds), "The measured WT reference is missing")
    variants = [(kd, row) for kd, row in measured if row["cdrh3"] != WT_CDRH3]
    ordered = sorted(variants, key=lambda pair: (pair[0], pair[1]["cdrh3"]))

    def pair(query_id, row, **details):
        vh, vl = row["vh_sequence"], row["vl_sequence"]
        require(all(sequence and set(sequence) <= set("ACDEFGHIKLMNPQRSTVWY") for sequence in (vh, vl)),
                f"Invalid antibody sequences for {query_id}")
        return dict(query_id=query_id, vh_sequence=vh, vl_sequence=vl,
                    vh_sha256=sha256(vh.encode()), vl_sha256=sha256(vl.encode()), **details)

    panel = [pair("control_wt", dict(vh_sequence=wt["H"], vl_sequence=wt["L"]),
                  cdrh3=WT_CDRH3, kd_nM=wt_kds[0], all_wt_kd_nM=wt_kds),
             pair("control_facs_low", low, cdrh3=low["cdrh3"], facs_class="low")]
    for identifier, (kd, row) in zip(("control_spr_strong", "control_spr_weak"), (ordered[0], ordered[-1])):
        panel.append(pair(identifier, row, cdrh3=row["cdrh3"], kd_nM=kd))

    with identity.open(newline="") as handle:
        candidate, = [row for row in csv.DictReader(handle) if row["arm"] == "denovo-fullseq-k1"
                      and row["design"] == "her2_denovo_full_sequence_1"]
    require((candidate["vh_sequence"], candidate["vl_sequence"]) ==
            (queries[DE_NOVO]["H"], queries[DE_NOVO]["L"]), "De novo query and identity record differ")
    require(candidate["canonical_cysteines_only"] == candidate["below_70_both_local_references"] == "True",
            "De novo local checks failed")
    require(int(candidate["paired_length"]) == len(candidate["vh_sequence"] + candidate["vl_sequence"]) <= 250
            and min(float(candidate[key]) for key in ("vh_confidence", "vl_confidence")) >= 0.1,
            "De novo length or numbering-confidence check failed")
    with scores.open(newline="") as handle:
        geometry, = [row for row in csv.DictReader(handle) if row["query_id"] == DE_NOVO]
    panel.append(pair(DE_NOVO, candidate, identity_design=candidate["design"],
                      design_sha256=candidate["source_sha256"], official_novelty=candidate["official_novelty"],
                      local_paired_imgt_identity={ref: float(candidate[f"{ref}_paired_imgt_identity"])
                                                  for ref in ("trastuzumab", "adalimumab")},
                      historical_footprint_recall=float(geometry["wt_footprint_recall"]),
                      historical_footprint_jaccard=float(geometry["wt_footprint_jaccard"])))
    with ranking.open(newline="") as handle:
        ranked = sorted(csv.DictReader(handle), key=lambda row: int(row["shortlist_pick"]))
    if optimizer_ids:
        by_cdr = {row["cdrh3"]: row for row in ranked}
        ranked = [by_cdr[cdr] for cdr in optimizer_ids]
    for row in ranked[:2]:
        require(row["strict_local_gate"] == "True" and row["vl_sequence"] == wt["L"]
                and wt["H"].replace(WT_CDRH3, row["cdrh3"]) == row["vh_sequence"], "Unexpected optimizer")
        panel.append(pair(f"opt_{row['cdrh3']}", row, cdrh3=row["cdrh3"],
                          shortlist_pick=int(row["shortlist_pick"]), official_validation=row["official_validation"]))
    require(len(panel) == len({row["query_id"] for row in panel}) == 7, "Expected seven unique queries")
    return panel, len(variants)


def json_text(data):
    return json.dumps(data, indent=2) + "\n"


def prepare(args):
    sources = {}
    for key, path, url in (
        ("p04626_fasta", args.uniprot_fasta, P04626), ("pdb_cif", args.pdb_cif, PDB_CIF),
        ("rcsb_entity", args.rcsb_entity, RCSB), ("pdbe_mapping", args.pdbe_mapping, PDBE),
    ):
        data = read_source(path, url)
        sources[key] = (data, url)
    canonical = canonical_sequence(sources["p04626_fasta"][0])
    rcsb_map, pdbe_map = validate_sifts(json.loads(sources["rcsb_entity"][0]),
                                        json.loads(sources["pdbe_mapping"][0]), canonical)
    full = gemmi.cif.read_string(sources["pdb_cif"][0]).sole_block()
    organizer = gemmi.cif.read(str(args.organizer_cif)).sole_block()
    historical = gemmi.cif.read(str(args.historical_template)).sole_block()
    reference_queries = json.loads(args.historical_queries.read_text())["queries"]
    atoms, audit = audit_coordinates(full, organizer, historical, canonical, args.local_pdb)
    deposited = [row for row in rows(full, "_entity_poly_seq.") if row["entity_id"] == "3"]
    require([int(row["num"]) for row in deposited] == list(range(1, DEPOSITED_END + 1))
            and "".join(residue_aa(row["mon_id"]) for row in deposited) == canonical[OFFSET:629],
            "Deposited polymer metadata differs from P04626")
    for query in reference_queries:
        require(all(chain_sequences(query)[chain] == target_sequence(canonical, bounds)
                    for chain, bounds in CROPS.items()), "Historical query crop sequence differs")
    for chain, bounds in CROPS.items():
        actual = "".join(residue_aa(atom["label_comp_id"]) for atom in protein_atoms(historical, chain)
                         if atom["label_atom_id"] == "CA")
        require(actual == target_sequence(canonical, bounds), "Historical template crop sequence differs")

    panel, variant_count = antibody_panel(reference_queries, args.facs, args.spr, args.identity,
                                           args.ranking, args.scores, args.optimizer_cdrh3)
    audit["design_hotspots"] = []
    design = json.loads(DESIGN_SPEC.read_text())[HOTSPOT_SPEC_KEY]
    for item in design["conditions"]["C_HOT"]["true"]:
        selector = item["select"]
        chain, _, residue, atom_name = selector.removeprefix("sel('").removesuffix("')").split("/")
        mature = int(residue)
        require(any(row["auth_seq_id"] == residue and row["auth_atom_id"] == atom_name
                    for row in protein_atoms(organizer, chain)), f"Missing design hotspot {selector}")
        audit["design_hotspots"].append(dict(selector=selector, mature_position=mature,
                                             canonical_position=mature + OFFSET,
                                             amino_acid=canonical[mature + OFFSET - 1],
                                             domain_iv_query_position=mature - 488, ecd_query_position=mature))
    audit["deposited_domain_iv_disulfides_canonical"] = [
        [int(row[f"ptnr{i}_auth_seq_id"]) + OFFSET for i in (1, 2)]
        for row in rows(full, "_struct_conn.") if row["conn_type_id"] == "disulf"
        and row["ptnr1_label_asym_id"] == row["ptnr2_label_asym_id"] == "C"
        and min(int(row[f"ptnr{i}_auth_seq_id"]) for i in (1, 2)) >= 489
    ]
    manifest = {
        "sources": {key: {"url": url, "sha256": sha256(data)} for key, (data, url) in sources.items()},
        "sifts": {"rcsb": rcsb_map, "pdbe": pdbe_map, "mature_to_canonical_offset": OFFSET},
        "canonical_sequence_sha256": sha256(canonical.encode()),
        "matched_non_wt_quantified_spr_rows": variant_count,
        "audit": audit,
        "panel": panel,
        "arms": {},
    }
    for path in (args.organizer_cif, args.historical_template, args.historical_queries, args.facs, args.spr,
                 args.identity, args.ranking, args.scores, DESIGN_SPEC, args.local_pdb):
        if path is not None:
            name = str(path.relative_to(ROOT)) if path.is_relative_to(ROOT) else path.name
            manifest["sources"][name] = {"sha256": sha256(path.read_bytes())}

    files = {"P04626.fasta": sources["p04626_fasta"][0].decode()}
    target_fasta = []
    for arm, chains in ARMS.items():
        template_name, query_name = f"{arm}-target.cif", f"{arm}-queries.json"
        files[template_name] = make_template(full, atoms, chains, canonical)
        target_chains = []
        chain_map = {}
        for chain, bounds in chains.items():
            first, last = bounds
            canonical_positions = list(range(first + OFFSET, last + OFFSET + 1))
            sequence = target_sequence(canonical, bounds)
            missing = set(range(first, last + 1)) - residue_numbers(atoms)
            chain_map[chain] = dict(length=len(sequence), sequence_sha256=sha256(sequence.encode()),
                                    query_position_to_canonical=canonical_positions,
                                    unresolved_deposited_mature=sorted(r for r in missing if r <= DEPOSITED_END),
                                    beyond_deposit_mature=sorted(r for r in missing if r > DEPOSITED_END),
                                    cysteine_canonical_positions=[r for r in canonical_positions
                                                                  if canonical[r - 1] == "C"])
            target_chains.append(dict(chain_id=chain, sequence=sequence, msa_path=None))
            target_fasta.append(f">{arm}|{chain}|P04626:{first + OFFSET}-{last + OFFSET}\n"
                                + textwrap.fill(sequence, width=70) + "\n")
        files[query_name] = json_text({"queries": [
            dict(query_id=row["query_id"], chains=target_chains + [
                dict(chain_id=chain, sequence=row[key], msa_path=None)
                for chain, key in (("H", "vh_sequence"), ("L", "vl_sequence"))
            ]) for row in panel
        ]})
        manifest["arms"][arm] = dict(template=template_name, queries=query_name, target_chains=chain_map)
    files["targets.fasta"] = "".join(target_fasta)
    manifest["artifact_sha256"] = {name: sha256(body.encode()) for name, body in files.items()}
    files["manifest.json"] = json_text(manifest)
    args.output.mkdir(parents=True, exist_ok=True)
    for name, body in files.items():
        (args.output / name).write_text(body)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--facs", type=Path, required=True)
    parser.add_argument("--spr", type=Path, required=True)
    parser.add_argument("--output", type=Path, default=CONFIG)
    parser.add_argument("--uniprot-fasta", type=Path, default=CONFIG / "P04626.fasta")
    parser.add_argument("--pdb-cif", type=Path, help="Local official 1N8Z mmCIF, otherwise download from RCSB")
    parser.add_argument("--rcsb-entity", type=Path, help="Local official RCSB entity JSON")
    parser.add_argument("--pdbe-mapping", type=Path, help="Local official PDBe SIFTS JSON")
    parser.add_argument("--local-pdb", type=Path, help="Independently verify deposited chain C atom coordinates")
    parser.add_argument("--organizer-cif", type=Path, default=CAMPAIGN / "her2_fv_1n8z.cif")
    parser.add_argument("--historical-template", type=Path, default=RESULTS / "input-specs/her2_target_refold.cif")
    parser.add_argument("--historical-queries", type=Path, default=RESULTS / "denovo-refold-queries.json")
    parser.add_argument("--identity", type=Path, default=RESULTS / "denovo_identity.csv")
    parser.add_argument("--ranking", type=Path, default=ROOT / "results/spr-candidate-ranking/shortlist.csv")
    parser.add_argument("--scores", type=Path, default=ROOT / "results/denovo-refold/scores.csv")
    parser.add_argument("--optimizer-cdrh3", nargs=2, metavar=("FIRST", "SECOND"))
    prepare(parser.parse_args())


if __name__ == "__main__":
    main()
