"""Run the fixed trastuzumab paired-Fv controls and published Absci panel on CPU."""

import argparse
import csv
import hashlib
import json
import math
import sys
from dataclasses import asdict, dataclass
from pathlib import Path

from Bio.PDB import PDBParser
from Bio.SeqUtils import seq1
from scipy.stats import spearmanr

COMMIT = "df9a75aaddb674a7c4c193005031fc0536d325fb"
TEMPLATES = ("1N8Z", "4HKZ")
VL = "DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRFSGSRSGTDFTLTISSLQPEDFATYYCQQHYTTPPTFGQGTKVEIK"
VH = "EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKGRFTISADTSKNTAYLQMNSLRAEDTAVYYCSRWGGDGFYAMDYWGQGTLVTVSS"
DOMAINS = {"A": VL, "B": VH}
BACKBONE = ("N", "CA", "C", "O")
ALPHABET = "ACDEFGHIKLMNPQRSTVWYX"  # upstream output order
STUDY_LOOP = "SRWGGDGFYAMDY"
LOOP_START = 97  # one-based source VH numbering
PANEL = ("1", "4", "5", "7", "11")  # equal-length entries from Supplementary Table 8
CONTROLS = (
    ("H91Y", 68.5, 72.5, "Morrison Table 1", (("A", 91, "H", "Y"),)),
    ("A34D+Y49S", 68.5, 62.5, "Morrison Table 1", (("A", 34, "A", "D"), ("A", 49, "Y", "S"))),
    ("R66G", 61.8, 66.7, "Okazaki Table 1", (("A", 66, "R", "G"),)),
)


@dataclass(frozen=True)
class Mapping:
    source_chain: str
    source_resnum: int
    source_icode: str
    wt: str
    model_resnum: int
    alternate_backbone: str


def backbone_atom(residue, name):
    atom = residue[name]
    if atom.is_disordered():
        # Highest occupancy, with lexical alternate code as the tie breaker.
        return min(atom.disordered_get_list(), key=lambda a: (-a.occupancy, a.altloc))
    return atom


def prepare_fv(source, destination):
    """Write both variable-domain backbones, retaining chain IDs and coordinates."""
    structure = PDBParser(PERMISSIVE=False, QUIET=True).get_structure(source.stem, source)
    if len(structure) != 1:
        raise ValueError(f"Expected one coordinate model in {source}")
    lines, mapping, ca = [], [], {}
    serial = 1
    for chain, expected in DOMAINS.items():
        residues = [r for r in structure[0][chain] if r.id[0] == " "][:len(expected)]
        sequence = "".join(seq1(r.resname) for r in residues)
        if sequence != expected:
            raise ValueError(f"{source} chain {chain}: variable-domain sequence mismatch")
        for residue in residues:
            _, number, insertion = residue.id
            if residue.is_disordered() == 2 or any(name not in residue for name in BACKBONE):
                raise ValueError(f"{source} {chain}:{number}{insertion}: ambiguous or incomplete backbone")
            position = len(mapping) + 1
            alternates = []
            for name in BACKBONE:
                atom = backbone_atom(residue, name)
                x, y, z = map(float, atom.coord)
                if not all(map(math.isfinite, (x, y, z))):
                    raise ValueError(f"{source} {chain}:{number}{insertion}: nonfinite backbone")
                if atom.altloc != " ":
                    alternates.append(f"{name}:{atom.altloc}({atom.occupancy:.2f})")
                if name == "CA":
                    ca[chain, number, insertion.strip()] = (x, y, z)
                lines.append(
                    f"ATOM  {serial:5d} {atom.fullname} {residue.resname} {chain}{position:4d}    "
                    f"{x:8.3f}{y:8.3f}{z:8.3f}{atom.occupancy:6.2f}{atom.bfactor:6.2f}          {atom.element:>2}"
                )
                serial += 1
            mapping.append(Mapping(chain, number, insertion.strip(), seq1(residue.resname),
                                   position, ";".join(alternates)))
        lines.append(f"TER   {serial:5d}      {residue.resname} {chain}{position:4d}")
        serial += 1
    destination.write_text("\n".join([*lines, "END", ""]))
    return mapping, ca


def lookup_mutations(mapping, substitutions):
    by_source = {(r.source_chain, r.source_resnum, r.source_icode): r for r in mapping}
    mutations = []
    for chain, number, wt, alt in substitutions:
        residue = by_source.get((chain, number, ""))
        if residue is None or residue.wt != wt or alt not in ALPHABET[:20] or wt == alt:
            raise ValueError(f"Invalid substitution or WT mapping: {chain}:{wt}{number}{alt}")
        mutations.append((residue, alt))
    if len({r.model_resnum for r, _ in mutations}) != len(mutations):
        raise ValueError("Repeated position in substitution list")
    return mutations


def read_csv(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def inflection(row, column):
    return float(row[column].replace("\u200b", ""))


def stress_panel(scratch):
    """Join both FLAb transitions by paired sequences and verify author IDs and labels."""
    authors = {row["Candidate"]: row for row in read_csv(scratch / "absci.csv")}
    transitions = {}
    for i in (1, 2):
        column = f"DSF - IP #{i} (oC)"
        rows = read_csv(scratch / f"dsf{i}.csv")
        table = {(row["heavy"], row["light"]): inflection(row, column) for row in rows}
        if len(table) != len(rows):
            raise ValueError(f"Duplicate FLAb sequence pair in transition {i}")
        transitions[i] = table
    if transitions[1].keys() != transitions[2].keys():
        raise ValueError("FLAb paired sequences differ between transitions")
    start = LOOP_START - 1
    end = start + len(STUDY_LOOP)
    if VH[start:end] != STUDY_LOOP or sorted(
        candidate for candidate, row in authors.items()
        if candidate.isdecimal() and 1 <= int(candidate) <= 11 and len(row["HCDR3"]) == len(STUDY_LOOP)
    ) != sorted(PANEL):
        raise ValueError("Unexpected WT loop or equal-length Absci panel")
    panel = {}
    for candidate in ("Trastuzumab", *PANEL):
        author = authors[candidate]
        loop = author["HCDR3"]
        heavy = VH[:start] + loop + VH[end:]
        measured = [transitions[i][heavy, VL] for i in (1, 2)]
        if measured != [inflection(author, f"DSF - IP #{i} (oC)") for i in (1, 2)]:
            raise ValueError(f"FLAb and author values differ for {candidate}")
        substitutions = [
            ("B", position, wt, alt)
            for position, (wt, alt) in enumerate(zip(STUDY_LOOP, loop, strict=True), LOOP_START)
            if wt != alt
        ]
        panel[candidate] = (measured, substitutions)
    if panel["Trastuzumab"][1]:
        raise ValueError("Unexpected substitutions in the published WT")
    return panel


def check_input(pdb, mapping):
    """Ensure the upstream parser sees the full pair in mapping order, without gaps."""
    if pdb["num_of_chains"] != 2 or pdb["seq"] != "".join(DOMAINS.values()):
        raise ValueError("Model parser did not preserve the paired WT sequence")
    for chain, sequence in DOMAINS.items():
        positions = [str(r.model_resnum) for r in mapping if r.source_chain == chain]
        if pdb[f"seq_chain_{chain}"] != sequence or pdb[f"resn_list_{chain}"] != positions:
            raise ValueError(f"Model numbering or sequence differs for chain {chain}")


def single_scores(upstream, pdb, mapping, config, model):
    ddg, sequence = upstream.run_single_ssm(pdb, config, model)
    wt = "".join(ALPHABET[index] for index in sequence.tolist())
    if wt != "".join(r.wt for r in mapping) or ddg.shape != (len(mapping), len(ALPHABET)):
        raise ValueError("Featurized WT sequence or score dimensions differ from mapping")
    scores, labels = upstream.format_output_single(ddg, sequence, threshold=100)
    if len(labels) != len(mapping) * 20 or not ddg.isfinite().all().item():
        raise ValueError("Incomplete or nonfinite upstream substitution matrix")
    return dict(zip(labels, map(float, scores), strict=True))


def joint_score(upstream, torch, pdb, config, model, mutations):
    """Use the upstream double head on the one preregistered pair."""
    pdb["mutation"] = upstream.Mutation([0], ["A"], ["A"], [0.0], "")
    x, seq, mask, _, chain_mask, chain_ids, residue_idx, *_ = upstream.tied_featurize_mut([pdb])
    if not mask.all().item():
        raise ValueError("Missing backbone in joint control")
    x = torch.nan_to_num(x, nan=0.0)
    hidden, embed, _, edges = model.prot_mpnn(x, seq, mask, chain_mask, residue_idx, chain_ids)
    positions = [residue.model_resnum - 1 for residue, _ in mutations]
    wildtypes = [ALPHABET.index(residue.wt) for residue, _ in mutations]
    if seq[0, positions].tolist() != wildtypes:
        raise ValueError("Joint WT sequence differs from mapping")
    alts = [ALPHABET.index(alt) for _, alt in mutations]
    batch = tuple(torch.tensor([values], device="cpu") for values in (positions, wildtypes, alts))
    score = float(upstream.run_double(hidden, embed, config, [batch], 1, model, x, mask, edges))
    if not math.isfinite(score):
        raise ValueError("Nonfinite joint control score")
    return score


def write_csv(path, rows):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def verify_scratch(scratch):
    """Check pinned inputs and patched entry points before importing upstream code."""
    evidence = Path(__file__).resolve().parents[1] / "results" / "thermal-pilot"
    for manifest in ("inputs.sha256", "scratch-derived.sha256"):
        for line in (evidence / manifest).read_text().splitlines():
            expected, relative = line.split(maxsplit=1)
            # Derived backbone files are created later by prepare_fv().
            if manifest == "scratch-derived.sha256" and not relative.endswith(".py"):
                continue
            path = scratch / relative
            if not path.is_file() or hashlib.sha256(path.read_bytes()).hexdigest() != expected:
                raise ValueError(f"Thermal input integrity check failed: {relative}")


def run(args):
    import torch

    scratch = args.scratch.resolve()
    output = args.output.resolve()
    verify_scratch(scratch)
    if output.exists() and any(output.iterdir()):
        raise ValueError("Use a new empty output directory to preserve reference results")
    source = scratch / f"ThermoMPNN-D-{COMMIT}"
    sys.path.insert(0, str(source))
    import v2_ssm as upstream

    torch.set_num_threads(2)
    torch.set_num_interop_threads(1)
    torch.manual_seed(0)
    output.mkdir(parents=True, exist_ok=True)
    panel = stress_panel(scratch)
    modes = ("single", "epistatic") if args.joint else ("single",)
    models = {}
    for mode in modes:
        config = upstream.get_config(mode)
        config.platform.thermompnn_dir = str(source)
        config.platform.accel = "cpu"
        models[mode] = (config, upstream.get_model(mode, config).eval())
        if any(p.device.type != "cpu" for p in models[mode][1].parameters()):
            raise ValueError(f"{mode} model is outside CPU")

    controls, stress, contributions, rank_stats = [], [], [], []
    with torch.inference_mode():
        for template in TEMPLATES:
            fv = scratch / f"{template}-fv.pdb"
            mapping, ca = prepare_fv(scratch / f"{template}.pdb", fv)
            write_csv(output / f"{template}-mapping.csv", [asdict(r) for r in mapping])
            pdb = upstream.load_pdb(str(fv), list(DOMAINS))
            check_input(pdb, mapping)
            scores = single_scores(upstream, pdb, mapping, *models["single"])

            def additive(panel_id, substitutions):
                mutations = lookup_mutations(mapping, substitutions)
                values = []
                for r, alt in mutations:
                    score = scores[f"{r.wt}{r.model_resnum}{alt}"]
                    values.append(score)
                    contributions.append({"template": template, "panel_id": panel_id,
                                          **asdict(r), "mutant": alt, "folding_ddg_kcal_mol": score})
                return math.fsum(values), mutations

            for name, wt_tm, mutant_tm, study, substitutions in CONTROLS:
                score, mutations = additive(name, substitutions)
                delta = round(mutant_tm - wt_tm, 2)
                distance, joint = None, None
                if len(mutations) == 2:
                    distance = math.dist(ca["A", 34, ""], ca["A", 49, ""])
                    if args.joint and 0 < distance <= 5.0:
                        joint = joint_score(upstream, torch, pdb, *models["epistatic"], mutations)
                controls.append({"template": template, "control": name, "study": study,
                                 "wt_tm_c": wt_tm, "mutant_tm_c": mutant_tm, "delta_tm_c": delta,
                                 "folding_ddg_kcal_mol": score, "direction_agrees": score * delta < 0,
                                 "ca_distance_angstrom": distance, "joint_ddg_kcal_mol": joint,
                                 "joint_direction_agrees": None if joint is None else joint * delta < 0})

            template_stress = []
            wt_first, wt_second = panel["Trastuzumab"][0]
            for candidate in PANEL:
                (first, second), substitutions = panel[candidate]
                score, _ = additive(f"Absci {candidate}", substitutions)
                delta = round(first - wt_first, 2)
                template_stress.append({"template": template, "absci_id": candidate,
                                        "substitutions": len(substitutions),
                                        "first_inflection_c": first, "second_inflection_c": second,
                                        "wt_first_inflection_c": wt_first, "wt_second_inflection_c": wt_second,
                                        "delta_first_inflection_c": delta, "folding_ddg_kcal_mol": score,
                                        "direction_agrees": None if candidate == "4" else score * delta < 0})
            stress.extend(template_stress)
            for subset, rows in (("all_five", template_stress),
                                 ("four_clear", [r for r in template_stress if r["absci_id"] != "4"])):
                severity = [-r["delta_first_inflection_c"] for r in rows]
                rank_stats.append({
                    "template": template, "subset": subset,
                    "proxy_spearman": float(spearmanr(severity, [r["folding_ddg_kcal_mol"] for r in rows]).statistic),
                    "count_spearman": float(spearmanr(severity, [r["substitutions"] for r in rows]).statistic),
                })
    write_csv(output / "controls.csv", controls)
    write_csv(output / "absci-stress.csv", stress)
    write_csv(output / "substitution-contributions.csv", contributions)
    (output / "rank-stats.json").write_text(json.dumps(rank_stats, indent=2) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--scratch", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--joint", action="store_true", help="Score the prespecified close double with the joint head")
    run(parser.parse_args())
