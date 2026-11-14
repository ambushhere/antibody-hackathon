#!/usr/bin/env python3
"""Audit the pinned FLAb release for paired antibody thermal training data."""

import argparse
import csv
import hashlib
import json
import math
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path

from rapidfuzz.distance import Levenshtein

FLAB_REVISION = "e05e3f0135783ffb7de412df6a1467e7a0031911"
ABSCI_REVISION = "cd2bd4b6932b5301d435555dc22084d407e4b92b"
AMINO_ACIDS = frozenset("ACDEFGHIKLMNPQRSTVWY")
TARGET_NAMES = ("trastuzumab", "herceptin", "4d5")
WT_VH = (
    "EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKGRFTISAD"
    "TSKNTAYLQMNSLRAEDTAVYYCSRWGGDGFYAMDYWGQGTLVTVSS"
)
WT_VL = (
    "DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRFSGSRSGTDF"
    "TLTISSLQPEDFATYYCQQHYTTPPTFGQGTKVEIK"
)
ABSCI_COLUMNS = {
    "first_inflection": "DSF - IP #1 (oC)",
    "second_inflection": "DSF - IP #2 (oC)",
    "pdi": "DLS - Cumulant PDI",
}


@dataclass(frozen=True)
class Table:
    source: str
    construct: str
    assay: str
    endpoint: str
    unit: str
    column: str
    family: str | None = None

    @property
    def paired(self):
        return self.construct not in {"single domain", "single domain collection"}


# Measured constructs and endpoints come from the source methods, not FLAb's
# description of which domains it stores. Unverified metadata stay explicit.
TABLES = {
    "garbinski2023_tm1": Table(
        "Garbinski", "unknown", "nanoDSF per FLAb", "Tm1", "unverified", "Tm1 (nanoDSF)"
    ),
    "jain2017biophysical_Tm": Table(
        "Jain 2017", "IgG1", "dye DSF", "Fab-assigned Tm", "degC", "Fab Tm by DSF (°C)"
    ),
    "jain2024assessment_Tm": Table(
        "Jain 2024", "IgG1 or Fab unresolved", "dye DSF", "Fab-assigned Tm",
        "degC", "Fab Tm by DSF (°C)"
    ),
    "sulea2019assisted_sdAb_tm": Table(
        "Sulea", "single domain", "unverified", "Tm", "degC", "Tm"
    ),
    "tresanco2023nbthermo_tm": Table(
        "NbThermo compilation", "single domain collection", "row-specific", "Tm",
        "degC", "Tm measurement"
    ),
}
for group, family in (
    ("C143", "C143"), ("MEDIUCA", "MEDI8852"), ("MEDI", "MEDI8852"),
    ("REGN10987", "REGN10987"), ("S309", "S309"), ("mAb114UCA", "mAb114"),
    ("mAb114", "mAb114"),
):
    TABLES[f"hie2023efficient_{group}_Tm"] = Table(
        "Hie", "purified Fab", "nanoDSF fluorescence ratio", "reported Tm, index unstated",
        "degC", "Tm", family
    )
for group, column in (
    ("multi", "Tm1 (°C) IgG Thermostability"),
    ("single", "Tm1 (C) IgG Thermostability"),
    ("top200", "Tm1 (C) Thermostability (IgG)"),
    ("top27", "Tm1 (C) IgG Thermostability"),
):
    TABLES[f"hutchinson2023enhancement_{group}tm1_igg"] = Table(
        "Hutchinson", "IgG", "intrinsic DSF barycentric mean", "Tm1", "degC", column, "D44.1"
    )
for family in ("adalimumab", "golimumab"):
    TABLES[f"rosace2023automated_tm1_{family}"] = Table(
        "Rosace", "scFv", "nanoDSF 350 nm", "apparent Tm, replicate 1", "degC",
        "Tm (C) rep1", family
    )
for suffix, endpoint in (
    ("DLS", "pdi"), ("DSF1", "first_inflection"), ("DSF2", "second_inflection"),
):
    TABLES[f"shanehsazzadeh2023unlocking_{suffix}"] = Table(
        "Shanehsazzadeh", "IgG1", "DLS" if suffix == "DLS" else "nanoDSF",
        endpoint, "dimensionless" if suffix == "DLS" else "degC",
        ABSCI_COLUMNS[endpoint], "trastuzumab/4D5"
    )


@dataclass(frozen=True)
class Observation:
    filename: str
    line: int
    pair: tuple[str, str]
    value: float | None
    assay: str
    target_reason: str | None

    @property
    def table(self):
        return TABLES[Path(self.filename).stem]

    @property
    def ref(self):
        return f"{self.filename}:{self.line}"

    @property
    def stratum(self):
        table = self.table
        return table.source, table.construct, self.assay, table.endpoint, table.unit


def file_hash(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_csv(path, required):
    # ASVS 1.5.2, 2.2.1 and 5.2.1: parse bounded UTF-8 CSV using a checked schema.
    if path.suffix != ".csv" or not path.is_file() or path.stat().st_size > 2_000_000:
        raise ValueError(f"{path}: expected a regular CSV file up to 2 MB")
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, strict=True)
        columns = reader.fieldnames or []
        if set(required) - set(columns) or len(set(columns)) != len(columns):
            raise ValueError(f"{path.name}: missing required columns or duplicate headers")
        for row in reader:
            if None in row or None in row.values():
                raise ValueError(f"{path.name}:{reader.line_num}: wrong column count")
            yield reader.line_num, row


def parse_number(value):
    try:
        number = float(value.strip().rstrip("\u200b").strip())
    except ValueError:
        return None
    return number if math.isfinite(number) else None


def target_reason(pair, lineage):
    if any(alias in lineage.casefold() for alias in TARGET_NAMES):
        return "known study lineage"
    vh, vl = pair
    vh_similarity = Levenshtein.normalized_similarity(vh, WT_VH)
    vl_similarity = Levenshtein.normalized_similarity(vl, WT_VL)
    if pair == (WT_VH, WT_VL):
        return "exact WT pair"
    if vl == WT_VL or vh_similarity >= 0.95 or min(vh_similarity, vl_similarity) >= 0.90:
        return "conservative WT similarity flag"
    return None


def read_table(path, table):
    required = {"heavy", table.column}
    if table.paired:
        required.add("light")
    if table.assay == "row-specific":
        required.update(("Tm Type", "Origin"))
    valid, invalid = [], []
    n = 0
    for line, row in read_csv(path, required):
        n += 1
        vh, vl = (row.get(chain, "").strip().upper() for chain in ("heavy", "light"))
        if not vh or bool(vl) != table.paired or set(vh + vl) - AMINO_ACIDS:
            invalid.append(line)
            continue
        pair = vh, vl
        valid.append(Observation(
            path.name, line, pair, parse_number(row[table.column]),
            row.get("Tm Type", table.assay),
            target_reason(pair, " ".join((table.family or "", row.get("Antibody", "")))),
        ))
    return n, valid, invalid


def pair_key(pair):
    return hashlib.sha256("|".join(pair).encode("ascii")).hexdigest()


def similarity(pair, other):
    """Return paired global edit similarity weighted by each domain's max length."""
    distance = sum(Levenshtein.distance(a, b) for a, b in zip(pair, other, strict=True))
    length = sum(max(len(a), len(b)) for a, b in zip(pair, other, strict=True))
    return 1 - distance / length


def duplicate_inventory(observations):
    by_pair = defaultdict(list)
    by_stratum = defaultdict(list)
    for row in observations:
        by_pair[row.pair].append(row)
        by_stratum[row.stratum, row.pair].append(row)

    duplicates, shared_pairs = [], set()
    for pair, rows in by_pair.items():
        if len(rows) > 1:
            cross_source = len({row.table.source for row in rows}) > 1
            if cross_source:
                shared_pairs.add(pair)
            duplicates.append({
                "pair_sha256": pair_key(pair),
                "cross_source": cross_source,
                "observations": [{"row": row.ref, "value": row.value} for row in rows],
            })

    # Collapse redundant copies only within an already separated source and
    # assay stratum. Conflicting measurements require the original row mapping.
    unique, conflicts = [], []
    for rows in by_stratum.values():
        measured = [row for row in rows if row.value is not None]
        if len({row.value for row in measured}) > 1:
            conflicts.append([row.ref for row in measured])
        elif measured:
            unique.append(measured[0])
    return duplicates, unique, conflicts, shared_pairs


def link_absci(path, observations):
    authors = {}
    for _, row in read_csv(path, {"Candidate", "HCDR3", *ABSCI_COLUMNS.values()}):
        loop = row["HCDR3"].strip().upper()
        if loop:
            if loop in authors:
                raise ValueError("Absci author table contains duplicate HCDR3s")
            authors[loop] = row

    loop = "SRWGGDGFYAMDY"
    prefix, suffix = WT_VH.split(loop)
    panel = {}
    for row in observations:
        if row.table.source != "Shanehsazzadeh":
            continue
        vh, vl = row.pair
        if vl != WT_VL or not vh.startswith(prefix) or not vh.endswith(suffix):
            raise ValueError(f"{row.ref}: unexpected trastuzumab scaffold")
        author = authors.get(vh[len(prefix):-len(suffix)])
        if not author or row.value is None or row.value != parse_number(author[row.table.column]):
            raise ValueError(f"{row.ref}: measurement does not match the author release")
        case = panel.setdefault(row.pair, {"case": author["Candidate"], "measurements": {}})
        endpoint = row.table.endpoint
        if endpoint in case["measurements"]:
            raise ValueError(f"{row.ref}: duplicate target observation")
        case["measurements"][endpoint] = row.value
    if {case["case"] for case in panel.values()} != {"Trastuzumab", *(str(i) for i in range(1, 13))}:
        raise ValueError("unexpected Absci case coverage")
    if any(set(case["measurements"]) != set(ABSCI_COLUMNS) for case in panel.values()):
        raise ValueError("incomplete Absci endpoint join")
    return panel


def summarize_strata(rows):
    by_stratum = defaultdict(list)
    for row in rows:
        by_stratum[row.stratum].append(row)
    return [
        {
            **dict(zip(("source", "construct", "assay", "endpoint", "unit"), key, strict=True)),
            "distinct_pairs_or_domains": len(group),
            "source_documented_families": dict(sorted(Counter(
                row.table.family for row in group if row.table.family is not None
            ).items())),
            "rows_without_documented_family": sum(row.table.family is None for row in group),
            "value_range": [min(row.value for row in group), max(row.value for row in group)],
        }
        for key, group in sorted(by_stratum.items())
    ]


def audit(flab_dir, absci_csv):
    paths = {f"{stem}.csv": flab_dir / f"{stem}.csv" for stem in sorted(TABLES)}
    found = {path.name for path in flab_dir.glob("*.csv")}
    if found != set(paths):
        raise ValueError(f"Unexpected FLAb table set. Missing {sorted(set(paths) - found)}, extra {sorted(found - set(paths))}")
    observed, tables = [], []
    for name, path in paths.items():
        table = TABLES[path.stem]
        n, rows, invalid = read_table(path, table)
        observed.extend(rows)
        tables.append({
            "file": name,
            "sha256": file_hash(path),
            "source": table.source,
            "construct": table.construct,
            "assay": table.assay,
            "endpoint": table.endpoint,
            "unit": table.unit,
            "label_column": table.column,
            "raw_rows": n,
            "valid_sequence_rows": len(rows),
            "numeric_value_rows": sum(row.value is not None for row in rows),
            "invalid_sequence_lines": invalid,
            "invalid_value_lines": [row.line for row in rows if row.value is None],
            "distinct_pairs_or_domains": len({row.pair for row in rows}),
            "lineage_holdout_rows": sum(row.target_reason is not None for row in rows),
        })
    duplicates, unique, conflicts, shared_pairs = duplicate_inventory(observed)
    non_target = [row for row in unique if not row.target_reason]
    panel = link_absci(absci_csv, observed)

    # Similarity is diagnostic. These pools have no role as training data.
    panels = defaultdict(list)
    for row in non_target:
        if row.table.paired:
            panels[row.table.source].append(row)
    neighbors = []
    for pair, case in panel.items():
        nearest = {}
        for source, rows in sorted(panels.items()):
            best = max(rows, key=lambda row: similarity(pair, row.pair))
            nearest[source] = {
                "row": best.ref,
                "paired_edit_similarity": round(similarity(pair, best.pair), 6),
                "vh_edit_similarity": round(Levenshtein.normalized_similarity(pair[0], best.pair[0]), 6),
                "vl_edit_similarity": round(Levenshtein.normalized_similarity(pair[1], best.pair[1]), 6),
            }
        neighbors.append({
            **case, "pair_sha256": pair_key(pair),
            "primary_supplement_case": case["case"] != "12",
            "nearest_observed_non_target_by_source": nearest,
        })

    return {
        "flab_revision": FLAB_REVISION,
        "absci_revision": ABSCI_REVISION,
        "absci_sha256": file_hash(absci_csv),
        "seed": 0,
        "random_operations": False,
        "model_fitted": False,
        "training_rows": 0,
        "nearest_train_similarity": None,
        "fit_decision": (
            "The scFv and confirmed purified Fab strata each have one independent study. "
            "No assay and construct matched trastuzumab fragment panel is present in FLAb."
        ),
        "raw_tables": tables,
        "raw_rows": sum(table["raw_rows"] for table in tables),
        "paired_rows": sum(row.table.paired for row in observed),
        "paired_tm_or_inflection_rows": sum(
            row.table.paired and row.table.endpoint != "pdi" and row.value is not None for row in observed
        ),
        "unique_valid_sequence_keys": len({row.pair for row in observed}),
        "holdout_observations": [
            {"row": row.ref, "reason": row.target_reason, "endpoint": row.table.endpoint, "value": row.value}
            for row in observed if row.target_reason
        ],
        "cross_source_duplicate_keys": len(shared_pairs),
        "conflicting_same_stratum_observations": conflicts,
        "non_target_strata_after_within_stratum_deduplication": summarize_strata(non_target),
        "non_target_paired_strata_without_cross_source_duplicate_keys": summarize_strata(
            row for row in non_target if row.table.paired and row.pair not in shared_pairs
        ),
        "target_panel": neighbors,
    }, duplicates


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--flab-dir", type=Path, required=True, help="pinned FLAb data/thermostability directory")
    parser.add_argument("--absci-csv", type=Path, required=True, help="pinned Absci author CSV")
    parser.add_argument("--output", type=Path, required=True, help="destination for the two JSON audit files")
    args = parser.parse_args()
    try:
        report, duplicates = audit(args.flab_dir, args.absci_csv)
    except (OSError, ValueError, csv.Error) as error:
        parser.error(str(error))
    args.output.mkdir(parents=True, exist_ok=True)
    for name, content in (("audit.json", report), ("duplicates.json", duplicates)):
        (args.output / name).write_text(json.dumps(content, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(f"Audited {report['raw_rows']} FLAb rows. Model fit skipped. See {args.output / 'audit.json'}")


if __name__ == "__main__":
    main()
