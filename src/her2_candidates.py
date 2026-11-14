"""Generate and shortlist trastuzumab CDRH3 optimization candidates."""

from __future__ import annotations

import argparse
import csv
import math
import random
import re
import sys
from collections import Counter
from collections.abc import Iterable, Iterator, Mapping, Sequence
from contextlib import contextmanager
from dataclasses import dataclass
from io import TextIOWrapper
from pathlib import Path
from typing import TextIO
from zipfile import BadZipFile, ZipFile


WT_CDRH3 = "WGGDGFYAMD"
CDRH3_LENGTH = len(WT_CDRH3)
AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"
PROPOSAL_ALPHABET = AMINO_ACIDS.replace("C", "")
CANONICAL = frozenset(AMINO_ACIDS)
PROPOSAL_RESIDUES = frozenset(PROPOSAL_ALPHABET)
ORIGINS = ("train_high_mutation", "wt_mutation", "train_high_position_frequency")

# Organizer FACS framework surrounding the reported ten-residue CDRH3.
VH_PREFIX = (
    "EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVK"
    "GRFTISADTSKNTAYLQMNSLRAEDTAVYYCSR"
)
VH_SUFFIX = "YWGQGTLVTVSS"
TRASTUZUMAB_VL = (
    "DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRFSG"
    "SRSGTDFTLTISSLQPEDFATYYCQQHYTTPPTFGQGTKVEIK"
)
CHAIN_COLUMNS = ("vh_sequence", "vl_sequence")
MISSING_EVIDENCE = {
    "structural_status": "missing",
    "thermal_status": "missing",
    "abguard_status": "unverified",
    "de_novo_novelty_status": "unverified",
    "submission_ready": "false",
}
CANDIDATE_COLUMNS = (
    "cdrh3",
    "origin",
    "parent_cdrh3",
    "seed",
    "wt_hamming_distance",
    "novelty_gate",
    "liability_motifs",
    *MISSING_EVIDENCE,
)
LIABILITY_PATTERNS = {
    "n_linked_sequon": re.compile(r"N[^P][ST]"),
    "deamidation_motif": re.compile(r"N[GST]"),
    "isomerization_motif": re.compile(r"D[GS]"),
    "oxidation_susceptible_residue": re.compile(r"[MW]"),
}


def hamming_distance(left: str, right: str) -> int:
    if len(left) != len(right):
        raise ValueError("Hamming distance requires sequences of equal length")
    return sum(a != b for a, b in zip(left, right))


def is_candidate_sequence(sequence: str) -> bool:
    return len(sequence) == CDRH3_LENGTH and PROPOSAL_RESIDUES.issuperset(sequence)


class NoveltyIndex:
    """Exact radius-one Levenshtein lookup for canonical amino acid strings."""

    def __init__(self, sequences: Iterable[str]):
        self.sequences = frozenset(sequences)
        if any(not sequence or not CANONICAL.issuperset(sequence) for sequence in self.sequences):
            raise ValueError("Novelty corpus contains a missing or noncanonical CDRH3")
        self.lengths = {len(sequence) for sequence in self.sequences}

    def has_neighbor(self, sequence: str) -> bool:
        """Detect an exact match, substitution, deletion or insertion."""
        known = self.sequences
        if sequence in known:
            return True
        size = len(sequence)
        if size - 1 in self.lengths and any(
            sequence[:i] + sequence[i + 1 :] in known for i in range(size)
        ):
            return True
        if size in self.lengths and any(
            sequence[:i] + residue + sequence[i + 1 :] in known
            for i in range(size)
            for residue in AMINO_ACIDS
        ):
            return True
        return size + 1 in self.lengths and any(
            sequence[:i] + residue + sequence[i:] in known
            for i in range(size + 1)
            for residue in AMINO_ACIDS
        )


@dataclass(frozen=True)
class TrastuzumabTemplate:
    vh_prefix: str
    vh_suffix: str
    vl_sequence: str

    @classmethod
    def from_facs_row(cls, row: Mapping[str, str]) -> TrastuzumabTemplate | None:
        sequence = row["cdrh3"]
        vh = row["vh_sequence"]
        vl = row["vl_sequence"]
        if len(sequence) != CDRH3_LENGTH or vh != VH_PREFIX + sequence + VH_SUFFIX:
            return None
        if vl != TRASTUZUMAB_VL:
            return None
        start = len(VH_PREFIX)
        return cls(vh[:start], vh[start + CDRH3_LENGTH :], vl)

    def chains(self, sequence: str) -> dict[str, str]:
        if not is_candidate_sequence(sequence):
            raise ValueError("VH export requires a canonical, cysteine-free ten-residue CDRH3")
        return {
            "vh_sequence": self.vh_prefix + sequence + self.vh_suffix,
            "vl_sequence": self.vl_sequence,
        }


@contextmanager
def open_csv(path: Path) -> Iterator[TextIO]:
    if path.suffix.lower() != ".zip":
        with path.open(encoding="utf-8-sig", newline="") as file:
            yield file
        return

    with ZipFile(path) as archive:
        members = [
            entry
            for entry in archive.infolist()
            if entry.filename.lower().endswith(".csv")
            and "__MACOSX" not in Path(entry.filename).parts
            and not Path(entry.filename).name.startswith("._")
        ]
        if len(members) != 1:
            raise ValueError(f"{path}: expected exactly one data CSV in the archive")
        with archive.open(members[0]) as raw, TextIOWrapper(
            raw, encoding="utf-8-sig", newline=""
        ) as file:
            yield file


def csv_rows(path: Path, required: Sequence[str]) -> Iterator[dict[str, str]]:
    # ASVS 2.2.1: enforce the schema before consuming decision inputs.
    with open_csv(path) as file:
        reader = csv.DictReader(file, strict=True)
        columns = reader.fieldnames or []
        missing = set(required).difference(columns)
        if missing or len(columns) != len(set(columns)):
            raise ValueError(f"{path}: expected unique headers including {', '.join(required)}")
        for row in reader:
            if None in row or any(value is None for value in row.values()):
                raise ValueError(f"{path}:{reader.line_num}: wrong number of CSV fields")
            yield row


def load_corpus(
    facs_path: Path, spr_path: Path, *, export_vh_vl: bool = False
) -> tuple[NoveltyIndex, TrastuzumabTemplate | None]:
    """Load every CDRH3 from both tables, including sequences with cysteine."""
    sequences = set()
    template = None
    for path in (facs_path, spr_path):
        find_template = path == facs_path and export_vh_vl
        columns = ("cdrh3", *CHAIN_COLUMNS) if find_template else ("cdrh3",)
        row_count = 0
        for row in csv_rows(path, columns):
            row_count += 1
            sequences.add(row["cdrh3"])
            if find_template and template is None:
                template = TrastuzumabTemplate.from_facs_row(row)
        if row_count == 0:
            raise ValueError(f"{path}: the novelty table is empty")
    if export_vh_vl and template is None:
        raise ValueError("No FACS row matches the fixed trastuzumab VH framework and VL")
    return NoveltyIndex(sequences), template


def load_high_parents(train_path: Path) -> set[str]:
    """Read high labels exclusively from the published training split."""
    return {
        row["seq"]
        for row in csv_rows(train_path, ("seq", "class"))
        if row["class"] == "high" and is_candidate_sequence(row["seq"])
    }


def liability_motifs(sequence: str) -> str:
    return "|".join(name for name, pattern in LIABILITY_PATTERNS.items() if pattern.search(sequence))


def screen_candidate(sequence: str, known: NoveltyIndex, seen: set[str]) -> int | None:
    """Return WT Hamming distance if all sequence gates pass."""
    if not is_candidate_sequence(sequence) or sequence in seen:
        return None
    wt_distance = hamming_distance(sequence, WT_CDRH3)
    if wt_distance < 3 or known.has_neighbor(sequence):
        return None
    return wt_distance


def mutate(
    parent: str, count: int, profile: Sequence[Sequence[int]], rng: random.Random
) -> str:
    residues = list(parent)
    for position in rng.sample(range(CDRH3_LENGTH), count):
        weights = [
            weight if residue != parent[position] else 0
            for residue, weight in zip(PROPOSAL_ALPHABET, profile[position])
        ]
        residues[position] = rng.choices(PROPOSAL_ALPHABET, weights=weights)[0]
    return "".join(residues)


def generate_candidates(
    high_parents: Iterable[str],
    known: NoveltyIndex,
    *,
    pool_size: int = 10_000,
    seed: int = 42,
    max_attempts: int | None = None,
    template: TrastuzumabTemplate | None = None,
) -> list[dict[str, str | int]]:
    if pool_size < 1 or (max_attempts is not None and max_attempts < 1):
        raise ValueError("Pool size and attempt limit must be positive")
    parents = sorted(set(high_parents))
    if not parents or any(not is_candidate_sequence(parent) for parent in parents):
        raise ValueError("Expected canonical, cysteine-free ten-residue training high parents")

    # Add one pseudocount per allowed residue to keep the profile usable on small splits.
    counts = [Counter(position) for position in zip(*parents)]
    profile = [[position[residue] + 1 for residue in PROPOSAL_ALPHABET] for position in counts]
    rng = random.Random(seed)
    candidates = []
    seen: set[str] = set()
    attempt_limit = max_attempts if max_attempts is not None else 20 * pool_size
    for attempt in range(attempt_limit):
        origin = ORIGINS[attempt % len(ORIGINS)]
        if origin == "train_high_mutation":
            parent = rng.choice(parents)
            sequence = mutate(parent, rng.randint(2, 4), profile, rng)
        elif origin == "wt_mutation":
            parent = WT_CDRH3
            sequence = mutate(parent, rng.randint(3, 4), profile, rng)
        else:
            parent = ""
            sequence = "".join(rng.choices(PROPOSAL_ALPHABET, weights=weights)[0] for weights in profile)

        wt_distance = screen_candidate(sequence, known, seen)
        if wt_distance is None:
            continue
        seen.add(sequence)
        candidate = {
            "cdrh3": sequence,
            "origin": origin,
            "parent_cdrh3": parent,
            "seed": seed,
            "wt_hamming_distance": wt_distance,
            "novelty_gate": "facs_spr_edit_distance_at_least_2",
            "liability_motifs": liability_motifs(sequence),
            **MISSING_EVIDENCE,
        }
        if template is not None:
            candidate.update(template.chains(sequence))
        candidates.append(candidate)
        if len(candidates) == pool_size:
            break
    return candidates


def load_scores(path: Path) -> dict[str, float]:
    scores = {}
    seen = set()
    for row in csv_rows(path, ("cdrh3", "predicted_score")):
        sequence = row["cdrh3"]
        if not is_candidate_sequence(sequence) or sequence in seen:
            raise ValueError(f"{path}: invalid or duplicate CDRH3 in scores")
        seen.add(sequence)
        raw_score = row["predicted_score"].strip()
        if raw_score:
            score = float(raw_score)
            if not math.isfinite(score):
                raise ValueError(f"{path}: predicted scores must be finite numbers or blank")
            scores[sequence] = score
    return scores


def greedy_select(
    candidates: Iterable[str],
    scores: Mapping[str, float],
    *,
    count: int,
    min_hamming: int,
) -> list[str]:
    """Select by descending score, keeping input order for tied scores."""
    if count < 1 or not 0 <= min_hamming <= CDRH3_LENGTH:
        raise ValueError(f"Expected count >= 1 and minimum Hamming distance in 0..{CDRH3_LENGTH}")
    pool = list(candidates)
    if any(not is_candidate_sequence(sequence) for sequence in pool) or len(set(pool)) != len(pool):
        raise ValueError("Candidate pool contains an invalid or duplicate CDRH3")
    scored = [sequence for sequence in pool if sequence in scores]
    if any(not math.isfinite(scores[sequence]) for sequence in scored):
        raise ValueError("Predicted scores must be finite")

    selected = []
    for sequence in sorted(scored, key=scores.__getitem__, reverse=True):
        if all(hamming_distance(sequence, earlier) >= min_hamming for earlier in selected):
            selected.append(sequence)
            if len(selected) == count:
                break
    return selected


def write_csv(path: Path, fields: Sequence[str], rows: Iterable[Mapping[str, object]]) -> None:
    def safe_cell(value: object) -> object:
        # ASVS 1.2.10: protect spreadsheet consumers of copied text cells.
        if isinstance(value, str) and (
            value.startswith(("\t", "\r", "\n", "\0"))
            or value.lstrip().startswith(("=", "+", "-", "@"))
        ):
            return "'" + value
        return value

    def write(file: TextIO) -> None:
        writer = csv.writer(file)
        writer.writerow(fields)
        for row in rows:
            writer.writerow(safe_cell(row.get(field, "")) for field in fields)

    if path == Path("-"):
        write(sys.stdout)
    else:
        with path.open("w", encoding="utf-8", newline="") as file:
            write(file)


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    generate = commands.add_parser("generate", help="propose and screen optimization candidates")
    generate.add_argument("--facs", type=Path, required=True, help="full FACS CSV or CSV ZIP")
    generate.add_argument("--spr", type=Path, required=True, help="full SPR CSV")
    generate.add_argument("--train", type=Path, required=True, help="published training CSV with seq,class")
    generate.add_argument("--output", type=Path, required=True, help="candidate CSV, or - for stdout")
    generate.add_argument("--pool-size", type=int, default=10_000)
    generate.add_argument("--seed", type=int, default=42)
    generate.add_argument("--max-attempts", type=int, help="default: 20 times the pool size")
    generate.add_argument("--export-vh-vl", action="store_true", help="export the two variable chains")

    select = commands.add_parser("select", help="shortlist a generated pool using external scores")
    select.add_argument("--candidates", type=Path, required=True, help="candidate CSV from generate")
    select.add_argument("--scores", type=Path, required=True, help="CSV with cdrh3,predicted_score")
    select.add_argument("--output", type=Path, required=True, help="shortlist CSV, or - for stdout")
    select.add_argument("--count", type=int, default=10)
    select.add_argument("--min-hamming", type=int, default=3)

    args = parser.parse_args(argv)
    try:
        if args.command == "generate":
            known, template = load_corpus(args.facs, args.spr, export_vh_vl=args.export_vh_vl)
            parents = load_high_parents(args.train)
            rows = generate_candidates(
                parents,
                known,
                pool_size=args.pool_size,
                seed=args.seed,
                max_attempts=args.max_attempts,
                template=template,
            )
            fields = (*CANDIDATE_COLUMNS, *CHAIN_COLUMNS) if template else CANDIDATE_COLUMNS
            write_csv(args.output, fields, rows)
            print(
                f"Generated {len(rows)}/{args.pool_size} candidates from {len(parents)} eligible "
                f"training high parents, screened against {len(known.sequences)} distinct CDRH3s.",
                file=sys.stderr,
            )
        else:
            pool = list(csv_rows(args.candidates, CANDIDATE_COLUMNS))
            scores = load_scores(args.scores)
            chosen = greedy_select(
                (row["cdrh3"] for row in pool),
                scores,
                count=args.count,
                min_hamming=args.min_hamming,
            )
            by_sequence = {row["cdrh3"]: row for row in pool}
            fields = CANDIDATE_COLUMNS
            if pool and all(column in pool[0] for column in CHAIN_COLUMNS):
                fields = (*fields, *CHAIN_COLUMNS)
            selected = [
                {**by_sequence[sequence], "predicted_score": scores[sequence], "selection_rank": rank}
                for rank, sequence in enumerate(chosen, 1)
            ]
            write_csv(args.output, (*fields, "predicted_score", "selection_rank"), selected)
            print(f"Selected {len(selected)}/{args.count} candidates.", file=sys.stderr)
    except (OSError, ValueError, csv.Error, BadZipFile) as error:
        parser.exit(2, f"{parser.prog}: {error}\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
