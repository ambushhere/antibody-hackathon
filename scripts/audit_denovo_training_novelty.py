"""Compare de novo designs with paired Fv entries in the pinned RFD4 distillation table.

This is a larger local reference than two antibodies. It is still not Adaptyv's database,
and the table describes predicted structures, not measured binders. A design that clears
this check has not been shown to clear the official validator.
"""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from biotite.sequence import ProteinSequence
from biotite.sequence.align import SubstitutionMatrix, align_optimal

PINNED_SHA256 = "d4264ccca50832c8ed9edbe6c4c252de9a017069a71890551f573977b1b45f63"


def variable_domains(annotation):
    """Return {'heavy': (sequence, cdr3), 'light': ...} from one cdr_annotations value."""
    out = {}
    for unit in json.loads(annotation).values():
        if unit.get("ig_type") != "antibody":
            continue
        for domain in unit.get("domains", []):
            kind = domain.get("chain_type")
            if kind not in ("heavy", "light", "kappa", "lambda"):
                continue
            kind = "heavy" if kind == "heavy" else "light"
            out.setdefault(kind, domain)
    return out


def reference_table(path):
    import pyarrow.parquet as pq

    frame = pq.read_table(path, columns=[
        "example_id", "source", "format", "ig_seqs", "ig_chains", "cdr_annotations",
    ]).to_pandas()
    frame = frame[frame.format == "paired_fv"]
    rows = []
    for record in frame.itertuples(index=False):
        domains = variable_domains(record.cdr_annotations)
        if set(domains) != {"heavy", "light"}:
            continue
        sequences = {}
        for kind, domain in domains.items():
            index = list(record.ig_chains).index("heavy" if kind == "heavy" else "light")
            whole = record.ig_seqs[index]
            start, end = domain["fw_start_seq_idx"], domain["fw_end_seq_idx"]
            sequences[kind] = whole[start:end + 1]
        if not all(10 < len(s) < 200 for s in sequences.values()):
            continue
        rows.append({"example_id": record.example_id, "source": record.source,
                     "vh": sequences["heavy"], "vl": sequences["light"]})
    return pd.DataFrame(rows)


def kmers(sequence, k=4):
    return {sequence[i:i + k] for i in range(len(sequence) - k + 1)}


def aligned_identity(left, right, matrix):
    alignment = align_optimal(ProteinSequence(left), ProteinSequence(right),
                              matrix, local=False, max_number=1)[0]
    trace = alignment.trace
    both = (trace[:, 0] >= 0) & (trace[:, 1] >= 0)
    symbols = alignment.sequences
    matches = sum(symbols[0][a] == symbols[1][b] for a, b in trace[both])
    return matches, len(trace)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--metadata", type=Path, required=True,
                        help="pinned distilled_examples.parquet")
    parser.add_argument("--designs", type=Path, required=True,
                        help="identity.csv from audit_denovo_identity")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--shortlist", type=int, default=400,
                        help="entries aligned exactly per design after the k-mer prefilter")
    args = parser.parse_args()

    reference = reference_table(args.metadata)
    designs = pd.read_csv(args.designs)
    designs = designs[
        designs.below_70_both_local_references
        & designs.canonical_cysteines_only
        & (designs.paired_length <= 250)
        & (designs.vh_confidence >= 0.1)
        & (designs.vl_confidence >= 0.1)
    ].reset_index(drop=True)

    matrix = SubstitutionMatrix.std_protein_matrix()
    # Prefilter on both chains: a nearest neighbour can be mediocre on VH and close on VL.
    reference_kmers = [(kmers(h), kmers(l)) for h, l in zip(reference.vh, reference.vl)]
    rows = []
    for record in designs.itertuples(index=False):
        query_h, query_l = kmers(record.vh_sequence), kmers(record.vl_sequence)
        overlap = np.array([
            (len(query_h & other_h) + len(query_l & other_l))
            / max(1, len(query_h | other_h) + len(query_l | other_l))
            for other_h, other_l in reference_kmers])
        order = np.argsort(-overlap)[:args.shortlist]
        best = None
        for index in order:
            entry = reference.iloc[index]
            mh, nh = aligned_identity(record.vh_sequence, entry.vh, matrix)
            ml, nl = aligned_identity(record.vl_sequence, entry.vl, matrix)
            paired = (mh + ml) / (nh + nl)
            if best is None or paired > best["training_paired_identity"]:
                best = {"training_paired_identity": paired,
                        "training_vh_identity": mh / nh,
                        "training_vl_identity": ml / nl,
                        "nearest_example_id": entry.example_id,
                        "nearest_source": entry.source}
        rows.append({"arm": record.arm, "design": record.design,
                     "trastuzumab_paired_imgt_identity": record.trastuzumab_paired_imgt_identity,
                     **best,
                     "below_70_training_reference": best["training_paired_identity"] < .7})
    result = pd.DataFrame(rows).sort_values("training_paired_identity")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(args.output, index=False)
    print(f"reference paired Fv entries: {len(reference)}")
    print(f"designs checked: {len(result)}; below 70%: {int(result.below_70_training_reference.sum())}")
    print(result.groupby("arm").training_paired_identity.describe()[["min", "50%", "max"]].to_string())


if __name__ == "__main__":
    main()
