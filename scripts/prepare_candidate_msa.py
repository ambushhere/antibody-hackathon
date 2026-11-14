"""Freeze the target-only MSA variant of the HER2 candidate panel."""

import argparse
import json
from pathlib import Path, PurePosixPath

from scripts.audit_template_controls import MSA_CONFIG, validate_target_msa
from scripts.prepare_candidate_screen import (
    CONFIG as PARENT, CONTROLS, QUERY_FILE, ROOT, TEMPLATE, json_text, panel_bundle, require, sha256,
)

CONFIG = ROOT / "configs/her2-candidate-msa"
ALIGNMENT = "target.a3m"
STAGED_DIR = PurePosixPath("/home/hackathon/almaz/work/candidate-msa-screen-20261114/input")
FALLBACK = "docs/research/candidate-fallback-protocol.md"
PARENT_HASHES = {
    QUERY_FILE: "c35ac8f72b33e20f7b157da82c172e7d0b51ee2590e20cccb198598bd8cc6f14",
    TEMPLATE: "7f0642f0319952d7e78aede7756a808df661c43b947509b32a7415f0b44a2cd7",
    "manifest.json": "0e99a647b10e1c92105b000cf5a1c9c999b86b8aaa4c811ad8591059ee8964b7",
}
CONTROL_SOURCE_HASHES = {
    ALIGNMENT: "dfe58225ad54a75494f9db6e647d8e9c6e9048bb8708addc729a12759a9e63b5",
    QUERY_FILE: "c1dbf19065dfd5e048e7dc15b14885fb8f23aa4607295c47e77a3c33a9b7f21c",
    "manifest.json": "20fd9bc372c8457e65628cd55c08c687e23d452fdce4af0c54662f7203af2300",
}
FALLBACK_HASH = "32d11c44de1df04287b2e5ccf958528416c20a24225ea4a77817ff0f37055319"


def target_msa_control_inputs(template_hash):
    # ASVS 11.4.1: Pin SHA-256 digests independently of the copied alignment manifest.
    for name, digest in CONTROL_SOURCE_HASHES.items():
        require(sha256((MSA_CONFIG / name).read_bytes()) == digest, f"Frozen target MSA source changed: {name}")
    original_controls = json.loads((CONTROLS / QUERY_FILE).read_bytes())
    return validate_target_msa(MSA_CONFIG, original_controls, template_hash)


def msa_bundle():
    parent, canonical = panel_bundle()
    require({name: sha256(payload) for name, payload in parent.items()} == PARENT_HASHES
            and all((PARENT / name).read_bytes() == payload for name, payload in parent.items()),
            "Frozen no-MSA candidate inputs changed")
    require(sha256((ROOT / FALLBACK).read_bytes()) == FALLBACK_HASH, "Frozen MSA fallback protocol changed")
    control_hashes, _ = target_msa_control_inputs(PARENT_HASHES[TEMPLATE])
    queries = json.loads(parent[QUERY_FILE])
    for query in queries["queries"]:
        query["chains"][0]["msa_path"] = str(STAGED_DIR / ALIGNMENT)
    artifacts = {QUERY_FILE: json_text(queries).encode(), TEMPLATE: parent[TEMPLATE],
                 ALIGNMENT: (MSA_CONFIG / ALIGNMENT).read_bytes()}
    manifest = json.loads(parent["manifest.json"])
    manifest["source_sha256"].update({
        FALLBACK: FALLBACK_HASH,
        **{str((PARENT / name).relative_to(ROOT)): digest for name, digest in PARENT_HASHES.items()},
        **{str((MSA_CONFIG / name).relative_to(ROOT)): digest for name, digest in CONTROL_SOURCE_HASHES.items()},
    })
    manifest["parent_input_sha256"] = PARENT_HASHES
    manifest["target_msa_control_input_sha256"] = control_hashes
    manifest["run"].update(msas="target_only", required_control_arm="complete_target_msa",
                           fallback_protocol=FALLBACK, target_msa_path=str(STAGED_DIR / ALIGNMENT))
    manifest["artifact_sha256"] = {name: sha256(payload) for name, payload in artifacts.items()}
    return {**artifacts, "manifest.json": json_text(manifest).encode()}, canonical


def prepare(output):
    artifacts, _ = msa_bundle()
    output.mkdir(parents=True, exist_ok=True)
    for name, payload in artifacts.items():
        path = output / name
        require(not path.is_symlink() and (not path.exists() or path.read_bytes() == payload),
                f"Refusing to replace frozen input {path}")
    for name, payload in artifacts.items():
        path = output / name
        if not path.exists():
            with path.open("xb") as stream:
                stream.write(payload)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=CONFIG)
    args = parser.parse_args()
    prepare(args.output)


if __name__ == "__main__":
    main()
