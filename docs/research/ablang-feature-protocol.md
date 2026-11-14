# Frozen AbLang2 SPR features

14 November 2026. The extraction specification implements the paired representation and single structural ablation in the [interaction protocol](interaction-protocol.md).

## Model and pairing

Use Oxford Protein Informatics Group's AbLang2 package at commit `586af3083c32b5fb2f0a1c855f2ad4f1cad15ec3` (package 0.1.1, BSD 3-Clause). The upstream `ablang2-paired` loader points to Tobias Olsen's Zenodo record 10185169, which labels the weights BSD 3-Clause Clear. The archive contains `model.pt` and `hparams.json`. The extractor verifies the installed Python source tree and both extracted files before a strict, tensor-only CPU state load using the upstream constructor arguments. The checkpoint has 12 encoder blocks, width 480, 20 attention heads and SwiGLU. [Official source and license](https://github.com/oxpig/AbLang2/tree/586af3083c32b5fb2f0a1c855f2ad4f1cad15ec3), [Official loader](https://github.com/oxpig/AbLang2/blob/586af3083c32b5fb2f0a1c855f2ad4f1cad15ec3/ablang2/load_model.py), [Zenodo weights and license](https://zenodo.org/records/10185169), [AbLang2 paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC11543624/)

| Artifact | SHA-256 |
|---|---|
| Installed AbLang2 Python tree | `9fbb5d5cc62e7d9dfb8f99dfa8df81e8bb04779a15c82496e00c570a93cb11c2` |
| AbLang2 source archive | `e83a6f24c23e93dee4e6d5bef2a75d942f05bb0da8fde8035efdd8fdd5ec7c12` |
| Zenodo archive | `7acb23b880e3b69b3fea8b7d6cacac0fda147daafd99edc411e0ff5c61252ac5` |
| `model.pt` | `56d6f07862a6f824f88c8707bbc03e4026c9db762be2d3041e9767e2e6f86386` |
| `hparams.json` | `33474b76d78e01beb24b3621cd44d42b0813333d21e5b806168580bd836c5c8d` |

Read the 158 `Primary novelty cohort` CDRH3 IDs in frozen OOF file order. For each ID, validate the full paired VH and VL against the original SPR sequence columns. Repeated rows may agree on the same pair. Each ID must have exactly one distinct pair, with the loop in the recorded trastuzumab heavy scaffold and with the matched WT VL. The WT pair also has to match the variable-domain prefixes of deposited 1N8Z chains B (VH, 120 aa) and A (VL, 107 aa). Only sequence columns and OOF cohort membership enter extraction. The frozen cohort already reflects the eligibility criteria in the original SPR benchmark.

Use `<VH>|<VL>`, the exact format built by the upstream `pretrained` residue-encoding wrapper for complete chains. Pass that explicit string to `ABtokenizer(..., pad=True, w_extra_tkns=False)`. The separator is token 25, and the start and end tokens are 0 and 22. The 232-element token arrays contain five boundary/separator tokens and 227 residues. No padding is introduced because all contexts have the same length. Verify the decoded pair and H3 slice against the expected sequences before each forward pass. Take `AbRep(tokens).last_hidden_states`, the final normalized hidden state at every token. Use CPU float32 with one intra-op and one inter-op thread, batch size four, seed zero, deterministic algorithms, eval and inference mode. There is no masking, dropout or SPR fitting during extraction. [Official complete-chain formatting](https://github.com/oxpig/AbLang2/blob/586af3083c32b5fb2f0a1c855f2ad4f1cad15ec3/ablang2/pretrained.py), [Tokenizer](https://github.com/oxpig/AbLang2/blob/586af3083c32b5fb2f0a1c855f2ad4f1cad15ec3/ablang2/models/ablang2/tokenizers.py), [Residue output](https://github.com/oxpig/AbLang2/blob/586af3083c32b5fb2f0a1c855f2ad4f1cad15ec3/ablang2/models/ablang2/ablang.py)

## Window mapping and pooling

The ten library positions follow the heavy scaffold's `CSR` and precede `YWGQGTLVTVSS`. The mapping checks the deposited polymer sequence scheme, author IDs, and resolved atoms. WT H3 positions outside window10 stay in the paired model context. In 1N8Z, the target is polymer entity 3, author and label chain C. Its deposited reference alignment and exact polymer sequence map author 1–607 to UniProt P04626 positions 23–629. [1N8Z deposition](https://www.rcsb.org/structure/1N8Z), [Canonical HER2 and SIFTS audit](her2-target-context.md)

| Library position (1-based) | WT | VH residue index (0-based) | AbRep output index (0-based) | 1N8Z heavy label and author residue | P04626 contacts |
|---:|:---:|---:|---:|:---:|---|
| 1 | W | 98 | 99 | B99 | 594, 595 |
| 2 | G | 99 | 100 | B100 | |
| 3 | G | 100 | 101 | B101 | |
| 4 | D | 101 | 102 | B102 | |
| 5 | G | 102 | 103 | B103 | 592, 601, 613, 615 |
| 6 | F | 103 | 104 | B104 | 615 |
| 7 | Y | 104 | 105 | B105 | 592, 593, 594, 595 |
| 8 | A | 105 | 106 | B106 | |
| 9 | M | 106 | 107 | B107 | |
| 10 | D | 107 | 108 | B108 | |

Primary `cdrh3_mean` is the arithmetic mean of the ten contextual window10 vectors. Structural ablation `wt_contact_mean` averages the same vectors at library positions **1, 5, 6 and 7**. A contact requires at least one interchain pair of observed polymer non-hydrogen/non-deuterium atoms with positive occupancy at distance strictly below 4.5 Å, in model 1 of the deposited WT 1N8Z mmCIF. The computation includes all positive occupancy deposited conformers and excludes non-polymer atoms. The extractor stops on ambiguous mapping, missing resolved H3 residues or an empty mask. The mask describes the observed WT geometry. Variant contact retention and affinity require separate evidence.

## Reproduction

Run from the repository root. Use the original frozen SPR CSV (the organizer's `trastuzumab_her2_spr_all_768.csv`, SHA-256 `4d6f53a6183f377670e9577b3233a7c4c14888965ac63824863219d13424fa9b`). The script also checks the OOF, 1N8Z, and canonical FASTA digests. If reproducing an existing feature artifact, use a fresh `--output-dir` as shown below.

```bash
set -e
ABLANG_RUN="$(mktemp -d /tmp/her2-ablang-replay.XXXXXX)"
SPR_CSV=/absolute/path/to/trastuzumab_her2_spr_all_768.csv
python3.12 -m venv "$ABLANG_RUN/venv"
source "$ABLANG_RUN/venv/bin/activate"
# Linux x86_64 CPU. On macOS omit this line and install from requirements directly.
python -m pip install --index-url https://download.pytorch.org/whl/cpu 'torch==2.8.0+cpu'
python -m pip install -r requirements-ablang.txt
curl -fL --retry 2 https://files.rcsb.org/download/1N8Z.cif -o "$ABLANG_RUN/1N8Z.cif"
curl -fL --retry 2 https://zenodo.org/records/10185169/files/ablang2-weights.tar.gz \
  -o "$ABLANG_RUN/ablang2-weights.tar.gz"
(
  cd "$ABLANG_RUN"
  printf '%s  %s\n' 7acb23b880e3b69b3fea8b7d6cacac0fda147daafd99edc411e0ff5c61252ac5 ablang2-weights.tar.gz \
    | shasum -a 256 -c -
)
mkdir "$ABLANG_RUN/weights"
tar -xzf "$ABLANG_RUN/ablang2-weights.tar.gz" -C "$ABLANG_RUN/weights" hparams.json model.pt
PYTHONDONTWRITEBYTECODE=1 python scripts/extract_ablang_features.py \
  --spr "$SPR_CSV" --cif "$ABLANG_RUN/1N8Z.cif" --weights-dir "$ABLANG_RUN/weights" \
  --output-dir "$ABLANG_RUN/features"
```

The output contains `features.npz` and `metadata.json`. Read the NPZ with `allow_pickle=False`. The Unicode `sequence_ids` array uses exact frozen OOF order. Both feature arrays have shape `(158, 480)` and dtype float32. Metadata records checksums, exact sequence inputs, model parameters, residue indices, closest WT heavy-atom distances and runtime versions. Weights and downloaded structure stay outside the repository. Deterministic settings apply to a fixed runtime and platform. Floating-point replay across CPU architectures or build variants can differ slightly.

The local macOS arm64 run encoded all 158 pairs in 8.295 seconds of inference on one CPU thread. The [feature archive](../../results/ablang-features/features.npz) has SHA-256 `13de9d1450fcaf7fa24b2faeadb911b46ea3ce18ec22530ba777d298d1afad27`. The [metadata](../../results/ablang-features/metadata.json) and frozen OOF agree on every ID and its order. Run the tests in the extraction environment with `python -m pip install pytest==8.4.2` and `PYTHONDONTWRITEBYTECODE=1 python -m pytest -p no:cacheprovider -q tests/test_ablang_features.py`. The tests use the tokenizer and synthetic hidden states, with no checkpoint load or network requests.
