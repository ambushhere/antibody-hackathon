# Setup and historical team operations

For selection replay without datasets, credentials or GPU access, start with the [README](../README.md#reproduce-the-final-selection). For archived raw evidence and deleted cache restoration, see [data retention](DATA_RETENTION.md).

## Full development environment

Clone the repository, install [uv](https://docs.astral.sh/uv/getting-started/installation/) and run `bash scripts/bootstrap.sh`. Activate `.venv` before Python commands. The bootstrap creates `.env` from the portable template only if absent.

## Event infrastructure (historical)

The access instructions below describe the 14 November 2026 event. Availability and credentials have not been revalidated after the hackathon. Selection replay does not require these services.

Get the team SSH key and Nebius credentials privately from a teammate. GitHub Actions secrets cannot be downloaded back to a laptop. Set AMASS_API_KEY, JAPANFOLD_API_KEY, NEBIUS_CREDS and VM_SSH_KEY in `.env`; keep credential files mode 0600 and outside Git.

```bash
source .venv/bin/activate
set -a
source .env
set +a
ssh -i "$VM_SSH_KEY" "$VM_HOST"
```

During the event, the VM contained `/home/hackathon/almaz/data/hackathon-data`, organizer docs at `/home/hackathon/almaz/data/hack/docs` and Python environment `/home/hackathon/almaz/work/.venv`. It was scheduled to shut down on 14 November; do not assume it is still available. Coordinate GPU use with the team. Never overwrite another run directory.

## GitHub Actions

`Tests` runs CPU unit tests with no secrets. `Team access check` is manually dispatched from Actions and checks credential presence plus authenticated VM connectivity. It does not start models or spend credits. The workflow expects repository secrets named AMASS_API_KEY, JAPANFOLD_API_KEY, NEBIUS_CREDENTIALS_JSON, VM_HOST, VM_SSH_KEY and VM_KNOWN_HOSTS. The Nebius secret contains API/storage credentials, not the browser login password. SSH host verification is enabled using the previously trusted host key.

## Codex MCP

Node.js/npm and Python are needed. RCSB is a public HTTP MCP; the other servers use the launcher, which reads credentials at runtime and passes them through environment variables, not process arguments. The bridge is pinned to mcp-remote 0.16.6.

```bash
codex mcp add rcsb --url https://mcp-beta.rcsb.org/mcp
codex mcp add amass -- python3 "$PWD/scripts/mcp_launcher.py" amass --env-file "$PWD/.env"
codex mcp add japanfold -- python3 "$PWD/scripts/mcp_launcher.py" japanfold --env-file "$PWD/.env"
codex mcp add nebius-scientific -- python3 "$PWD/scripts/mcp_launcher.py" nebius --credentials "$NEBIUS_CREDS"
```

Start a new Codex session after configuration. On the maintainer machine all four endpoints accepted MCP initialize; Amass, Nebius and JapanFold also passed initialize through the actual stdio launcher. New tools are not dynamically added to an already running chat.

Apheris on the VM is Docker/CLI, not a Foundry Hub endpoint. Use RUNBOOK.md and the organizer `18-vm.md`. Nebius also offers the organizer `forge.py` wrapper; the participant installer installs that into Claude skills, independently of MCP.

## Design campaign

Read [CAMPAIGN.md](CAMPAIGN.md) for rule discrepancies, benchmark design and evidence boundaries. The RFD4 runner records seed/image provenance and refuses existing output folders. `prepare_denovo.py` constructs an exploratory HER2 CDR-redesign input, but database novelty is not established. `her2_structure_audit.py` measures explicit antibody-target geometry and excludes VH–VL contacts.

The installed abguard 0.1.0 package and its bundled reference were later found and archived. All 30 final rows passed that build. See the [validation report](research/abguard-validation.md) and [archive restoration guide](../results/archives/README.md). Installed-tool eligibility does not establish final organizer acceptance or biological performance.
