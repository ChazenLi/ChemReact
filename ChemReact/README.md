# ChemReact

ChemReact is a publishable chemistry skill for synthesis/retrosynthesis workflow automation.
It combines RDKit analysis, route auditing contracts, visualization generation, and markdown reporting
into a single self-contained closed-loop pipeline for host LLM tools (Claude Code, OpenCode, Cursor).

## Features

- End-to-end pipeline from `target_smiles` to report and route visuals.
- Contract-driven JSON validation for:
  - `strategy.json`
  - `routes.json`
  - `vis_plan.json`
- Visualization outputs:
  - target molecule
  - route overview grid
  - orthogonal retrosynthesis tree (`target -> pre1 + pre2`)
  - per-step reaction images
- Host adapters for CC/OpenCode/Cursor.

## Minimal Requirements

- Python 3.10+
- RDKit
- Pillow

Install (pip environment):

```bash
pip install -r skills/ChemReact/requirements.txt
```

If RDKit installation fails with pip on your platform, install RDKit via conda first, then run:

```bash
pip install Pillow
```

## Repository Structure

- `skills/ChemReact/scripts/run_closed_loop.py`: main pipeline entry
- `skills/ChemReact/internal/`: bundled runtime modules (no external folder dependency)
- `skills/ChemReact/adapters/*.py`: host integration adapters
- `skills/ChemReact/schemas/*.json`: data contracts
- `skills/ChemReact/samples/*.json`: runnable sample inputs
- `skills/ChemReact/USAGE.md`: command cookbook
- `skills/ChemReact/SKILL.md`: skill metadata/instructions

## Quick Start

1. Generate schema + visualization template:

```bash
python skills/ChemReact/scripts/run_closed_loop.py \
  --output-dir outputs/run_001 \
  --emit-vis-plan-template-only
```

2. Run full pipeline:

```bash
python skills/ChemReact/scripts/run_closed_loop.py \
  --target-smiles "CCCCC1=NC(Cl)=C(CO)N1CC2=CC=C(C3=CC=CC=C3C4=NNN=N4)C=C2" \
  --routes-file skills/ChemReact/samples/routes.json \
  --strategy-file skills/ChemReact/samples/strategy.json \
  --vis-plan-file skills/ChemReact/samples/vis_plan.json \
  --output-dir outputs/run_001
```

3. Read outputs:

- `outputs/run_001/RETRO_REPORT.md`
- `outputs/run_001/images/*.png`
- `outputs/run_001/run_summary.json`

## Host Integration

Claude Code:

```bash
python skills/ChemReact/adapters/claudecode_adapter.py \
  --request-file skills/ChemReact/samples/host_request.json
```

OpenCode:

```bash
python skills/ChemReact/adapters/opencode_adapter.py \
  --request-file skills/ChemReact/samples/host_request.json
```

Cursor:

```bash
python skills/ChemReact/adapters/cursor_adapter.py \
  --request-file skills/ChemReact/samples/host_request.json
```

## FAQ

Q: Must I provide `vis_plan.json`?
A: No. Omit it and ChemReact falls back to score-based route selection.

Q: Why did validation fail before rendering?
A: One of your input JSON files does not satisfy contract checks. Use `--validate-only` first.

Q: Can this run multi-step retrosynthesis trees?
A: Current format is route + step list (single-track per route). Multi-step DAG support is planned.

Q: Which file is best for downstream automation?
A: `run_summary.json` is the canonical machine-readable output.


## Versioning and Release

Use semantic versioning (`vMAJOR.MINOR.PATCH`).
See `skills/ChemReact/RELEASE.md` for tag/release checklist and commands.


