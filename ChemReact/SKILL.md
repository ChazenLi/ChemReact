---
name: chemreact
description: Integrate RDKit analysis and retrosynthesis visualization into an executable closed loop. Use when the user needs one skill that consumes target SMILES plus route/audit JSON, generates visualization assets and markdown report, and can be plugged into host LLM tools such as Claude Code, OpenCode, and Cursor.
---

# ChemReact

Execute a deterministic self-contained pipeline:
1. Validate and profile target molecule.
2. Build a visualization plan prompt for post-route-design rendering.
3. Render route assets (target, step, grid, tree).
4. Generate a markdown report and run summary JSON.

## Required Inputs

- `target_smiles`: target molecule SMILES.
- `routes_file`: audited routes JSON array.

## Optional Inputs

- `strategy_file`: strategy analysis JSON.
- `vis_plan_file`: visualization plan JSON produced by LLM.
- `output_dir`: output directory.
- `top_k`: fallback route count when no visualization plan is provided.
- `emit_vis_plan_template`: output path for visualization plan template JSON.

## Execute

Run:

```bash
python skills/ChemReact/scripts/run_closed_loop.py \
  --target-smiles "CCO" \
  --routes-file path/to/audited_routes.json \
  --strategy-file path/to/strategy.json \
  --output-dir out
```

If `vis_plan_file` is omitted, use score-based fallback selection.

Generate only contracts/template:

```bash
python skills/ChemReact/scripts/run_closed_loop.py \
  --output-dir out \
  --emit-vis-plan-template-only
```

## Outputs

- `run_summary.json`: machine-readable pipeline result.
- `visualization_prompt.txt`: prompt for Visualization Specialist persona.
- `vis_plan.template.json`: strict template for host LLM to fill.
- `schemas/*.schema.json`: machine-checkable input contracts.
- `RETRO_REPORT.md`: final report.
- `images/*.png`: generated route images.

## Host Integration

- Read `references/host-integration.md` to wire this skill into Claude Code/OpenCode/Cursor workflows.
- Read `references/data-contract.md` for strict input JSON contracts.

