# ChemReact README

ChemReact is a host-oriented retrosynthesis skill with a strict closed loop:
`analyze molecule -> plan routes -> strict audit -> repair loop -> report/visualization`.

## What This Skill Is
- A workflow skill for host orchestrators (`OpenCode`, `ClaudeCode`, `Cursor`).
- Not a standalone built-in planner engine.
- Uses host execution wiring through `host_planner_command`.

## Non-Negotiable Rules
1. Analysis-first is mandatory.
2. Planner generation must use structured chemistry context from analysis.
3. Strict audit is always enforced in operational runs.
4. If repair iterations hit max limit, output must still be produced (report + visualizations for PASS/REJECT routes).
5. No planner-option exploration, no provider menu chatter, no fallback/mock command in strict workflow.

## Required Reading Order (for LLM and Integrators)
1. `SKILL.md`
2. `USAGE.md`
3. `workflow.md`

## Main Entry
```bash
python skills/ChemReact/scripts/run_strict_workflow.py --target-smiles "<TARGET_SMILES>" --output-dir "<OUTPUT_DIR>"
```

## Core Outputs
- `<OUTPUT_DIR>/run_summary.json`
- `<OUTPUT_DIR>/RETRO_REPORT.md`
- `<OUTPUT_DIR>/images/route_*_overview.png`
- `<OUTPUT_DIR>/images/route_*_tree.png`
- Planner loop artifacts: `planner_request.json`, `planner_prompt.txt`, audited/repair round files.

## Pre-Upload Checklist
- All docs in this folder are consistent with code behavior.
- `host_planner_command` contract is documented and enforced.
- Analysis-first requirement is visible in docs and planner request contract.
- Repair-limit behavior is documented: still output report/visualization.
- Tests pass:
  - `pytest -q skills/ChemReact/tests`
