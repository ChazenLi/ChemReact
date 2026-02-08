# ChemReact Usage

## Scope
- Host orchestration skill for `OpenCode`, `ClaudeCode`, `Cursor`.
- Use host natural-language request to trigger retrosynthesis.
- No model/provider configuration guidance in this file.
- Read in order before execution: `README.md` -> `USAGE.md` -> `workflow.md`.

## Quick Start
```bash
python skills/ChemReact/scripts/run_strict_workflow.py --target-smiles "<TARGET_SMILES>" --output-dir "<OUTPUT_DIR>"
```

## Host Planner Setup
Set a real configured host command before running:
```bash
set CHEMREACT_HOST_PLANNER_COMMAND=<HOST_REAL_PLANNER_COMMAND_WITH_PLACEHOLDERS>
```

Required placeholders:
- `{prompt_file}`
- `{output_file}`
- `{request_file}`
- `{stage}`
- `{round}`

Strict mode rejects:
- placeholder commands
- fallback/mock commands

## Mandatory Runtime Behavior
1. First run chemistry analysis tools on the target molecule.
2. Build planner context from analysis output.
3. Execute configured `host_planner_command` for planning/repair stages.
4. Continue strict loop: plan -> audit -> repair -> re-audit.
5. If repair max iterations is reached, still output report and visualizations for all audited routes.

## Host Natural Language Template
```text
Please run a strict retrosynthesis closed-loop analysis for this molecule.

Target smiles: <TARGET_SMILES>
Output directory: <OUTPUT_DIR>

Requirements:
1) Analyze target chemistry first using tools and build structured planner context.
2) Execute configured host_planner_command for planner and repair stages.
3) Use strict loop (plan -> audit -> repair -> re-audit) with no planner-option exploration.
4) Always produce outputs. If repair max iterations is reached, still output report and route visualizations for PASS/REJECT routes.
5) Report planner_loop and audit key fields from run_summary.json.
```

## Mode Note: --planner-request-only
`--planner-request-only` only emits planning artifacts and exits:
- `planner_request.json`
- `planner_routes.template.json`
- `planner_prompt.txt`
It does not run audit/repair/report/visualization.

## Troubleshooting
- `host planner command is placeholder/fallback`: replace with real configured host command.
- `Missing host planner command`: set `CHEMREACT_HOST_PLANNER_COMMAND` or pass `--host-planner-command`.
- Zero recommended routes with no repair-limit reached: inspect `run_summary.json.audit` and round-level audit details.
- All routes rejected at repair limit: outputs should still be present by design (report + route images).
