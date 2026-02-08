# ChemReact Skill

ChemReact is a strict retrosynthesis closed-loop skill for host software (`OpenCode` / `ClaudeCode` / `Cursor`).

## Positioning
- Input: target molecule (`target_smiles`).
- Execution: host runs a real configured `host_planner_command`.
- Process: tool-based chemistry analysis -> planner generation -> strict audit -> repair loop -> report/visualization.
- Output: traceable structured artifacts for both machine and human review.

## Required Execution Order
1. Analyze target chemistry first (structure/features/functional groups/ring topology).
2. Build planner request with structured chemistry context.
3. Run planner invocation via `host_planner_command`.
4. Run strict audit and repair loop (up to max iterations).
5. Always output report and visualizations; if max repair limit is reached, include both PASS and REJECT routes.

## Input Contract
- `target_smiles`
- `output_dir`
- `host_planner_command` (or `CHEMREACT_HOST_PLANNER_COMMAND`)

## Host Command Contract
`host_planner_command` must:
- be real (no placeholder/fallback/mock)
- include placeholders:
  - `{prompt_file}`
  - `{output_file}`
  - `{request_file}`
  - `{stage}`
  - `{round}`
- be executable by host integration as one complete command string

## Forbidden Assistant Output
Do not output:
- "Let me check available planner options"
- provider/API option menus
- API-key setup detours
- fallback experiments in strict mode

Required behavior:
- Execute the configured command path and continue strict loop execution.

## Single Command Entry
```bash
python skills/ChemReact/scripts/run_strict_workflow.py --target-smiles "<TARGET_SMILES>" --output-dir "<OUTPUT_DIR>"
```

## Completion Criteria
Required artifacts:
- `<OUTPUT_DIR>/RETRO_REPORT.md`
- `<OUTPUT_DIR>/run_summary.json`
- `<OUTPUT_DIR>/images/route_*_overview.png`
- `<OUTPUT_DIR>/images/route_*_tree.png`

## Key Summary Fields
- `planner_loop.rounds`
- `planner_loop.max_iterations`
- `planner_loop.termination_reason`
- `planner_loop.reached_repair_limit`
- `audit.recommended_routes`
- `audit.rejected_routes`
