# ChemReact Workflow

## Main Flow (Host-Native)
1. Host natural-language request provides target molecule.
2. Pre-planning chemistry analysis runs first (structure/features/functional groups/ring-topology context).
3. Planner request is built with structured `target_chemistry`.
4. Configured `host_planner_command` is executed to generate route candidates.
5. Strict audit validates chemistry and route consistency.
6. Failed routes enter repair loop and are re-audited (max 5 rounds).
7. Outputs are generated:
   - Normal case: report + visualizations for recommended routes (strict report still includes audited detail).
   - Repair-limit case: report + visualizations still generated for all audited routes (PASS and REJECT).

## Strict Audit Coverage (Summary)
- schema validity and route structure checks
- parsable reaction steps (`reaction_smiles` / `rxn_smiles` / `smirks`)
- atom conservation and byproduct plausibility
- reaction class vs bond-change consistency
- ring/heterocycle topology continuity
- final-step product matching target
- portfolio style constraints in auto-propose mode

## Intermediate (Non-Final) Modes
- `planner_request_only`
- `planner_generation_pending`
- `planner_repair_requested`
- `planner_repair_pending`

## Completion Conditions
Final run must include:
- `RETRO_REPORT.md`
- `run_summary.json`
- `images/route_*_overview.png`
- `images/route_*_tree.png`

## Entry Command
```bash
python skills/ChemReact/scripts/run_strict_workflow.py --target-smiles "<TARGET_SMILES>" --output-dir "<OUTPUT_DIR>"
```

## Host Command Contract
`host_planner_command` must:
- be real (no fallback/mock/placeholder)
- include placeholders: `{prompt_file}` `{output_file}` `{request_file}` `{stage}` `{round}`
- be executed directly by host integration without planner-option exploration
