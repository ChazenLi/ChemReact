# Data Contract

Runtime contract enforcement:
- `run_closed_loop.py` validates `routes_file`, `strategy_file`, and `vis_plan_file` before rendering.
- Validation errors are returned in `run_summary.json` with `success=false`.
- Schema files are emitted to `schemas/*.schema.json` on every run.
- Host adapter request contract: `schemas/host_request.schema.json`.

## `routes_file` JSON

Top-level type: array of route objects.

Minimum route fields:

```json
[
  {
    "route_id": 1,
    "score": 8.7,
    "audit_verdict": "PASS",
    "critical_issues": [],
    "precursors": ["CCO", "O=C=O"],
    "steps": [
      {
        "reaction_class": "Example",
        "reagents": ["NaOH"],
        "conditions": "RT, 1h",
        "reaction_smiles": "CCCl>>CCO"
      }
    ]
  }
]
```

Notes:
- `precursors` and `steps[*].reaction_smiles` are optional but recommended for image generation.
- `reaction_smiles` may also be named `rxn_smiles` or `smirks`.

## `strategy_file` JSON

Top-level object:

```json
{
  "analysis": {
    "core_skeleton": "Biaryl",
    "complexity_features": ["feature1", "feature2"],
    "strategy_type": "Convergent"
  }
}
```

## `vis_plan_file` JSON

Use the schema required by `get_visualization_specialist_prompt` in:
- `retroskill/prompts_personas.py`

Core fields:
- `selected_route_ids`
- `target_image.highlight_atoms`
- `routes[*].route_id`
- `routes[*].overview_grid.enabled`
- `routes[*].tree_view.enabled`
- `routes[*].step_views[*].step_index`

Template generation:

```bash
python ChemReact/scripts/run_closed_loop.py \
  --output-dir out \
  --emit-vis-plan-template-only
```
