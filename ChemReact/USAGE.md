# ChemReact Usage

## 1. Generate Template and Schemas Only

Use this when host LLM has not produced a visualization plan yet.

```bash
python ChemReact/scripts/run_closed_loop.py \
  --output-dir ChemReact/sample_output_template \
  --emit-vis-plan-template-only
```

Outputs:
- `vis_plan.template.json`
- `schemas/routes.schema.json`
- `schemas/strategy.schema.json`
- `schemas/vis_plan.schema.json`

## 2. Validate Inputs Only

```bash
python ChemReact/scripts/run_closed_loop.py \
  --target-smiles "CCO" \
  --routes-file ChemReact/samples/routes.json \
  --strategy-file ChemReact/samples/strategy.json \
  --vis-plan-file ChemReact/samples/vis_plan.json \
  --output-dir ChemReact/sample_output_validate \
  --validate-only
```

## 3. Run Full Closed Loop

```bash
python ChemReact/scripts/run_closed_loop.py \
  --target-smiles "CCCCC1=NC(Cl)=C(CO)N1CC2=CC=C(C3=CC=CC=C3C4=NNN=N4)C=C2" \
  --routes-file ChemReact/samples/routes.json \
  --strategy-file ChemReact/samples/strategy.json \
  --vis-plan-file ChemReact/samples/vis_plan.json \
  --output-dir ChemReact/sample_output
```

Main outputs:
- `RETRO_REPORT.md`
- `images/*.png`
- `visualization_prompt.txt`
- `run_summary.json`

## 4. Run via Host Adapter

Use shared request payload:
- `ChemReact/samples/host_request.json`

Claude Code:

```bash
python ChemReact/adapters/claudecode_adapter.py \
  --request-file ChemReact/samples/host_request.json
```

OpenCode:

```bash
python ChemReact/adapters/opencode_adapter.py \
  --request-file ChemReact/samples/host_request.json
```

Cursor:

```bash
python ChemReact/adapters/cursor_adapter.py \
  --request-file ChemReact/samples/host_request.json
```
