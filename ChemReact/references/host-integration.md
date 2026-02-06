# Host Integration

## Shared Pattern

Use a two-turn closed loop in any host LLM tool:

1. Route design + audit phase
- Ask model to produce strategy and audited routes JSON.

2. Visualization + reporting phase
- Feed strategy/routes JSON to `run_closed_loop.py`.
- Optionally call model with `visualization_prompt.txt` and pass returned JSON as `--vis-plan-file`.

## Claude Code

- Keep outputs under a workspace folder (for example `outputs/run_001`).
- Direct run:

```bash
python skills/ChemReact/scripts/run_closed_loop.py \
  --target-smiles "<SMILES>" \
  --routes-file outputs/run_001/routes.json \
  --strategy-file outputs/run_001/strategy.json \
  --output-dir outputs/run_001
```

- Adapter run:

```bash
python skills/ChemReact/adapters/claudecode_adapter.py \
  --request-file skills/ChemReact/samples/host_request.json
```

## OpenCode

- Store route JSON in project files and execute the same command.
- Use `run_summary.json` as structured state for the next turn.
- Adapter run:

```bash
python skills/ChemReact/adapters/opencode_adapter.py \
  --request-file skills/ChemReact/samples/host_request.json
```

## Cursor

- Add a task/command preset for `run_closed_loop.py`.
- Bind target/routing JSON files as arguments from workspace paths.
- Adapter run:

```bash
python skills/ChemReact/adapters/cursor_adapter.py \
  --request-file skills/ChemReact/samples/host_request.json
```

## Recommended Operational Guards

- Validate JSON before run.
- Keep route IDs stable across iterations.
- Version outputs by timestamped directory names.
- Validate adapter request payload against `schemas/host_request.schema.json`.

