---
name: render_synthesis_report
description: Generate a visual report (Markdown + PNG) from a completed RetroGraph.
---

# Skill: render_synthesis_report

## Purpose

After the retrosynthesis plan is complete, render:
1. A **Markdown report** summarizing the synthesis routes
2. **PNG images** of the synthesis tree and individual reactions

This is the **final output** of a successful retrosynthesis session.

## When to Use

- When RetroGraph is complete (all nodes resolved)
- When user requests visualization of current progress
- At the end of `finalize` command

## Input Parameters

| Parameter | Type | Required | Description |
|---|---|---|---|
| `retro_graph` | `object` | Yes | Serialized RetroGraph (dict or RetroGraph instance) |
| `output_dir` | `string` | No | Output directory (default: `"outputs/report"`) |

## Output

```json
{
  "success": true,
  "data": {
    "output_dir": "outputs/report",
    "report_content": "# Synthesis Report\n\n## Target: ..."
  }
}
```

## Underlying Implementation

- Module: `tools/output/report_generator.py` — `generate_synthesis_report()`
- Data Model: `tools/models/output_models.py` — `RetroGraph`

## CLI Example

```bash
python scripts/run_retro.py skill render_report --args '{"retro_graph": {"target_smiles": "...", "nodes": {...}, "edges": {...}}, "output_dir": "outputs/my_report"}'
```

## Report Contents

The generated Markdown report includes:
- Target molecule info (SMILES, MW, SA Score)
- Recommended synthesis route(s)
- Step-by-step reaction details
- Precursor availability summary
- Alternative routes (if any)

## Output Directory Structure

```
outputs/report/
├── SYNTHESIS_REPORT.md    # Markdown report
├── retro_graph.json       # Graph data
└── images/
    ├── target.png
    ├── step_001.png
    └── synthesis_tree.png
```

## Integration with Orchestrator

Called by `finalize` command:
```
finalize → render_report → SYNTHESIS_REPORT.md
```
