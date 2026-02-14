---
name: multi_turn_render
description: Render complete synthesis route from multi-turn analysis session with standardized visualization.
---

# Skill: multi_turn_render

## Purpose

Generate **complete visualization** after multi-turn LLM analysis:
- Reads session history and route graph
- Generates all molecule images
- Creates synthesis tree image
- Produces SYNTHESIS_REPORT.md

## When to Use

- After multiple analysis rounds with LLM
- When route_graph has been populated in session
- To generate final visual report from a session

## Input Parameters

| Parameter | Type | Required | Description |
|---|---|---|---|
| `session_id` | `string` | Yes | Session ID from session_manager |
| `base_dir` | `string` | No | Base output directory (default: `outputs`) |

## Output

```json
{
  "success": true,
  "data": {
    "session_id": "session_20260209_140000_abc12345",
    "session_dir": "outputs/session_20260209_140000_abc12345",
    "report_path": "outputs/session_20260209_140000_abc12345/SYNTHESIS_REPORT.md",
    "images_generated": 12,
    "total_rounds": 5
  }
}
```

## Session Directory Structure

```
outputs/{session_id}/
├── session_meta.json      # Session metadata and history
├── images/
│   ├── target.png         # Target molecule
│   ├── step_001_product.png
│   ├── step_001_reaction.png
│   ├── step_002_product.png
│   ├── step_002_reaction.png
│   └── synthesis_tree.png  # Complete tree
├── route_graph.json       # Route data
└── SYNTHESIS_REPORT.md    # Final report
```

## Underlying Implementation

- Module: `tools/output/report_generator.py`
- Module: `tools/output/visualizer.py`

## CLI Example

```bash
python scripts/run_retro.py skill multi_turn_render --args '{"session_id": "session_20260209_140000_abc12345"}'
```

## Prerequisites

Before calling this skill:
1. Session must exist (created via session_manager)
2. `route_graph` must be populated in session
3. Analysis rounds should be recorded

## Error Handling

```json
{
  "success": false,
  "error": "No route graph in session. Run route analysis first."
}
```

## Image Naming Convention

| Image | Naming |
|-------|--------|
| Target | `target.png` |
| Step N product | `step_NNN_product.png` |
| Step N reaction | `step_NNN_reaction.png` |
| Synthesis tree | `synthesis_tree.png` |
