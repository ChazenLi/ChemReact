---
name: check_precursor_availability
description: Check commercial availability or synthetic complexity of precursor molecules.
---

# Skill: check_precursor_availability

## Purpose

Determine if a precursor molecule is:
1. **Purchasable** (SA Score < 2.2)
2. **Easily synthesizable** (2.2 ≤ SA Score < 3.5)
3. **Complex** (SA Score ≥ 3.5, requires further decomposition)

This skill helps decide when to **terminate** a branch of the retrosynthesis tree.

## When to Use

- After `break_bond` or `propose_disconnection` generates precursor SMILES
- To determine if a node should be marked `TERMINAL` in RetroGraph
- During recursion decision to evaluate precursor complexity

## Input Parameters

| Parameter | Type | Required | Description |
|---|---|---|---|
| `smiles_list` | `list[string]` or `string` | Yes | List of precursor SMILES (JSON array or comma-separated) |

## Output

```json
{
  "success": true,
  "data": {
    "results": [
      {
        "smiles": "Oc1ccccc1C(=O)O",
        "sa_score": 1.45,
        "classification": "purchasable",
        "is_terminal": true
      },
      {
        "smiles": "CC(=O)Cl",
        "sa_score": 1.12,
        "classification": "purchasable",
        "is_terminal": true
      }
    ]
  }
}
```

## Classification Thresholds (SAThresholds)

| Classification | SA Score Range | Action |
|----------------|----------------|--------|
| `purchasable` | SA < 2.2 | Mark as terminal (可购买) |
| `easily_synthesizable` | 2.2 ≤ SA < 3.5 | Optional decomposition (易合成) |
| `complex` | SA ≥ 3.5 | Continue decomposition (需继续分解) |

## Underlying Implementation

- Module: `tools/chem/sa_scorer.py` — `classify_sa()`
- Constants: `tools/common/constants.py` — `SAThresholds`

## CLI Example

```bash
python scripts/run_retro.py skill check_availability --args '{"smiles_list": ["Oc1ccccc1C(=O)O", "CC(=O)Cl"]}'
```

## Integration with Orchestrator

The Orchestrator uses this skill during `recursion_decision`:
1. Validate reaction passes
2. Call `check_availability` on precursors
3. Use SA Score classification to recommend terminal/decompose
