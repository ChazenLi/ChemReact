---
name: propose_disconnection
description: Propose strategic retrosynthetic disconnection points for a target molecule.
---

# Skill: propose_disconnection

## Purpose

Given a target molecule, identify **strategic bonds** to disconnect and propose precursor molecules. This is the core "planning" step in retrosynthesis.

## When to Use

- After `analyze_molecule` has characterized the target
- When expanding the retrosynthesis DAG frontier
- When you want automatic bond selection (vs. manual `break_bond`)

## Input Parameters

| Parameter | Type | Required | Description |
|---|---|---|---|
| `target_smiles` | `string` | Yes | SMILES of the molecule to disconnect |
| `strategy_hint` | `string` | No | Optional hint (e.g., "use cross-coupling", "avoid stereocenters") |
| `max_proposals` | `int` | No | Maximum number of proposals (default: 5) |

## Output

```json
{
  "success": true,
  "data": {
    "target_smiles": "CC(=O)Oc1ccccc1C(=O)O",
    "disconnections": [
      {
        "id": "disc_001",
        "strategy": "ester_hydrolysis",
        "bond_description": "C-O ester bond",
        "precursors": ["Oc1ccccc1C(=O)O", "CC(=O)O"],
        "reaction_class": "Fischer Esterification (reverse)",
        "confidence": 0.85,
        "rationale": "Simple ester disconnection to salicylic acid and acetic acid."
      },
      {
        "id": "disc_002",
        "strategy": "aryl_functionalization",
        "bond_description": "C-C aryl bond",
        "precursors": ["c1ccccc1", "CC(=O)Cl"],
        "reaction_class": "Friedel-Crafts Acylation",
        "confidence": 0.60
      }
    ]
  }
}
```

## Underlying Implementation

- Module: `tools/guidance/disconnection_proposer.py` — `propose_disconnections()`
- Constants: `tools/common/constants.py` — `DEFAULT_MAX_PROPOSALS = 5`

## CLI Example

```bash
python scripts/run_retro.py skill propose_disconnection --args '{"target_smiles": "CC(=O)Oc1ccccc1C(=O)O"}'
```

## Important Notes

- This skill **proposes**, it does not **validate**
- All proposed reactions must be passed to `validate_reaction`
- The `confidence` score is heuristic; do not use as a hard threshold
- For precise control, use `break_bond` with specific atom indices

## Integration with Orchestrator

```
analyze_molecule → propose_disconnection → validate_reaction → ...
       ↓                    ↓                     ↓
   atom_bond_map      disconnection list      is_valid check
```
