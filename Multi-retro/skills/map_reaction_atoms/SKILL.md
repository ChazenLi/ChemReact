---
name: map_reaction_atoms
description: Add atom-to-atom mapping to a reaction SMILES using RXNMapper.
---

# Skill: map_reaction_atoms

## Purpose

Given an unmapped reaction SMILES, compute atom-atom mapping (AAM) to enable:
- Precise bond change analysis
- Reaction mechanism visualization
- Template extraction

## When to Use

- When `break_bond` or `propose_disconnection` generates a reaction without mapping
- Before detailed validation if mapping is needed
- For reaction visualization and documentation

## Input Parameters

| Parameter | Type | Required | Description |
|---|---|---|---|
| `reaction_smiles` | `string` | Yes | Unmapped reaction SMILES |

## Output

```json
{
  "success": true,
  "data": {
    "input_smiles": "CC.O>>CCO",
    "mapped_smiles": "[CH3:1][CH3:2].[OH2:3]>>[CH3:1][CH2:2][OH:3]",
    "confidence": 0.95
  }
}
```

If mapping fails:

```json
{
  "success": false,
  "error": "Mapping failed: no valid correspondence found."
}
```

## Underlying Implementation

- Module: `tools/chem/atom_mapper.py` — `map_reaction()`
- Uses: RXNMapper (transformer-based model)

## CLI Example

```bash
python scripts/run_retro.py skill map_atoms --args '{"reaction_smiles": "CC(=O)O.Oc1ccccc1>>CC(=O)Oc1ccccc1"}'
```

## Notes

- Atom mapping is computationally intensive; cache results when possible
- The `confidence` score indicates model certainty
- RXNMapper requires rxnmapper package to be installed

## Integration with Orchestrator

Optional step in the workflow:
```
break_bond → (map_atoms) → validate_reaction → ...
```
