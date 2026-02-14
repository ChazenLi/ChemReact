---
name: repair_reaction
description: Attempt to repair a chemically invalid reaction SMILES.
---

# Skill: repair_reaction

## Purpose

When `validate_reaction` returns `is_valid: false`, this skill attempts to **fix** the reaction by:
1. Adding missing byproducts (e.g., H2O, HCl)
2. Correcting reagent stoichiometry
3. Suggesting alternative reaction conditions

> **WARNING**: This skill uses heuristics + LLM guidance. Repaired reactions **must** be re-validated via `validate_reaction`.

## When to Use

- After `validate_reaction` fails
- Before giving up on a proposed disconnection
- Maximum retry: `MAX_REPAIR_RETRIES = 3`

## Input Parameters

| Parameter | Type | Required | Description |
|---|---|---|---|
| `reaction_smiles` | `string` | Yes | The invalid reaction SMILES |
| `issues` | `list` | Yes | Error codes/messages from `validate_reaction` |

## Output

```json
{
  "success": true,
  "data": {
    "reaction_smiles": "CC.O>>CCO",
    "issues": [
      {"code": "ATOM_CONSERVATION_FAIL", "message": "Oxygen not balanced"}
    ],
    "guidance": [
      "Add explicit water molecule to balance oxygen",
      "Consider: CC.[OH2:1]>>[CH3:2][CH2:3][OH:1]"
    ],
    "repair_strategy": "add_byproduct",
    "categories_affected": ["atom_balance"],
    "revalidation": {
      "is_valid": true,
      "verdict": "PASS"
    }
  }
}
```

If repair is not possible:

```json
{
  "success": false,
  "error": "Unable to repair: fundamental disconnection logic is incorrect."
}
```

## Underlying Implementation

- Module: `tools/guidance/repair_advisor.py` — `generate_repair_guidance()`
- Module: `tools/chem/reaction_validator.py` — `validate_reaction()` (for re-validation)
- Constants: `tools/common/constants.py` — `MAX_REPAIR_RETRIES = 3`

## CLI Example

```bash
python scripts/run_retro.py skill repair_reaction --args '{"reaction_smiles": "CC>>CCO", "issues": ["ATOM_CONSERVATION_FAIL"]}'
```

## Integration with Orchestrator

```
Orchestrator -> validate_reaction(rxn)
       ↓
   FAIL, issues=[...]
       ↓
Orchestrator -> repair_reaction(rxn, issues)
       ↓
   repaired_rxn
       ↓
Orchestrator -> validate_reaction(repaired_rxn)
       ↓
   PASS or retry
```

## Retry Limit

After `MAX_REPAIR_RETRIES = 3` failed attempts:
- Use `abandon_bond` to try alternative disconnection
- Use `retry_other_bond` to select different strategic bond
