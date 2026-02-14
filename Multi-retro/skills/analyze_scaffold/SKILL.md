---
name: analyze_scaffold
description: Analyze molecular scaffold for strategic disconnections, ring systems, and synthesis strategy feasibility.
---

# Skill: analyze_scaffold

## Purpose

Analyze the **core scaffold** structure of a molecule to identify:
- Ring systems and fusion patterns
- Strategic disconnection points (BRICS bonds)
- Convergent vs linear synthesis feasibility

**CRITICAL**: This skill analyzes the molecular **backbone/skeleton**, not just functional groups. Use this BEFORE functional group analysis for complex molecules.

## When to Use

- Before starting retrosynthesis on complex polycyclic molecules
- When deciding between linear vs convergent synthesis strategies
- To identify ring-forming reactions (Diels-Alder, RCM, etc.)

## Input Parameters

| Parameter | Type | Required | Description |
|---|---|---|---|
| `smiles` | `string` | Yes | SMILES of the molecule to analyze |

## Output

```json
{
  "success": true,
  "data": {
    "scaffold_type": "fused_polycyclic",
    "murcko_scaffold": "c1ccc2c(c1)Cc1ccccc1C2",
    "ring_count": 3,
    "fused_ring_count": 2,
    "ring_systems": [
      {"ring_id": 0, "size": 6, "is_aromatic": true, "fused_with": [1]}
    ],
    "strategic_bonds": [
      {"atom1_idx": 5, "atom2_idx": 6, "priority": "HIGH", "reason": "BRICS:L1-L3"}
    ],
    "convergent_feasibility": 0.75,
    "linear_feasibility": 0.35,
    "complexity_score": 6.5
  }
}
```

## Key Fields

| Field | Meaning |
|-------|---------|
| `scaffold_type` | `acyclic`, `simple_monocyclic`, `bicyclic_fused`, `fused_polycyclic`, `spiro`, `bridged_polycyclic`, `macrocyclic` |
| `strategic_bonds` | Bonds to disconnect for simplifying the scaffold |
| `convergent_feasibility` | 0-1 score, higher = better for convergent synthesis |
| `linear_feasibility` | 0-1 score, higher = better for linear synthesis |

## Underlying Implementation

- Module: `tools/guidance/strategy_advisor.py` — `analyze_scaffold()`

## CLI Example

```bash
python scripts/run_retro.py skill analyze_scaffold --args '{"smiles": "c1ccc2c(c1)C3CCCc4ccccc4C3C2"}'
```

## Workflow Position

```
analyze_scaffold → plan_ring_strategy → analyze_selectivity → decide_strategy
       ↓                    ↓                    ↓                   ↓
   Scaffold info      Ring plan          FG conflicts        Final decision
```
