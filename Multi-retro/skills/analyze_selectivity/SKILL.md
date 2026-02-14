---
name: analyze_selectivity
description: Analyze chemoselectivity issues including competing functional groups, reactivity conflicts, and protection requirements.
---

# Skill: analyze_selectivity

## Purpose

Identify **chemoselectivity issues** that affect synthesis planning:
- Multiple similar FGs competing for same reaction
- Reactivity conflicts between different FGs
- Protection requirements and timing

## When to Use

- When target has multiple similar functional groups (e.g., 3 phenols)
- Before planning reactions that could have selectivity issues
- To determine protection strategy

## Input Parameters

| Parameter | Type | Required | Description |
|---|---|---|---|
| `smiles` | `string` | Yes | SMILES of the molecule to analyze |

## Output

```json
{
  "success": true,
  "data": {
    "functional_groups": [
      {"fg_type": "phenol", "count": 3, "atom_indices": [5, 12, 18]}
    ],
    "competing_groups": [
      {
        "fg_type": "phenol",
        "count": 3,
        "reaction_context": "nucleophilicity, acidity",
        "differentiation_strategy": "steric/electronic",
        "recommendation": "需要正交保护策略，或分步官能团化"
      }
    ],
    "protection_requirements": [
      {
        "fg_type": "alcohol",
        "reason": "selectivity",
        "recommended_pg": "TBS",
        "install_timing": "early",
        "remove_timing": "late"
      }
    ],
    "reactivity_conflicts": [
      {
        "fg1_type": "aldehyde",
        "fg2_type": "primary_amine",
        "conflict_reaction": "reductive_amination",
        "severity": "high",
        "resolution": "醛与胺会直接反应形成亚胺，需控制条件或保护"
      }
    ],
    "overall_complexity": "high"
  }
}
```

## Key Selectivity Issues

| Issue | Example | Resolution |
|-------|---------|------------|
| Multiple phenols | 3 -OH on different rings | Selective protection or use steric difference |
| Aldehyde + Amine | -CHO and -NH2 | Protect aldehyde or control pH |
| Alkene + Aldehyde | Michael acceptor | Selectively reduce one first |

## Functional Groups Detected

- Alcohols, Phenols
- Primary/Secondary Amines
- Aldehydes, Ketones
- Carboxylic Acids, Esters, Amides
- Alkenes, Alkynes
- Aryl/Alkyl Halides

## Underlying Implementation

- Module: `tools/guidance/selectivity_analyzer.py` — `analyze_selectivity()`

## CLI Example

```bash
python scripts/run_retro.py skill analyze_selectivity --args '{"smiles": "Oc1ccc(O)c(O)c1C=O"}'
```

## Protection Group Orthogonality

| PG | Orthogonal To |
|----|--------------|
| Boc | Fmoc, Cbz, TBS |
| Fmoc | Boc, Cbz, TBS |
| TBS | Boc, Fmoc, Ac |
| Bn | TBS, TIPS, Boc |
