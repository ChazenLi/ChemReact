---
name: plan_ring_strategy
description: Plan ring formation strategy including construction order and retro-ring-opening disconnections.
---

# Skill: plan_ring_strategy

## Purpose

Determine the **ring construction strategy** for synthesis:
- Which rings to build first (aromatic → non-aromatic)
- Ring-forming reactions to use (Diels-Alder, RCM, aldol cyclization)
- Retro-ring-opening disconnection points

## When to Use

- After `analyze_scaffold` for polycyclic molecules
- When target has 2+ rings that need strategic construction order
- To identify key ring-forming steps

## Input Parameters

| Parameter | Type | Required | Description |
|---|---|---|---|
| `smiles` | `string` | Yes | SMILES of the molecule to analyze |

## Output

```json
{
  "success": true,
  "data": {
    "scaffold_analysis": {...},
    "ring_strategy": {
      "total_rings": 3,
      "construction_order": [
        {"ring_id": 0, "construction_method": "Diels-Alder / aromatization", "step_order": 1, "is_aromatic": true}
      ],
      "ring_opening_disconnections": [
        {"ring_id": 1, "opening_reaction": "retro-Diels-Alder", "priority": "HIGH", "fragments_description": "diene + dienophile"}
      ],
      "ring_forming_reactions_used": ["Diels-Alder / aromatization", "Nazarov / aldol cyclization"],
      "early_stage_rings": [0, 1],
      "late_stage_rings": [2],
      "timing_rationale": "芳香环优先构建 (稳定性); 复杂环/大环后期构建 (避免副反应)"
    }
  }
}
```

## Ring-Forming Reactions Recognized

| Reaction | Ring Size | Use Case |
|----------|-----------|----------|
| Diels-Alder | 6 | Carbocycle, diene + dienophile |
| RCM | 5-8 | Macrocycles, dienes |
| Aldol cyclization | 5-6 | Enolate chemistry |
| Robinson annulation | 6 | Michael + aldol |
| Nazarov | 5 | Penta-2,4-dienone |
| Pictet-Spengler | 6 | Tetrahydroisoquinoline |
| Fischer Indole | 5 | Indoles from hydrazines |

## Underlying Implementation

- Module: `tools/guidance/strategy_advisor.py` — `plan_ring_strategy()`

## CLI Example

```bash
python scripts/run_retro.py skill plan_ring_strategy --args '{"smiles": "C1CCC2=C(C1)C=CC(=O)C2"}'
```

## Decision Logic

1. **Aromatic rings first** - More stable, set foundation
2. **Small rings → Large rings** - Easier to form first
3. **Fused rings** - Consider shared bond construction order
4. **Macrocycles last** - RCM at late stage
