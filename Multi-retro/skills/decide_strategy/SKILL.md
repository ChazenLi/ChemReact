---
name: decide_strategy
description: Comprehensive strategy decision (linear vs convergent) based on scaffold, ring, and selectivity analysis.
---

# Skill: decide_strategy

## Purpose

**MASTER SKILL** that runs all 3 analyzers and makes a synthesis strategy decision:
1. Runs `analyze_scaffold` → Scaffold analysis
2. Runs `plan_ring_strategy` → Ring strategy
3. Runs `analyze_selectivity` → Chemoselectivity

Then decides: **Linear**, **Convergent**, or **Mixed** synthesis strategy.

## When to Use

- **First step** for any complex retrosynthesis task
- When you need a comprehensive analysis in one call
- To get LLM-ready strategy recommendation

## Input Parameters

| Parameter | Type | Required | Description |
|---|---|---|---|
| `smiles` | `string` | Yes | Target SMILES to analyze |

## Output

```json
{
  "success": true,
  "data": {
    "recommended_strategy": "convergent",
    "rationale": "多环结构，适合汇聚合成分割片段",
    "scores": {
      "convergent_feasibility": 0.67,
      "linear_feasibility": 0.35
    },
    "scaffold_summary": {
      "type": "fused_polycyclic",
      "ring_count": 5,
      "complexity": 7.2
    },
    "ring_strategy_summary": {
      "total_rings": 5,
      "reactions_used": ["Diels-Alder / aromatization", "Nazarov / aldol cyclization"]
    },
    "selectivity_summary": {
      "competing_groups": 2,
      "overall_complexity": "medium",
      "protection_needed": 1
    },
    "full_analyses": {
      "scaffold": {...},
      "ring_plan": {...},
      "selectivity": {...}
    }
  }
}
```

## Decision Logic

| Condition | Strategy | Rationale |
|-----------|----------|-----------|
| `conv > lin` AND `rings >= 2` | **Convergent** | 多环结构，适合汇聚合成分割片段 |
| `complexity == high` OR `competing >= 2` | **Mixed** | 官能团复杂，需要混合策略 + 选择性保护 |
| Otherwise | **Linear** | 结构相对简单，线性合成效率高 |

## Strategy Types

| Type | Description |
|------|-------------|
| **Linear** | Step-by-step, A→B→C→D. Best for simple molecules. |
| **Convergent** | Build fragments separately, then join. Best for complex polycyclics. |
| **Mixed** | Combination. Use when selectivity issues require staged protection. |

## Underlying Implementation

- Module: `tools/guidance/strategy_advisor.py` — `decide_strategy()`

## CLI Example

```bash
python scripts/run_retro.py skill decide_strategy --args '{"smiles": "OC12C(C(C3=CC=CC=C3)C(C(OC)=O)C2O)(C4=CC=C(OC)C=C4)CC5=CC(OCC6=CC=CC=C6)=CC(OC)=C51"}'
```

## Workflow Position

```
                decide_strategy (ONE CALL)
                       |
            ┌──────────┼──────────┐
            ▼          ▼          ▼
    analyze_scaffold  plan_ring  analyze_selectivity
            |          |          |
            └──────────┼──────────┘
                       ▼
              Strategy Decision
                       ▼
            propose_disconnection
```

## IMPORTANT FOR LLM

- **Use this skill FIRST** for any complex molecule
- The `full_analyses` field contains all raw data for follow-up reasoning
- `recommended_strategy` is a suggestion, LLM should verify with chemical reasoning
