---
name: get_global_strategy
description: Get global strategy context for LLM prompts with PG policy and method preferences.
---

# Skill: get_global_strategy

## Purpose

获取**全局策略上下文**，为 LLM 提供合成策略的宏观指导：
- 优先方法偏好 (PREFERRED_METHODS)
- 保护基操作限制
- 最大合成步骤数约束

此 Skill 通常在复杂分子的逆合成分析开始时调用，为后续决策提供策略框架。

## When to Use

- 复杂分子的逆合成分析开始时
- 需要了解全局方法偏好时
- 需要设置保护基策略约束时

## Input Parameters

| Parameter | Type | Required | Description |
|---|---|---|---|
| `target_smiles` | `string` | Yes | 目标分子 SMILES |
| `max_pg_operations` | `int` | No | 最大保护基操作数 (默认: 2) |
| `max_steps` | `int` | No | 最大合成步骤数 (默认: 10) |
| `constraints` | `object` | No | 额外约束 |

## Output

```json
{
  "success": true,
  "data": {
    "target_smiles": "CC(=O)Oc1ccccc1C(=O)O",
    "strategy_context": {
      "preferred_methods": [
        "amide_coupling",
        "cross_coupling",
        "esterification",
        "nucleophilic_substitution"
      ],
      "pg_policy": {
        "max_operations": 2,
        "preferred_pg": ["Boc", "TBS", "Ac"],
        "orthogonal_sets": [["Boc", "Fmoc"], ["TBS", "TIPS"]]
      },
      "step_constraints": {
        "max_steps": 10,
        "preferred_step_range": [3, 7]
      }
    },
    "preferred_methods": {
      "amide_coupling": {"priority": 1, "conditions": "EDC/HOBt, room temp"},
      "cross_coupling": {"priority": 1, "conditions": "Pd-catalyzed, base"},
      "esterification": {"priority": 2, "conditions": "DCC/DMAP or acid chloride"},
      "nucleophilic_substitution": {"priority": 2, "conditions": "SN2, base"}
    },
    "constraints_applied": {
      "max_pg_operations": 2,
      "max_steps": 10
    }
  }
}
```

## Key Fields

### preferred_methods

LLM 应优先考虑的合成方法及其优先级:
- `priority=1`: 高优先级 (酰胺偶联、交叉偶联)
- `priority=2`: 中优先级 (酯化、亲核取代)
- `priority=3`: 低优先级 (氧化还原等)

### pg_policy

保护基策略约束:
- `max_operations`: 最大保护/脱保护操作数
- `preferred_pg`: 推荐的保护基
- `orthogonal_sets`: 正交保护基组合

## Underlying Implementation

- Module: `tools/legacy/planner/global_strategy_manager.py`
- Constants: `PREFERRED_METHODS` 字典

## CLI Example

```bash
python scripts/run_retro.py skill get_global_strategy --args '{"target_smiles": "CC(=O)Oc1ccccc1C(=O)O"}'
```

## Integration with Workflow

```
get_global_strategy (可选)
        ↓
  策略上下文 (方法偏好 + PG 约束)
        ↓
decide_strategy (综合分析)
        ↓
analyze_molecule → break_bond → validate_reaction → ...
```

## Notes

- 此 Skill 是**可选的**，主要用于复杂分子的策略规划
- 返回的 `preferred_methods` 是建议性的，LLM 可以根据具体情况调整
- `pg_policy` 可以指导后续 `analyze_selectivity` 的保护基推荐
