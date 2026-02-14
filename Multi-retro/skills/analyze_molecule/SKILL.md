---
name: analyze_molecule
description: Analyze target molecule structure, functional groups, SA score, and identify strategic bonds for retrosynthetic disconnection.
---

# Skill: analyze_molecule

## Purpose

**核心分析 Skill** — 分析目标分子的结构特征，为逆合成策略提供决策依据。

返回信息包括：
- **atom_bond_map**: 原子-键映射，包含每个键的 `chemical_label`、`strategic_priority`、`retro_hint`
- **retro_guidance**: 逆合成指导文本，按优先级排序的断键建议
- **functional_groups**: 官能团识别结果
- **sa_score**: 合成可及性评分 (1-10，越低越容易)

## When to Use

- 逆合成任务的**第一步**
- 评估新前体分子时
- 需要识别战略断键位置时

## Input Parameters

| Parameter | Type | Required | Description |
|---|---|---|---|
| `smiles` | `string` | Yes | 待分析分子的 SMILES |

## Output

```json
{
  "success": true,
  "data": {
    "smiles": "CC(=O)Oc1ccccc1C(=O)O",
    "canonical_smiles": "CC(=O)Oc1ccccc1C(=O)O",
    "molecular_weight": 180.16,
    "sa_score": 1.85,
    "atom_bond_map": {
      "canonical_smiles": "CC(=O)Oc1ccccc1C(=O)O",
      "atoms": [
        {"idx": 0, "symbol": "C", "aromatic": false, "in_ring": false, "degree": 1, "num_hs": 3, "neighbors": [1]}
      ],
      "bonds": [
        {
          "idx": 0,
          "atom1_idx": 0,
          "atom2_idx": 1,
          "bond_type": "SINGLE",
          "in_ring": false,
          "breakable": true,
          "chemical_label": "烷基 C-C 键",
          "potential_reactions": ["Grignard加成", "烷基化"],
          "retro_hint": "断裂 -> 格氏试剂 + 羰基 或 两个自由基片段",
          "strategic_priority": 3
        },
        {
          "idx": 5,
          "atom1_idx": 6,
          "atom2_idx": 7,
          "bond_type": "SINGLE",
          "in_ring": false,
          "breakable": true,
          "chemical_label": "酯键 C(=O)-O",
          "potential_reactions": ["酯水解/Fischer酯化", "转酯化"],
          "retro_hint": "断裂 -> 羧酸 + 醇 (经典逆合成断裂)",
          "strategic_priority": 1
        }
      ],
      "ring_info": [
        {"ring_id": 0, "size": 6, "is_aromatic": true, "ring_type": "苯环", "formation_hints": []}
      ],
      "retro_guidance": "【推荐断键位置 (按优先级排序)】\n  1. atom 6-7 (酯键 C(=O)-O), priority=1\n     反应: 酯水解/Fischer酯化, 转酯化\n     提示: 断裂 -> 羧酸 + 醇 (经典逆合成断裂)\n\n【建议操作顺序】\n  Step 1: break_bond atom1_idx=6 atom2_idx=7 (酯键 C(=O)-O)"
    },
    "functional_groups": {
      "ester": 1,
      "carboxyl": 1
    },
    "protecting_groups": {}
  }
}
```

## Key Fields for LLM Decision

### strategic_priority (断键优先级)

| 值 | 含义 | LLM 操作建议 |
|---|---|---|
| `1` | **高优先级** | 优先选择此键断裂 (酰胺、酯、磺酰胺、联芳键等) |
| `2` | **中优先级** | 可选断裂 (醚键、N-烷基、C=C 等) |
| `3` | **低优先级** | 最后考虑 (普通烷基 C-C) |

### chemical_label (键类型)

常见高价值键类型：
- `酯键 C(=O)-O` → 酯水解
- `酰胺键 C(=O)-N` → 酰胺偶联
- `磺酰胺键 SO2-N` → 磺酰胺形成
- `联芳 C-C 键` → Suzuki 偶联
- `二芳胺键 Ar-NH-Ar` → Buchwald-Hartwig
- `α-羰基 C-C 键` → Aldol/Michael

### retro_guidance

系统生成的**操作建议文本**，包含：
1. 按 `strategic_priority` 排序的断键列表
2. 每个键对应的反应类型
3. 建议的操作顺序

**LLM 应该优先参考 `retro_guidance` 中的建议。**

## Underlying Implementation

- Module: `tools/chem/structure_analyzer.py` — `analyze_molecule()`
- Data Model: `tools/models/chem_models.py` — `MoleculeInfo`, `AtomBondMap`, `BondInfo`

## CLI Example

```bash
python scripts/run_retro.py skill analyze_molecule --args '{"smiles": "CC(=O)Oc1ccccc1C(=O)O"}'
```

## LLM Workflow After This Skill

```
analyze_molecule 返回 atom_bond_map 和 retro_guidance
        ↓
LLM 读取 retro_guidance，找到 strategic_priority=1 的键
        ↓
调用 break_bond --args '{"smiles": "...", "atom1_idx": X, "atom2_idx": Y}'
        ↓
调用 validate_reaction 验证反应
        ↓
如果失败，尝试下一个优先级的键
```

## Error Handling

```json
{
  "success": false,
  "error": "Invalid SMILES: RDKit cannot parse"
}
```

## Integration with Orchestrator

Orchestrator 在以下场景调用此 Skill:
1. `plan` 创建初始路线后
2. `break_bond` 或 `propose_disconnection` 之前
3. 递归分解时评估新前体
