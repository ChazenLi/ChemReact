---
name: validate_reaction
description: Validate reaction for atom balance and chemical plausibility. Hard gate - all reactions must pass before adding to RetroGraph.
---

# Skill: validate_reaction

## Purpose

**硬门控** — 验证反应的化学正确性。所有 `break_bond` 或 `propose_disconnection` 生成的反应都必须通过此验证。

检查项：
1. **原子守恒** — 反应物和产物的原子数是否平衡
2. **副产物推断** — 自动推断可能的副产物 (H2O, HCl 等 18 种模式)
3. **选择性提示** — 基于反应类型的注意事项

> **CRITICAL**: 此 Skill 只验证，不修复。验证失败必须调用 `repair_reaction` 或换键。

## When to Use

- `break_bond` 后**必须调用**
- `propose_disconnection` 返回的每个反应**必须验证**
- `repair_reaction` 修复后**必须重新验证**

## Input Parameters

| Parameter | Type | Required | Description |
|---|---|---|---|
| `reaction_smiles` | `string` | Yes | 反应 SMILES，格式 `前体.前体>>产物` |
| `reaction_class` | `string` | No | 预期反应类型 (用于选择性提示) |

## Output

### 成功案例 (PASS)

```json
{
  "success": true,
  "data": {
    "reaction_smiles": "Oc1ccccc1C(=O)O.CC(=O)Cl>>CC(=O)Oc1ccccc1C(=O)O",
    "is_valid": true,
    "verdict": "PASS",
    "issues": [],
    "byproducts": []
  }
}
```

### 需要副产物平衡 (CONDITIONAL)

```json
{
  "success": true,
  "data": {
    "is_valid": false,
    "verdict": "CONDITIONAL",
    "issues": [
      {"code": "ATOM_CONSERVATION_FAIL", "message": "Atom imbalance: {'H': 2, 'O': 1}"}
    ],
    "byproducts": [
      {"formula": "H2O", "name": "Water", "smiles": "O", "confidence": 0.8}
    ]
  }
}
```

**LLM 操作**: 副产物已推断，需要手动添加到反应式中。

### 完全失败 (FAIL)

```json
{
  "success": true,
  "data": {
    "is_valid": false,
    "verdict": "FAIL",
    "issues": [
      {"code": "ATOM_CONSERVATION_FAIL", "message": "Atom imbalance: {'C': -1, 'O': 2}"}
    ],
    "byproducts": []
  }
}
```

**LLM 操作**: 调用 `repair_reaction` 或换键重试。

## Three-Tier Verdict

| Verdict | is_valid | 含义 | LLM 操作 |
|---------|----------|------|----------|
| `PASS` | `true` | 原子平衡，无问题 | 直接进入下一步 |
| `CONDITIONAL` | `false` | 不平衡但可推断副产物 | 根据 `byproducts` 添加副产物后重试 |
| `FAIL` | `false` | 无法解释的不平衡 | 调用 `repair_reaction` 或换键 |

## Underlying Implementation

- Module: `tools/chem/reaction_validator.py` — `validate_reaction()`
- 18 种副产物模式: H2, H2O, HCl, HBr, HI, HF, NH3, CO2, CO, SO2, N2, NaCl, NaBr, KBr, LiCl, CH3OH, C2H5OH, AcOH
- 12 种选择性提示: esterification, amide_coupling, cross_coupling, diels_alder, etc.

## CLI Example

```bash
python scripts/run_retro.py skill validate_reaction --args '{"reaction_smiles": "CC.O>>CCO"}'
```

## LLM Workflow

```
break_bond 返回 reaction_smiles
        ↓
validate_reaction
        ↓
   ┌────┴────┐
   ↓         ↓
 PASS    FAIL/CONDITIONAL
   ↓         ↓
下一步   检查 byproducts
             ↓
        有副产物? ─Yes→ 手动添加后重试
             ↓No
        repair_reaction 或换键
```

## Common Issues and Solutions

| Issue Code | 含义 | 解决方案 |
|------------|------|----------|
| `ATOM_CONSERVATION_FAIL` | 原子不平衡 | 检查 `byproducts`，添加缺失的副产物 |
| `REACTION_FORMAT_INVALID` | 格式错误 | 检查 SMILES 格式，确保有 `>>` |
| `BOND_ANALYSIS_FAILED` | 键分析失败 | 通常是警告，不影响验证 |

## Integration with Orchestrator

```
break_bond → validate_reaction → (PASS) → check_availability → ...
                  ↓
              (FAIL) → repair_reaction → validate_reaction → ...
                              ↓
                          (MAX_RETRIES=3) → 换键重试
```
