---
name: break_bond
description: Break a specific bond by atom indices to generate precursor SMILES with retrosynthetic guidance.
---

# Skill: break_bond

## Purpose

**精确断键** — 按指定的原子索引断开分子中的键，生成前体 SMILES。

与 `propose_disconnection` 不同，此 Skill 需要 LLM 明确指定断键位置（来自 `analyze_molecule` 的 `strategic_bonds`）。

返回：
- 前体 SMILES 列表
- 反应 SMILES（可直接用于 `validate_reaction`）
- `retro_analysis_guide`（检测到的反应模式和修复建议）

## When to Use

- `analyze_molecule` 返回了 `strategic_bonds` 后
- 需要精确控制断键位置时（而非自动提议）
- `propose_disconnection` 没有找到合适的提议时

## Input Parameters

| Parameter | Type | Required | Description |
|---|---|---|---|
| `smiles` | `string` | Yes | 目标分子 SMILES |
| `atom1_idx` | `int` | Yes | 断裂键的第一个原子索引 (0-based) |
| `atom2_idx` | `int` | Yes | 断裂键的第二个原子索引 (0-based) |

## Output

```json
{
  "success": true,
  "data": {
    "target_smiles": "CC(=O)Oc1ccccc1C(=O)O",
    "bond_broken": {
      "atom1_idx": 6,
      "atom1_symbol": "C",
      "atom2_idx": 7,
      "atom2_symbol": "O",
      "bond_type": "SINGLE",
      "in_ring": false
    },
    "fragments": ["Oc1ccccc1C(=O)O", "CC(=O)O"],
    "num_fragments": 2,
    "reaction_smiles": "Oc1ccccc1C(=O)O.CC(=O)O>>CC(=O)Oc1ccccc1C(=O)O",
    "retro_analysis_guide": {
      "bond_type": "C(=O)-O (ester bond)",
      "retro_pattern": "ester hydrolysis",
      "broken_atoms": {
        "atom1": {"idx": 6, "symbol": "C", "in_ring": false, "is_aromatic": false},
        "atom2": {"idx": 7, "symbol": "O", "in_ring": false, "is_aromatic": false}
      },
      "fragment_profiles": [
        {"index": 0, "smiles": "Oc1ccccc1C(=O)O", "valid": true, "functional_groups": ["OH", "C=O"]}
      ],
      "retro_hint": "Amide retrosynthesis: after C(=O)-N cleavage, the C=O end should become an acyl chloride...",
      "llm_task": "Detected ester hydrolysis (C(=O)-O). Use fragment_profiles to determine nucleophile/electrophile roles..."
    },
    "next_step": "Use validate_reaction to verify: --reaction_smiles \"Oc1ccccc1C(=O)O.CC(=O)O>>CC(=O)Oc1ccccc1C(=O)O\""
  }
}
```

## Detected Retro Patterns

当检测到已知反应模式时，返回 `retro_analysis_guide`:

| Pattern | 键类型 | 前体修复提示 |
|---------|--------|-------------|
| `amide_coupling` | C(=O)-N 酰胺键 | C 端 → 酰氯 C(=O)Cl；N 端 → 胺 |
| `sulfonamide_formation` | S(=O)2-N | S 端 → 磺酰氯 -SO2Cl；N 端 → 胺 |
| `ester_hydrolysis` | C(=O)-O 酯键 | C 端 → 酰氯或酸；O 端 → 醇 |
| `Heck_reaction` | Ar-C (aryl-vinyl) | 芳基端 → Ar-Br/I |
| `N_alkylation` | N-C (非酰胺) | N 端 → N-H；C 端 → 卤代物 |
| `diarylamine_formation` | Ar-N-Ar | 芳基端 → Ar-X；胺端 → Ar-NH2 |

## Underlying Implementation

- Module: `tools/skills/break_bond.py` — `BreakBondSkill`
- Uses: `tools/chem/mol_parser.py` — SMILES 解析和标准化
- Uses: `tools/guidance/retro_guide.py` — 逆合成指导构建

## CLI Example

```bash
python scripts/run_retro.py skill break_bond --args '{"smiles": "CC(=O)Oc1ccccc1C(=O)O", "atom1_idx": 6, "atom2_idx": 7}'
```

## LLM Workflow

```
analyze_molecule → 读取 atom_bond_map.bonds 或 retro_guidance
        ↓
找到 strategic_priority=1 的键: atom1_idx=X, atom2_idx=Y
        ↓
break_bond --args '{"smiles": "...", "atom1_idx": X, "atom2_idx": Y}'
        ↓
检查 retro_analysis_guide (如有)
        ↓
validate_reaction --args '{"reaction_smiles": "..."}'
        ↓
PASS → check_availability
FAIL → repair_reaction 或换键
```

## Error Handling

### 原子索引越界

```json
{
  "success": false,
  "error": "Atom index out of range. Molecule has 13 atoms (0-12). Got atom1_idx=15"
}
```

### 键不存在

```json
{
  "success": false,
  "error": "No bond between atom 3 (C) and atom 7 (O)"
}
```

### 芳香键不可断

```json
{
  "success": false,
  "error": "Bond 5-6 is aromatic. Breaking aromatic bonds is chemically unreasonable."
}
```

### 断键未产生碎片

```json
{
  "success": false,
  "error": "Bond break did not produce valid fragments. This may happen with ring bonds that don't split the molecule."
}
```

## Important Notes

1. **原子索引基于 canonical_smiles** — 使用 `analyze_molecule` 返回的 `canonical_smiles` 对应的索引
2. **芳香键不可断** — 系统会拒绝断裂芳香键
3. **必须验证** — 每次断键后都必须调用 `validate_reaction`
4. **检查 retro_analysis_guide** — 可能需要根据提示修正前体 SMILES

## Integration with Orchestrator

```
analyze_molecule
       ↓
   strategic_bonds 列表
       ↓
break_bond (选择 priority=1 的键)
       ↓
validate_reaction
       ↓
  ┌────┴────┐
PASS       FAIL
  ↓          ↓
check_    repair_reaction
availability    ↓
             validate_reaction
                   ↓
              (MAX_REPAIR_RETRIES=3)
                   ↓
              换键重试
```
