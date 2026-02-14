# skills/ - Skill 文档目录

每个 Skill 有独立文件夹，包含 `SKILL.md` 定义文档。实际实现在 `tools/skills/`。

## Skill 列表 (14 个)

| Skill | 用途 | 文档 |
|-------|------|------|
| `analyze_molecule` | 分子结构分析 + atom_bond_map + strategic_bonds | [SKILL.md](analyze_molecule/SKILL.md) |
| `analyze_scaffold` | 骨架分析、战略断裂点、环系统 | [SKILL.md](analyze_scaffold/SKILL.md) |
| `analyze_selectivity` | 化学选择性、竞争官能团、保护策略 | [SKILL.md](analyze_selectivity/SKILL.md) |
| `break_bond` | 精确断键 (按原子索引) | [SKILL.md](break_bond/SKILL.md) |
| `check_precursor_availability` | 前体可用性检查 (基于 SA Score) | [SKILL.md](check_precursor_availability/SKILL.md) |
| `decide_strategy` | 综合策略决策 (线性/汇聚/混合) | [SKILL.md](decide_strategy/SKILL.md) |
| `get_global_strategy` | 全局策略上下文 (PG 策略 + 方法偏好) | [SKILL.md](get_global_strategy/SKILL.md) |
| `map_reaction_atoms` | 原子映射 (RXNMapper) | [SKILL.md](map_reaction_atoms/SKILL.md) |
| `multi_turn_render` | 多轮会话渲染 | [SKILL.md](multi_turn_render/SKILL.md) |
| `plan_ring_strategy` | 环形成策略规划 | [SKILL.md](plan_ring_strategy/SKILL.md) |
| `propose_disconnection` | 断裂点自动提议 (规则匹配) | [SKILL.md](propose_disconnection/SKILL.md) |
| `render_synthesis_report` | 生成合成报告 | [SKILL.md](render_synthesis_report/SKILL.md) |
| `repair_reaction` | 反应修复 (基于 issues) | [SKILL.md](repair_reaction/SKILL.md) |
| `validate_reaction` | 反应验证 (硬门控) | [SKILL.md](validate_reaction/SKILL.md) |

## 核心工作流 Skill

逆合成分析的核心流程使用以下 5 个 Skill:

```
analyze_molecule → break_bond → validate_reaction → check_availability → (递归)
       ↓               ↓              ↓                   ↓
  strategic_bonds   fragments     is_valid          SA Score
       ↓               ↓              ↓                   ↓
   atom1_idx      reaction_smiles  PASS/FAIL      terminal/decompose
   atom2_idx
```

## 辅助分析 Skill

用于复杂分子的深度分析:

| Skill | 输入 | 输出 | 使用时机 |
|-------|------|------|----------|
| `analyze_scaffold` | smiles | 骨架类型、convergent_feasibility | 多环复杂分子 |
| `analyze_selectivity` | smiles | 竞争官能团、保护需求 | 多官能团分子 |
| `decide_strategy` | smiles | recommended_strategy | 任意目标（综合分析） |
| `plan_ring_strategy` | smiles | 环构建顺序、环开环方案 | 多环分子 |

## 源代码位置

所有 Skill 实现位于 `tools/skills/`:

```
tools/skills/
├── base.py                    # BaseSkill, SkillResult
├── analyze_molecule.py        # AnalyzeMoleculeSkill
├── analyze_scaffold.py        # AnalyzeScaffoldSkill
├── analyze_selectivity.py     # AnalyzeSelectivitySkill
├── break_bond.py              # BreakBondSkill
├── check_availability.py      # CheckAvailabilitySkill
├── decide_strategy.py         # DecideStrategySkill
├── get_global_strategy.py     # GetGlobalStrategySkill
├── map_atoms.py               # MapAtomsSkill
├── plan_ring_strategy.py      # PlanRingStrategySkill
├── propose_disconnection.py   # ProposeDisconnectionSkill
├── render_report.py           # RenderReportSkill
├── repair_reaction.py         # RepairReactionSkill
└── validate_reaction.py       # ValidateReactionSkill
```

## SKILL.md 统一格式

每个 SKILL.md 包含:

1. **YAML Frontmatter** — name, description
2. **Purpose** — Skill 目的和作用
3. **When to Use** — 触发场景
4. **Input Parameters** — 参数表
5. **Output** — JSON 输出示例
6. **Key Fields for LLM** — 关键字段解释 (核心 Skill)
7. **Underlying Implementation** — 底层模块路径
8. **CLI Example** — 命令行示例
9. **LLM Workflow** — 使用流程图 (核心 Skill)
10. **Error Handling** — 错误处理说明

## 快速参考

```bash
# 列出所有可用 Skill
python scripts/run_retro.py skill list

# 调用单个 Skill
python scripts/run_retro.py skill analyze_molecule --args '{"smiles": "CCO"}'
python scripts/run_retro.py skill break_bond --args '{"smiles": "CCO", "atom1_idx": 0, "atom2_idx": 1}'
python scripts/run_retro.py skill validate_reaction --args '{"reaction_smiles": "C.C>>CC"}'
```
