# MultiRetro Skill — 逆合成分析工具 (LLM API 参考)

> **版本**: v5.0 | **架构**: `tools/` + `orchestrator/`

## 概述

MultiRetro 是一个多步逆合成规划系统，帮助 LLM 设计有机分子的合成路线。所有操作通过统一 CLI 入口 `scripts/run_retro.py` 完成。

系统提供 **13 个 Skill**，覆盖分子分析、断键、验证、修复、可用性检查、策略决策和报告生成。Skill 执行结果统一封装为 `SkillResult`:

```json
{
  "success": true,
  "data": { ... },
  "error": ""
}
```

### 两种使用模式

| 模式 | 子命令 | 适用场景 |
|------|--------|----------|
| **交互式 Skill 调用** | `skill <name> --args '{...}'` | LLM 逐步决策、精细控制每一步 |
| **任务驱动流水线** | `plan` → `run` → `decide` → `finalize` | 结构化生命周期、多路线对比、LLM 驱动暂停-恢复决策 |

---

## 推荐交互式工作流 (步骤 1-8)

```bash
cd Multi-retro

# 步骤 1: 策略决策 — 确定线性/汇聚/混合策略
python scripts/run_retro.py skill decide_strategy --args '{"smiles":"TARGET_SMILES"}'

# 步骤 2: 分子分析 — 获取 atom_bond_map + 官能团 + SA Score
python scripts/run_retro.py skill analyze_molecule --args '{"smiles":"TARGET_SMILES"}'

# 步骤 3: 精确断键 — 根据 atom_bond_map 中的 strategic_bonds 选择键位
python scripts/run_retro.py skill break_bond --args '{"smiles":"TARGET_SMILES","atom1_idx":X,"atom2_idx":Y}'

# 步骤 4: 验证反应 (硬门控! 必须通过才能继续)
python scripts/run_retro.py skill validate_reaction --args '{"reaction_smiles":"前体>>产物"}'

# 步骤 5: (如验证失败) 修复反应
python scripts/run_retro.py skill repair_reaction --args '{"reaction_smiles":"...","issues":["issue1","issue2"]}'

# 步骤 6: 检查前体可用性
python scripts/run_retro.py skill check_availability --args '{"smiles_list":["前体SMILES1","前体SMILES2"]}'

# 步骤 7: 对每个非终端前体递归执行步骤 2-6

# 步骤 8: 生成合成报告 (必须!)
python scripts/run_retro.py skill render_report --args '{"retro_graph":{...},"output_dir":"outputs/report"}'
```

### 工作流说明

- **步骤 1** 可选但推荐: `decide_strategy` 综合分析分子复杂度，推荐 `linear`/`convergent`/`mixed` 策略
- **步骤 2** 返回的 `atom_bond_map` 包含 `strategic_bonds` 列表，每个键有 `strategic_priority` 评分和 `chemical_label` 语义标注
- **步骤 3** 使用 `break_bond` 而非 `propose_disconnection`，可精确控制断裂位点
- **步骤 4** 是硬门控: `is_valid=false` 时必须修复或换键，不可跳过
- **步骤 5** 最多重试 3 次 (`MAX_REPAIR_RETRIES=3`)，耗尽后应尝试其他断裂位点
- **步骤 6** 基于 SA Score 判断: `< 2.2` 可购买、`< 3.5` 易合成、`≥ 3.5` 需继续分解
- **步骤 7** 递归深度上限 `MAX_DEPTH=7`，分子量 `< 120` 自动标记为终端节点
- **步骤 8** 必须调用，生成完整的合成路线报告

### 辅助 Skill (可在任意步骤穿插使用)

| Skill | 用途 |
|-------|------|
| `analyze_scaffold` | 分析分子骨架、环系统、战略断裂点 |
| `analyze_selectivity` | 化学选择性分析、竞争官能团、保护策略 |
| `map_atoms` | 为反应式添加原子-原子映射 (RXNMapper) |
| `plan_ring_strategy` | 环形成策略规划、构建顺序、逆环开裂 |
| `propose_disconnection` | 自动提议断裂点 (当不需要精确控制时) |
| `get_global_strategy` | 获取全局策略上下文 (PG 策略 + 方法偏好) |

---

## Windows Shell 转义问题 + batch 推荐

SMILES 字符串包含大量 shell 特殊字符: `()[]>>=/#@\`。在 Windows CMD/PowerShell 中直接传递 JSON 参数极易出错。

**推荐方案**: 使用 `batch` 子命令，将任务写入 JSON 文件，彻底绕过 shell 转义:

```bash
# 创建任务文件 tasks.json
python scripts/run_retro.py batch tasks.json
```

`tasks.json` 格式:
```json
[
  {"skill": "analyze_molecule", "args": {"smiles": "CC(=O)Oc1ccccc1C(=O)O"}},
  {"skill": "break_bond", "args": {"smiles": "CC(=O)Oc1ccccc1C(=O)O", "atom1_idx": 3, "atom2_idx": 4}},
  {"skill": "validate_reaction", "args": {"reaction_smiles": "CC(=O)Cl.Oc1ccccc1C(=O)O>>CC(=O)Oc1ccccc1C(=O)O"}}
]
```

### 各平台转义对比

| 平台 | 问题 | 建议 |
|------|------|------|
| Windows CMD | `>` 被解释为重定向 | 使用 `batch` |
| PowerShell | 单引号内的 `"` 需要额外转义 | 使用 `batch` |
| Linux/macOS Bash | 单引号内基本安全 | 可直接使用 `--args '...'` |

---

## 任务驱动模式完整说明

### LLM 驱动暂停-恢复协议

```
plan → 创建路线 + 种子任务 (STRATEGY, ANALYZE, DISCONNECT)
  ↓
run → 执行任务 → 遇到决策点 → 任务状态变为 AWAITING_DECISION → 暂停 (返回 status=paused)
  ↓
输出 DecisionContext JSON (stdout)
  ↓
LLM 解析上下文 → 可选: 调用 explore 探索工具 → 生成 DecisionInstruction
  ↓
decide --decision '{...}' → 恢复执行 → 下一个决策点或完成
  ↓
finalize → 构建 RetroGraph + 合成报告
```

### 完整工作流命令

```bash
# 1. 规划 — 创建 session 目录和初始路线
python scripts/run_retro.py plan --target_smiles "CC(=O)Oc1ccccc1C(=O)O" --output_dir outputs/retro

# 2. 执行 — 运行任务直到遇到决策点暂停
python scripts/run_retro.py run --session_dir outputs/retro --route route_1

# 3. (可选) 探索 — 在决策前调用探索工具获取更多信息
python scripts/run_retro.py explore --session_dir outputs/retro --route route_1 \
  --tool assess_feasibility --params '{"reaction_smiles":"A>>B"}'

# 4. 提交决策
python scripts/run_retro.py decide --session_dir outputs/retro --route route_1 \
  --decision '{"task_id":"1","action":"choose_strategy","params":{"chosen_strategy":"linear"}}'

# 5. 重复步骤 2-4 直到所有任务完成

# 6. 查看状态
python scripts/run_retro.py status --session_dir outputs/retro

# 7. 多路线对比 (如有多条路线)
python scripts/run_retro.py compare --session_dir outputs/retro

# 8. 收尾 — 生成 RetroGraph 和合成报告
python scripts/run_retro.py finalize --session_dir outputs/retro --route route_1
```

### 9 个 CLI 子命令速查

| 子命令 | 说明 | 必需参数 |
|--------|------|----------|
| `skill` | 调用单个 Skill | `skill_name` (位置参数), `--args` (JSON) |
| `plan` | 创建逆合成规划 | `--target_smiles` |
| `run` | 执行任务 | `--session_dir`, `--route` |
| `decide` | 提交决策 | `--session_dir`, `--route`, `--decision` |
| `explore` | 探索工具调用 | `--session_dir`, `--route`, `--tool` |
| `compare` | 多路线对比 | `--session_dir` |
| `finalize` | 收尾生成报告 | `--session_dir`, `--route` |
| `status` | 查看 session 状态 | `--session_dir` |
| `batch` | 批量执行 JSON 任务文件 | `task_file` (位置参数) |


---

## 5 类决策点详解

任务驱动模式中，Orchestrator 在特定任务完成后自动检测决策点，暂停执行并输出 `DecisionContext`。LLM 分析上下文后提交 `DecisionInstruction`。

### 1. strategy_selection (策略选择)

**触发时机**: `STRATEGY` 类型任务完成后 (`decide_strategy` Skill 执行完毕)。

**上下文包含**: 推荐策略 (`recommended_strategy`)、策略评分 (`scores`)、汇聚可行性 (`convergent_feasibility`)、线性可行性 (`linear_feasibility`)、骨架摘要、选择性摘要。

**可选操作**:

| action | 说明 | params 格式 |
|--------|------|-------------|
| `choose_strategy` | 选择合成策略 | `{"chosen_strategy": "linear" \| "convergent" \| "mixed"}` |
| `use_default` | 使用 Skill 推荐的 `recommended_strategy` | 无需 params |

**回退行为** (`use_default`): 直接采用 `result.recommended_strategy`，默认 `"linear"`。

---

### 2. disconnection_decision (断键决策)

**触发时机**: `ANALYZE` 类型任务完成后 (`analyze_molecule` Skill 执行完毕)。

**上下文包含**: `atom_bond_map` (含 `strategic_bonds` 列表)、`retro_guidance` 语义指导、官能团列表、分子摘要。

**可选操作**:

| action | 说明 | params 格式 |
|--------|------|-------------|
| `select_bond` | 选择断裂位点 | `{"atom1_idx": int, "atom2_idx": int}` |
| `custom_reaction` | 提供自定义反应式 | `{"reaction_smiles": str, "precursors": [str, ...]}` |
| `use_default` | 按 `strategic_priority` 降序选第一个键 | 无需 params |

**回退行为** (`use_default`): 将 `strategic_bonds` 按 `strategic_priority` 降序排列，选择第一个键。若无战略键信息，由 `ProposeDisconnectionSkill` 自动选择。

---

### 3. validation_judgment (验证判断)

**触发时机**: `VALIDATE` 类型任务完成后 (`validate_reaction` Skill 执行完毕)。

**上下文包含**: `is_valid` 布尔值、`verdict` 判定文本、`issues` 问题列表、`atom_balance` 原子平衡详情、`reaction_smiles`、`feasibility_assessment` 可行性评估。

**可选操作**:

| action | 说明 | params 格式 |
|--------|------|-------------|
| `accept` | 接受验证结果 | 无需 params |
| `repair` | 创建修复子任务 | 无需 params |
| `retry_other_bond` | 放弃当前断键，回到 DISCONNECT 重新选择 | 无需 params |
| `re_analyze_disconnection` | 回到 ANALYZE 重新分析断裂方案 | 无需 params |
| `accept_with_note` | 带备注接受 (可附加 `reaction_conditions`) | 通过 `DecisionInstruction.reaction_conditions` 传递 |
| `use_default` | `is_valid=true` → accept; `is_valid=false` → repair | 无需 params |

**回退行为** (`use_default`): `is_valid=true` 时自动 `accept`; `is_valid=false` 时自动创建 `repair` 子任务。

---

### 4. repair_judgment (修复判断)

**触发时机**: `REPAIR` 类型任务完成后 (`repair_reaction` Skill 执行完毕)。

**上下文包含**: `repaired_smiles` 修复后的反应式、`repaired_precursors` 修复后的前体、`original_issues` 原始问题、`repair_strategy` 修复策略、`changes_made` 变更说明、`original_reaction_smiles` 原始反应式。

**可选操作**:

| action | 说明 | params 格式 |
|--------|------|-------------|
| `accept_repair` | 接受自动修复结果 | 无需 params |
| `use_custom_repair` | 使用自定义修复 | `{"repaired_reaction_smiles": str, "repaired_precursors": [str, ...]}` |
| `reject_and_retry` | 拒绝修复，重试 | 无需 params |
| `abandon_bond` | 放弃该断裂路径 | 无需 params |
| `use_default` | 有 `repaired_smiles` → accept_repair; 否则 → abandon_bond | 无需 params |

**回退行为** (`use_default`): 若 `result.repaired_smiles` 非空则 `accept_repair`; 否则 `abandon_bond`。

---

### 5. recursion_decision (递归决策)

**触发时机**: `VALIDATE` 通过后展开前体时 (由 `check_recursion_decision_point` 触发)。

**上下文包含**: `precursors` 列表 (每个前体含 `smiles`、`sa_score`、`molecular_weight`、`suggested_action`)、`current_depth` 当前递归深度。

**可选操作**:

| action | 说明 | params 格式 |
|--------|------|-------------|
| `recursion_decision` | 为每个前体指定操作 | 见下方 |
| `use_default` | 基于 SA Score 阈值自动决策 | 无需 params |

`recursion_decision` 的 params 格式:
```json
{
  "precursor_decisions": [
    {"smiles": "CCO", "action": "mark_terminal", "reason": "SA=1.5, MW=46.1, 可购买"},
    {"smiles": "c1ccc(C(=O)Cl)cc1", "action": "decompose", "reason": "SA=4.2, 需继续分解"},
    {"smiles": "CC(=O)O", "action": "mark_optional", "reason": "SA=2.8, 易合成但可选分解"}
  ]
}
```

每个前体的 `action` 可选值:
- `mark_terminal` — 标记为终端节点 (可购买/足够简单)
- `decompose` — 继续递归分解
- `mark_optional` — 标记为可选分解 (SA 介于 2.2~3.5)

**回退行为** (`use_default`):
- `MW < 120` 或 `SA < 2.2` → `mark_terminal`
- `2.2 ≤ SA < 3.5` → `mark_optional`
- `SA ≥ 3.5` → `decompose`

---

## DecisionContext JSON 格式

`DecisionContext` 是 Orchestrator 暂停时输出给 LLM 的完整上下文:

```json
{
  "decision_type": "disconnection_decision",
  "task_id": "1.2",
  "context": {
    "target_smiles": "CC(=O)Oc1ccccc1C(=O)O",
    "atom_bond_map": {"strategic_bonds": [...]},
    "retro_guidance": "建议优先断裂酯键...",
    "functional_groups": [...],
    "summary": {...},
    "route_summary": {
      "route_id": "abc12345",
      "target_smiles": "CC(=O)Oc1ccccc1C(=O)O",
      "completed_steps": 2,
      "total_tasks": 5,
      "current_max_depth": 1
    }
  },
  "available_actions": [
    {"action_id": "select_bond", "description": "选择断裂位点", "params": ["atom1_idx", "atom2_idx"]},
    {"action_id": "custom_reaction", "description": "提供自定义反应式", "params": ["reaction_smiles", "precursors"]},
    {"action_id": "use_default", "description": "使用规则驱动默认选择"}
  ],
  "decision_history": [
    {
      "decision_type": "strategy_selection",
      "task_id": "1",
      "action": "choose_strategy",
      "reasoning_preview": "选择线性策略..."
    }
  ],
  "exploration_tools": [
    {"tool_id": "assess_feasibility", "description": "正向合成可行性评估", "required_params": ["reaction_smiles"], "optional_params": ["reaction_class"]},
    {"tool_id": "analyze_selectivity", "description": "化学选择性分析", "required_params": ["smiles"], "optional_params": []}
  ],
  "exploration_budget": 5
}
```

### DecisionContext 字段说明

| 字段 | 类型 | 说明 |
|------|------|------|
| `decision_type` | `str` | 5 种决策类型之一 |
| `task_id` | `str` | 触发决策的任务 ID |
| `context` | `Dict[str, Any]` | 决策相关的完整上下文 (因决策类型而异) |
| `available_actions` | `List[Dict]` | 可选操作列表 |
| `decision_history` | `List[Dict]` | 历史决策摘要 |
| `exploration_tools` | `List[Dict]` | 当前决策类型可用的探索工具 |
| `exploration_budget` | `int` | 剩余探索调用次数 (默认 5) |

## DecisionInstruction JSON 格式

`DecisionInstruction` 是 LLM 提交给 Orchestrator 的决策指令:

```json
{
  "task_id": "1.2",
  "action": "select_bond",
  "params": {"atom1_idx": 3, "atom2_idx": 4},
  "reasoning": "选择酯键断裂，生成酰氯和醇前体",
  "exploration_log": ["[assess_feasibility] 可行性=high, 置信度=0.78", "[analyze_selectivity] 竞争基团=1, 冲突=0"],
  "reaction_conditions": {
    "reaction_type": "Ester hydrolysis",
    "solvent": "THF",
    "base": "NaOH",
    "temperature": "RT"
  }
}
```

### DecisionInstruction 字段说明

| 字段 | 类型 | 必需 | 说明 |
|------|------|------|------|
| `task_id` | `str` | 是 | 对应 DecisionContext 中的 task_id |
| `action` | `str` | 是 | 选择的操作 (必须在 `available_actions` 中) |
| `params` | `Dict[str, Any]` | 否 | 操作参数 (格式见各决策点说明) |
| `reasoning` | `Optional[str]` | 否 | LLM 的决策理由 |
| `exploration_log` | `Optional[List[str]]` | 否 | 探索工具调用记录 |
| `reaction_conditions` | `Optional[Dict[str, Any]]` | 否 | 反应条件 (用于 `accept_with_note`) |


---

## 6 种探索工具

在决策点暂停期间，LLM 可通过 `explore` 子命令调用探索工具获取额外信息，辅助决策。探索操作是**只读的**，不会改变任务状态。每个决策点有 `MAX_EXPLORATION_CALLS=5` 次调用预算。

### 探索工具列表

#### 1. assess_feasibility — 正向合成可行性评估

- **说明**: 检查官能团兼容性、副反应风险、保护基需求、机理合理性
- **必需参数**: `reaction_smiles` (str)
- **可选参数**: `reaction_class` (str)
- **适用决策**: `disconnection_decision`, `validation_judgment`, `repair_judgment`
- **局限性**: 基于静态规则 (官能团矩阵 + 模式匹配)，非反应机理动态模拟。置信度上限 0.85。LLM 需独立验证关键决策。

```bash
python scripts/run_retro.py explore --session_dir outputs/retro --route route_1 \
  --tool assess_feasibility --params '{"reaction_smiles":"CC(=O)Cl.Oc1ccccc1>>CC(=O)Oc1ccccc1"}'
```

#### 2. analyze_selectivity — 化学选择性分析

- **说明**: 识别竞争官能团、反应性冲突、保护策略
- **必需参数**: `smiles` (str)
- **可选参数**: 无
- **适用决策**: `strategy_selection`, `disconnection_decision`, `validation_judgment`, `recursion_decision`
- **局限性**: 基于 SMARTS 模式匹配，不考虑具体构象和溶剂效应。冲突对列表有限 (12 种)。LLM 需评估具体化学环境。

```bash
python scripts/run_retro.py explore --session_dir outputs/retro --route route_1 \
  --tool analyze_selectivity --params '{"smiles":"CC(=O)Oc1ccccc1C(=O)O"}'
```

#### 3. analyze_molecule — 分子分析

- **说明**: 获取 atom_bond_map、官能团、SA Score
- **必需参数**: `smiles` (str)
- **可选参数**: 无
- **适用决策**: `disconnection_decision`, `recursion_decision`
- **局限性**: 结构分析工具，输出客观数据，局限性较小。

```bash
python scripts/run_retro.py explore --session_dir outputs/retro --route route_1 \
  --tool analyze_molecule --params '{"smiles":"c1ccc(C(=O)Cl)cc1"}'
```

#### 4. validate_reaction — 反应验证

- **说明**: 原子平衡 + 化学合理性检查 + 正向可行性
- **必需参数**: `reaction_smiles` (str)
- **可选参数**: `reaction_class` (str)
- **适用决策**: `disconnection_decision`, `repair_judgment`
- **局限性**: 原子平衡检查是精确的；化学合理性判断基于规则。副产物推断覆盖 18 种模式，可能遗漏非常规副产物。

```bash
python scripts/run_retro.py explore --session_dir outputs/retro --route route_1 \
  --tool validate_reaction --params '{"reaction_smiles":"A>>B"}'
```

#### 5. check_availability — 前体可用性检查

- **说明**: 查询前体是否可购买 (基于 SA Score 估算)
- **必需参数**: `smiles_list` (List[str] 或 JSON 数组字符串)
- **可选参数**: 无
- **适用决策**: `recursion_decision`, `validation_judgment`
- **局限性**: 基于 SA Score 估算，非实时数据库查询。

```bash
python scripts/run_retro.py explore --session_dir outputs/retro --route route_1 \
  --tool check_availability --params '{"smiles_list":["CCO","c1ccccc1"]}'
```

#### 6. recommend_protecting_groups — 保护基推荐

- **说明**: 基于当前活跃 PG 和官能团冲突推荐正交保护基
- **必需参数**: `functional_groups` (List[str] 或 str)
- **可选参数**: 无
- **适用决策**: `disconnection_decision`, `validation_judgment`, `repair_judgment`
- **局限性**: 返回建议性推荐，不会自动在路线中插入保护/脱保护步骤。LLM 需自行评估: (1) 是否真的需要保护; (2) 保护基与后续所有步骤的兼容性; (3) 脱保护时机和条件。

```bash
python scripts/run_retro.py explore --session_dir outputs/retro --route route_1 \
  --tool recommend_protecting_groups --params '{"functional_groups":["hydroxyl","amine"]}'
```

### 探索工具适用性矩阵

| 探索工具 | strategy_selection | disconnection_decision | validation_judgment | repair_judgment | recursion_decision |
|----------|:-:|:-:|:-:|:-:|:-:|
| assess_feasibility | | ✓ | ✓ | ✓ | |
| analyze_selectivity | ✓ | ✓ | ✓ | | ✓ |
| analyze_molecule | | ✓ | | | ✓ |
| validate_reaction | | ✓ | | ✓ | |
| check_availability | | | ✓ | | ✓ |
| recommend_protecting_groups | | ✓ | ✓ | ✓ | |

---

## 错误处理指南

### RetroError 层级

所有错误继承自 `RetroError` 基类:

```python
class RetroError(Exception):
    def __init__(self, code: str, message: str, severity: ErrorSeverity = ErrorSeverity.MEDIUM):
        ...
    def to_dict(self) -> Dict[str, Any]: ...
```

| 异常类 | 说明 | 典型触发场景 |
|--------|------|-------------|
| `RetroError` | 基类 | 所有逆合成错误的父类 |
| `ValidationError` | 反应或输入验证失败 | SMILES 格式错误、原子平衡失败、反应式不合理 |
| `ExecutionError` | Skill 或工具执行失败 | RDKit 调用异常、文件 I/O 错误、超时 |
| `SkillError` | Skill 逻辑错误 | 参数不合法、前置条件未满足、Skill 内部断言失败 |
| `ChemistryError` | RDKit 或化学操作失败 | 分子解析失败、SMARTS 匹配异常、描述符计算错误 |

### ErrorSeverity 严重级别

| 级别 | 值 | 说明 |
|------|-----|------|
| `LOW` | `"low"` | 可忽略的警告 |
| `MEDIUM` | `"medium"` | 需要关注但不阻塞 (默认级别) |
| `HIGH` | `"high"` | 影响结果质量，建议处理 |
| `CRITICAL` | `"critical"` | 阻塞性错误，必须处理 |

### 常见失败场景排查

#### SMILES 解析失败

- **症状**: `ChemistryError` — "Failed to parse SMILES"
- **排查**: 检查 SMILES 格式是否正确 (括号匹配、原子价态、芳香性标记)
- **处理**: 使用 `analyze_molecule` Skill 验证 SMILES 有效性

#### 原子平衡失败

- **症状**: `ValidationError` — `is_valid=false`, `issues` 中包含原子平衡问题
- **排查**: 检查反应式两侧原子数是否匹配，注意隐式氢
- **处理**: 使用 `repair_reaction` Skill 自动修复，或手动调整反应式

#### 修复重试耗尽

- **症状**: 连续 3 次 (`MAX_REPAIR_RETRIES=3`) 修复失败
- **排查**: 原始反应式可能存在根本性问题
- **处理**: 放弃当前断键 (`abandon_bond`)，尝试其他断裂位点

#### 探索预算耗尽

- **症状**: `explore` 命令返回预算不足错误
- **排查**: 已使用 5 次 (`MAX_EXPLORATION_CALLS=5`) 探索调用
- **处理**: 基于已有信息做出决策，或使用 `use_default` 回退

#### 任务总数超限

- **症状**: 任务数达到 50 (`MAX_TOTAL_TASKS=50`)
- **排查**: 递归分解过深或分支过多
- **处理**: 简化策略，将更多前体标记为终端节点，减少递归深度

#### 递归深度超限

- **症状**: 递归深度达到 7 (`MAX_DEPTH=7`)
- **排查**: 分子过于复杂或策略不当
- **处理**: 当前前体自动标记为终端节点，不再继续分解

---

## SA Score 分类阈值 + 资源限制参考

### SA Score 分类阈值 (来源: `SAThresholds`)

| 分类 | 阈值条件 | 说明 |
|------|----------|------|
| **purchasable** (可购买) | `SA < 2.2` | 商业可得，标记为终端节点 |
| **easily_synthesizable** (易合成) | `2.2 ≤ SA < 3.5` | 合成难度低，可选择性分解 |
| **complex** (复杂) | `SA ≥ 3.5` | 需要继续逆合成分解 |

### 资源限制参考 (来源: `RetroLimits`)

| 常量 | 值 | 说明 |
|------|-----|------|
| `MAX_DEPTH` | `7` | 递归分解最大深度 |
| `MW_TERMINAL` | `120.0` | 分子量低于此值自动标记为终端节点 |
| `MAX_REPAIR_RETRIES` | `3` | 单次修复最大重试次数 |
| `MAX_EXPLORATION_CALLS` | `5` | 每个决策点的探索调用预算 |
| `MAX_GRAPH_DEPTH` | `20` | RetroGraph 最大深度 |
| `MAX_TOTAL_TASKS` | `50` | 单条路线最大任务数 |

### 其他常量

| 常量 | 值 | 说明 |
|------|-----|------|
| `CONFIDENCE_THRESHOLD` | `0.3` | 置信度阈值 |
| `DEFAULT_MAX_BYPRODUCTS` | `3` | 默认最大副产物数 |
| `DEFAULT_MAX_PROPOSALS` | `5` | 默认最大断裂提议数 |


---

## 每个 Skill 的简要接口

> 详细接口规范请参阅各 Skill 目录下的 `SKILL.md` 文件 (`skills/<name>/SKILL.md`)。

### analyze_molecule — 分子结构分析

- **类名**: `AnalyzeMoleculeSkill`
- **说明**: Analyze target molecule structure, functional groups, properties
- **输入**: `smiles` (str, 必需) — 待分析分子的 SMILES
- **输出**: `AnalysisResult` — 含 `molecule_info` (MoleculeInfo: smiles, canonical_smiles, molecular_weight, sa_score, atom_bond_map, functional_groups, protecting_groups)
- **CLI**: `python scripts/run_retro.py skill analyze_molecule --args '{"smiles":"CCO"}'`

### analyze_scaffold — 分子骨架分析

- **类名**: `AnalyzeScaffoldSkill`
- **说明**: Analyze molecular scaffold for strategic disconnections and ring systems
- **输入**: `smiles` (str, 必需) — 分子 SMILES
- **输出**: 骨架分析结果 (环系统、战略断裂点)
- **CLI**: `python scripts/run_retro.py skill analyze_scaffold --args '{"smiles":"c1ccccc1"}'`

### analyze_selectivity — 化学选择性分析

- **类名**: `AnalyzeSelectivitySkill`
- **说明**: Analyze chemoselectivity, competing FGs, and protection requirements
- **输入**: `smiles` (str, 必需) — 分子 SMILES
- **输出**: 选择性分析报告 (竞争官能团、反应性冲突、保护策略)
- **CLI**: `python scripts/run_retro.py skill analyze_selectivity --args '{"smiles":"CC(=O)Oc1ccccc1O"}'`

### break_bond — 精确断键

- **类名**: `BreakBondSkill`
- **说明**: Break a bond by atom indices to generate precursor SMILES automatically
- **输入**:
  - `smiles` (str, 必需) — 目标分子 SMILES
  - `atom1_idx` (int, 必需) — 断裂键的第一个原子索引
  - `atom2_idx` (int, 必需) — 断裂键的第二个原子索引
- **输出**: `DisconnectionResult` — 含 `fragments`, `reaction_smiles`, `retro_analysis_guide`
- **CLI**: `python scripts/run_retro.py skill break_bond --args '{"smiles":"CCO","atom1_idx":0,"atom2_idx":1}'`

### check_availability — 前体可用性检查

- **类名**: `CheckAvailabilitySkill`
- **说明**: Check precursor availability based on SA score
- **输入**: `smiles_list` (List[str] | str, 必需) — SMILES 列表 (JSON 数组或逗号分隔)
- **输出**: `AvailabilityResult` — 含 `precursors` 列表 (每个含可用性信息)
- **CLI**: `python scripts/run_retro.py skill check_availability --args '{"smiles_list":["CCO","c1ccccc1"]}'`

### decide_strategy — 综合策略决策

- **类名**: `DecideStrategySkill`
- **说明**: Comprehensive strategy decision (linear/convergent) based on all analyses
- **输入**: `smiles` (str, 必需) — 目标分子 SMILES
- **输出**: 策略推荐 (含 `recommended_strategy`, `scores`, `rationale`)
- **CLI**: `python scripts/run_retro.py skill decide_strategy --args '{"smiles":"CC(=O)Oc1ccccc1C(=O)O"}'`

### get_global_strategy — 全局策略上下文

- **类名**: `GetGlobalStrategySkill`
- **说明**: Get global strategy context for LLM prompts with PG policy and method preferences
- **输入**:
  - `target_smiles` (str, 必需) — 目标分子 SMILES
  - `max_pg_operations` (int, 可选, 默认 2) — 最大保护基操作数
  - `max_steps` (int, 可选, 默认 10) — 最大合成步骤数
  - `constraints` (Dict[str, Any], 可选) — 额外约束
- **输出**: 全局策略上下文 (PG 策略 + 方法偏好)
- **CLI**: `python scripts/run_retro.py skill get_global_strategy --args '{"target_smiles":"CCO"}'`

### map_atoms — 原子映射

- **类名**: `MapAtomsSkill`
- **说明**: Add atom-to-atom mapping using RXNMapper
- **输入**: `reaction_smiles` (str, 必需) — 反应 SMILES
- **输出**: 含原子映射的反应 SMILES
- **CLI**: `python scripts/run_retro.py skill map_atoms --args '{"reaction_smiles":"CCO>>CC=O"}'`

### plan_ring_strategy — 环策略规划

- **类名**: `PlanRingStrategySkill`
- **说明**: Plan ring formation strategy including construction order and retro-ring-opening
- **输入**: `smiles` (str, 必需) — 分子 SMILES
- **输出**: 环形成策略 (构建顺序、逆环开裂方案)
- **CLI**: `python scripts/run_retro.py skill plan_ring_strategy --args '{"smiles":"c1ccccc1"}'`

### propose_disconnection — 断裂提议

- **类名**: `ProposeDisconnectionSkill`
- **说明**: Propose retrosynthetic disconnection points
- **输入**:
  - `target_smiles` (str, 必需) — 目标分子 SMILES
  - `strategy_hint` (str, 可选) — 策略过滤提示
  - `max_proposals` (int, 可选, 默认 5) — 最大提议数
- **输出**: 断裂提议列表
- **CLI**: `python scripts/run_retro.py skill propose_disconnection --args '{"target_smiles":"CCO"}'`

### render_report — 合成报告生成

- **类名**: `RenderReportSkill`
- **说明**: Generate synthesis report from completed RetroGraph
- **输入**:
  - `retro_graph` (dict | RetroGraph, 必需) — 完成的 RetroGraph
  - `output_dir` (str, 可选, 默认 `"outputs/report"`) — 输出目录
- **输出**: 合成报告 (Markdown + JSON)
- **CLI**: `python scripts/run_retro.py skill render_report --args '{"retro_graph":{...},"output_dir":"outputs/report"}'`

### repair_reaction — 反应修复

- **类名**: `RepairReactionSkill`
- **说明**: Repair reaction by analysing issues and suggesting fixes
- **输入**:
  - `reaction_smiles` (str, 必需) — 待修复的反应 SMILES
  - `issues` (List[str | Dict[str, str]], 必需) — 问题列表 (字符串或 `{code, message}` 字典)
- **输出**: `RepairResult` — 含 `repaired_smiles`, `repair_actions`
- **CLI**: `python scripts/run_retro.py skill repair_reaction --args '{"reaction_smiles":"A>>B","issues":["atom_balance"]}'`

### validate_reaction — 反应验证

- **类名**: `ValidateReactionSkill`
- **说明**: Validate reaction for atom balance and chemical plausibility
- **输入**:
  - `reaction_smiles` (str, 必需) — 反应 SMILES
  - `reaction_class` (str, 可选) — 反应类型提示
- **输出**: `ValidationResult` — 含 `is_valid`, `verdict`, `issues`, `byproducts`
- **CLI**: `python scripts/run_retro.py skill validate_reaction --args '{"reaction_smiles":"CCO>>CC=O"}'`

---

## 附录: Skill 名称速查表

| CLI 名称 (name 属性) | Python 类名 | 源文件 |
|----------------------|-------------|--------|
| `analyze_molecule` | `AnalyzeMoleculeSkill` | `tools/skills/analyze_molecule.py` |
| `analyze_scaffold` | `AnalyzeScaffoldSkill` | `tools/skills/analyze_scaffold.py` |
| `analyze_selectivity` | `AnalyzeSelectivitySkill` | `tools/skills/analyze_selectivity.py` |
| `break_bond` | `BreakBondSkill` | `tools/skills/break_bond.py` |
| `check_availability` | `CheckAvailabilitySkill` | `tools/skills/check_availability.py` |
| `decide_strategy` | `DecideStrategySkill` | `tools/skills/decide_strategy.py` |
| `get_global_strategy` | `GetGlobalStrategySkill` | `tools/skills/get_global_strategy.py` |
| `map_atoms` | `MapAtomsSkill` | `tools/skills/map_atoms.py` |
| `plan_ring_strategy` | `PlanRingStrategySkill` | `tools/skills/plan_ring_strategy.py` |
| `propose_disconnection` | `ProposeDisconnectionSkill` | `tools/skills/propose_disconnection.py` |
| `render_report` | `RenderReportSkill` | `tools/skills/render_report.py` |
| `repair_reaction` | `RepairReactionSkill` | `tools/skills/repair_reaction.py` |
| `validate_reaction` | `ValidateReactionSkill` | `tools/skills/validate_reaction.py` |
