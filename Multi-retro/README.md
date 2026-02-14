# Multi-retro — LLM 驱动的多路线逆合成分析平台

Multi-retro 是一个面向 LLM 的逆合成分析 Skill 项目，通过 13 个专业 Skill 和暂停-恢复决策协议，支持 LLM 自主规划、执行和优化多条合成路线。

---

## 快速开始

### 环境安装

```bash
pip install -r requirements.txt
```

### 基本用法

```bash
# 1. 分析目标分子
python scripts/run_retro.py skill analyze_molecule --args '{"smiles": "CC(=O)Oc1ccccc1C(=O)O"}'

# 2. 创建逆合成计划
python scripts/run_retro.py plan --target_smiles "CC(=O)Oc1ccccc1C(=O)O"

# 3. 执行任务
python scripts/run_retro.py run --session_dir outputs/retro --route route_001

# 4. 查看状态
python scripts/run_retro.py status --session_dir outputs/retro
```

---

## 统一 CLI 入口

所有操作通过 `scripts/run_retro.py` 的 9 个子命令完成，替代了旧的多脚本架构。

### 子命令速查表

| 子命令 | 说明 | 关键参数 |
|--------|------|----------|
| `skill` | 调用单个 Skill | `skill_name` (位置参数), `--args` |
| `plan` | 创建逆合成计划 | `--target_smiles` (必需), `--output_dir` |
| `run` | 执行会话中的任务 | `--session_dir` (必需), `--route` (必需) |
| `decide` | 提交决策（恢复暂停的会话） | `--session_dir`, `--route`, `--decision` (均必需) |
| `explore` | 决策点处调用探索工具 | `--session_dir`, `--route`, `--tool`, `--params` |
| `compare` | 多路线对比 | `--session_dir` (必需) |
| `finalize` | 生成最终报告 | `--session_dir` (必需), `--route` (必需) |
| `status` | 查看会话状态 | `--session_dir` (必需) |
| `batch` | 从 JSON 文件批量执行 | `task_file` (位置参数) |

### 子命令详细用法

#### skill — 调用单个 Skill

```bash
# 列出所有可用 Skill
python scripts/run_retro.py skill list

# 调用指定 Skill
python scripts/run_retro.py skill analyze_molecule --args '{"smiles": "c1ccccc1"}'
```

参数:
- `skill_name` (位置参数, 必需): Skill 名称或 `list`
- `--args` (可选, 默认 None): Skill 参数的 JSON 字符串

#### plan — 创建逆合成计划

```bash
python scripts/run_retro.py plan --target_smiles "CC(=O)Oc1ccccc1C(=O)O" --output_dir outputs/retro
```

参数:
- `--target_smiles` (必需): 目标分子 SMILES
- `--output_dir` (可选, 默认 `outputs/retro`): 输出目录

#### run — 执行任务

```bash
python scripts/run_retro.py run --session_dir outputs/retro --route route_001
```

参数:
- `--session_dir` (必需): 会话目录路径
- `--route` (必需): 路线 ID

#### decide — 提交决策

```bash
python scripts/run_retro.py decide --session_dir outputs/retro --route route_001 \
  --decision '{"task_id": "task_003", "action": "accept", "params": {}}'
```

参数:
- `--session_dir` (必需): 会话目录路径
- `--route` (必需): 路线 ID
- `--decision` (必需): 决策 JSON 字符串，或 `@filepath` 从文件加载

#### explore — 探索工具调用

```bash
python scripts/run_retro.py explore --session_dir outputs/retro --route route_001 \
  --tool assess_feasibility --params '{"reaction_smiles": "CC>>C.C"}'
```

参数:
- `--session_dir` (必需): 会话目录路径
- `--route` (必需): 路线 ID
- `--tool` (必需): 探索工具 ID
- `--params` (可选, 默认 None): 工具参数 JSON 字符串

#### compare — 多路线对比

```bash
python scripts/run_retro.py compare --session_dir outputs/retro
```

参数:
- `--session_dir` (必需): 会话目录路径

#### finalize — 生成报告

```bash
python scripts/run_retro.py finalize --session_dir outputs/retro --route route_001
```

参数:
- `--session_dir` (必需): 会话目录路径
- `--route` (必需): 路线 ID

#### status — 查看状态

```bash
python scripts/run_retro.py status --session_dir outputs/retro
```

参数:
- `--session_dir` (必需): 会话目录路径

#### batch — 批量执行

```bash
python scripts/run_retro.py batch tasks.json
```

参数:
- `task_file` (位置参数, 必需): JSON 任务文件路径

---

## LLM 决策协议

Multi-retro 采用暂停-恢复决策协议，在关键节点暂停执行并请求 LLM 做出决策。

### 决策循环

```
plan → run → [暂停: 返回 DecisionContext] → LLM 探索 (可选) → decide → run → ...
```

当 `run` 遇到需要决策的任务时，执行暂停并返回 `DecisionContext`。LLM 可以：
1. 使用 `explore` 子命令调用探索工具收集更多信息（预算: MAX_EXPLORATION_CALLS=5）
2. 构造 `DecisionInstruction` 并通过 `decide` 子命令提交

### 5 类决策点

| 决策类型 | 触发时机 | 可选操作 |
|---------|---------|---------|
| `strategy_selection` | 分析完成后选择合成策略 | `linear`, `convergent`, `use_default` |
| `disconnection_decision` | 选择断键位点 | `select_bond`, `use_default`, `skip` |
| `validation_judgment` | 反应验证后判断 | `accept`, `repair`, `retry`, `use_default` |
| `repair_judgment` | 自动修复后判断 | `accept`, `reject`, `use_default` |
| `recursion_decision` | 前体展开决策 | `expand`, `terminate`, `use_default` |

### DecisionContext JSON 格式

LLM 收到的决策上下文:

```json
{
  "decision_type": "disconnection_decision",
  "task_id": "task_003",
  "context": {
    "target_smiles": "...",
    "analysis_result": { "..." : "..." }
  },
  "available_actions": [
    {"action": "select_bond", "params": {"atom1_idx": "int", "atom2_idx": "int"}},
    {"action": "use_default", "params": {}},
    {"action": "skip", "params": {}}
  ],
  "decision_history": [],
  "exploration_tools": [
    {"tool_id": "assess_feasibility", "required_params": ["reaction_smiles"]},
    {"tool_id": "analyze_molecule", "required_params": ["smiles"]}
  ],
  "exploration_budget": 5
}
```

### DecisionInstruction JSON 格式

LLM 提交的决策指令:

```json
{
  "task_id": "task_003",
  "action": "select_bond",
  "params": {
    "atom1_idx": 2,
    "atom2_idx": 5
  },
  "reasoning": "选择酯键断裂，生成两个简单前体",
  "exploration_log": ["assess_feasibility: 可行性 0.82"],
  "reaction_conditions": {
    "reagent": "NaOH",
    "solvent": "THF"
  }
}
```

### 6 种探索工具

在决策点处，LLM 可调用以下探索工具收集信息（只读，不改变任务状态）:

| 工具 ID | 说明 | 必需参数 | 适用决策类型 |
|---------|------|---------|-------------|
| `assess_feasibility` | 正向合成可行性评估 | `reaction_smiles` | disconnection, validation, repair |
| `analyze_selectivity` | 化学选择性分析 | `smiles` | strategy, disconnection, validation, recursion |
| `analyze_molecule` | 分子分析 (atom_bond_map, 官能团, SA Score) | `smiles` | disconnection, recursion |
| `validate_reaction` | 反应验证 (原子平衡 + 合理性) | `reaction_smiles` | disconnection, repair |
| `check_availability` | 前体可用性检查 | `smiles_list` | recursion, validation |
| `recommend_protecting_groups` | 保护基推荐 | `functional_groups` | disconnection, validation, repair |

### 决策回退机制

每种决策类型都支持 `use_default` 操作，当 LLM 无法做出决策时自动使用默认策略:
- `strategy_selection` → 默认选择 `linear` 策略
- `disconnection_decision` → 选择 `strategic_priority` 最高的可断键
- `validation_judgment` → 有效则 `accept`，无效则 `repair`
- `repair_judgment` → 修复成功则 `accept`，否则 `reject`
- `recursion_decision` → SA Score < 2.2 或 MW < 120 则 `terminate`，否则 `expand`

---

## 13 个 Skill 速查表

| Skill 名称 | 说明 | 关键输入 |
|------------|------|---------|
| `analyze_molecule` | 分析目标分子结构、官能团、性质 | `smiles` |
| `analyze_scaffold` | 分析分子骨架的战略断裂点和环系统 | `smiles` |
| `analyze_selectivity` | 分析化学选择性、竞争官能团、保护需求 | `smiles` |
| `break_bond` | 按原子索引断键生成前体 SMILES | `smiles`, `atom1_idx`, `atom2_idx` |
| `check_availability` | 基于 SA Score 检查前体可用性 | `smiles_list` |
| `decide_strategy` | 综合策略决策 (线性/汇聚) | `smiles` |
| `get_global_strategy` | 获取全局策略上下文 (PG 策略、方法偏好) | `target_smiles` |
| `map_atoms` | 使用 RXNMapper 添加原子映射 | `reaction_smiles` |
| `plan_ring_strategy` | 规划环形成策略 (构建顺序、逆环开裂) | `smiles` |
| `propose_disconnection` | 提出逆合成断裂点 | `target_smiles` |
| `render_report` | 从 RetroGraph 生成合成报告 | `retro_graph` |
| `repair_reaction` | 分析问题并修复反应 | `reaction_smiles`, `issues` |
| `validate_reaction` | 验证反应的原子平衡和化学合理性 | `reaction_smiles` |

> 每个 Skill 的详细接口规范见 `skills/<skill_folder>/SKILL.md`，完整 API 参考见 [SKILL.md](SKILL.md)。

---

## 项目架构

```
Multi-retro/
├── README.md                          # 项目入口 (本文件)
├── SKILL.md                           # Skill API 完整参考 (LLM 面向)
├── requirements.txt                   # Python 依赖
│
├── orchestrator/                      # 编排器
│   ├── __init__.py
│   └── orchestrator.py                # RetroOrchestrator: plan/run_all/decide/finalize
│
├── tools/                             # 核心工具库
│   ├── __init__.py
│   ├── chem/                          # 化学计算模块
│   │   ├── atom_mapper.py             # 原子映射 (RXNMapper)
│   │   ├── functional_groups.py       # 官能团识别
│   │   ├── mol_parser.py              # SMILES 解析
│   │   ├── reaction_classifier.py     # 反应分类
│   │   ├── reaction_validator.py      # 反应验证
│   │   ├── sa_scorer.py               # SA Score 计算
│   │   └── structure_analyzer.py      # 结构分析
│   ├── common/                        # 公共定义
│   │   ├── constants.py               # 常量 (SAThresholds, RetroLimits)
│   │   ├── errors.py                  # 错误层级 (RetroError → 4 子类)
│   │   └── status.py                  # 枚举 (TaskStatus, TaskType, RouteStatus)
│   ├── guidance/                      # 指导与建议模块
│   │   ├── disconnection_proposer.py  # 断裂点提议
│   │   ├── feasibility_assessor.py    # 可行性评估
│   │   ├── repair_advisor.py          # 修复建议
│   │   ├── retro_guide.py             # 逆合成指导
│   │   ├── selectivity_analyzer.py    # 选择性分析
│   │   └── strategy_advisor.py        # 策略建议
│   ├── models/                        # 数据模型
│   │   ├── chem_models.py             # AtomInfo, BondInfo, RingInfo, MoleculeInfo
│   │   ├── output_models.py           # RetroGraphNode, RetroGraphEdge, RetroGraph
│   │   ├── skill_models.py            # AnalyzeArgs, ValidateArgs, RepairArgs 等
│   │   └── workflow_models.py         # RetroTask, RetroRoute, DecisionContext 等
│   ├── output/                        # 输出与可视化
│   │   ├── graph_builder.py           # 图构建
│   │   ├── journal.py                 # 日志记录
│   │   ├── report_generator.py        # 报告生成
│   │   ├── status_renderer.py         # 状态渲染
│   │   └── visualizer.py              # 可视化
│   ├── skills/                        # Skill 实现 (13 个)
│   │   ├── base.py                    # BaseSkill, SkillResult
│   │   ├── analyze_molecule.py
│   │   ├── analyze_scaffold.py
│   │   ├── analyze_selectivity.py
│   │   ├── break_bond.py
│   │   ├── check_availability.py
│   │   ├── decide_strategy.py
│   │   ├── get_global_strategy.py
│   │   ├── map_atoms.py
│   │   ├── plan_ring_strategy.py
│   │   ├── propose_disconnection.py
│   │   ├── render_report.py
│   │   ├── repair_reaction.py
│   │   └── validate_reaction.py
│   └── workflow/                      # 工作流管理
│       ├── decision_manager.py        # DecisionManager: 5 类决策构建/应用/回退
│       ├── exploration.py             # ExplorationSession: 6 种探索工具
│       ├── route_manager.py           # RouteManager: 路线创建/对比/持久化
│       ├── skill_dispatch.py          # dispatch_skill: Skill 分发
│       ├── snapshot.py                # save/load_snapshot: 快照持久化
│       ├── subtask_manager.py         # 子任务创建与展开
│       └── task_state.py              # transition_status: 状态转换规则
│
├── scripts/                           # CLI 入口
│   └── run_retro.py                   # 统一 CLI (9 个子命令)
│
├── schemas/                           # JSON Schema 定义
│   ├── host_request.schema.json       # 宿主请求格式
│   ├── routes.schema.json             # 路线数据格式
│   ├── strategy.schema.json           # 策略数据格式
│   └── vis_plan.schema.json           # 可视化计划格式
│
├── skills/                            # Skill 文档 (每个 Skill 一个文件夹)
│   ├── analyze_molecule/SKILL.md
│   ├── analyze_scaffold/SKILL.md
│   ├── analyze_selectivity/SKILL.md
│   ├── check_precursor_availability/SKILL.md
│   ├── decide_strategy/SKILL.md
│   ├── map_reaction_atoms/SKILL.md
│   ├── multi_turn_render/SKILL.md
│   ├── plan_ring_strategy/SKILL.md
│   ├── propose_disconnection/SKILL.md
│   ├── render_synthesis_report/SKILL.md
│   ├── repair_reaction/SKILL.md
│   └── validate_reaction/SKILL.md
│
├── outputs/                           # 运行输出目录
│
└── tests/                             # 测试套件
    ├── conftest.py
    ├── chem/                          # 化学模块测试
    ├── docs/                          # 文档属性测试
    ├── guidance/                      # 指导模块测试
    ├── output/                        # 输出模块测试
    ├── skills/                        # Skill 测试
    └── workflow/                      # 工作流测试
```

---

## 枚举与常量速查表

### TaskStatus (8 个成员)

任务生命周期状态，定义于 `tools/common/status.py`:

| 成员 | 值 | 说明 |
|------|-----|------|
| `PENDING` | `"pending"` | 等待执行 |
| `IN_PROGRESS` | `"in_progress"` | 执行中 |
| `AWAITING_DECISION` | `"awaiting_decision"` | 等待 LLM 决策 |
| `VALIDATED` | `"validated"` | 已验证通过 |
| `COMPLETED` | `"completed"` | 已完成 |
| `FAILED` | `"failed"` | 执行失败 |
| `SKIPPED` | `"skipped"` | 已跳过 |
| `BLOCKED` | `"blocked"` | 被阻塞 |

### TaskType (7 个成员)

任务类型，映射到具体 Skill，定义于 `tools/common/status.py`:

| 成员 | 值 |
|------|-----|
| `STRATEGY` | `"strategy"` |
| `ANALYZE` | `"analyze"` |
| `DISCONNECT` | `"disconnect"` |
| `VALIDATE` | `"validate"` |
| `REPAIR` | `"repair"` |
| `AVAILABILITY` | `"availability"` |
| `REPORT` | `"report"` |

### RouteStatus (5 个成员)

合成路线状态，定义于 `tools/common/status.py`:

| 成员 | 值 |
|------|-----|
| `PLANNING` | `"planning"` |
| `COMPLETED` | `"completed"` |
| `ABANDONED` | `"abandoned"` |
| `PARTIAL` | `"partial"` |
| `FAILED` | `"failed"` |

### ErrorSeverity (4 个成员)

错误严重程度，定义于 `tools/common/errors.py`:

| 成员 | 值 |
|------|-----|
| `LOW` | `"low"` |
| `MEDIUM` | `"medium"` |
| `HIGH` | `"high"` |
| `CRITICAL` | `"critical"` |

### SAThresholds — SA Score 分类阈值

定义于 `tools/common/constants.py`:

| 常量 | 值 | 含义 |
|------|-----|------|
| `PURCHASABLE_MAX` | `2.2` | SA Score < 2.2 → 可购买 |
| `EASILY_SYNTHESIZABLE_MAX` | `3.5` | SA Score < 3.5 → 易合成 |
| `COMPLEX_MIN` | `3.5` | SA Score ≥ 3.5 → 复杂分子 |

### RetroLimits — 递归与资源限制

定义于 `tools/common/constants.py`:

| 常量 | 值 | 含义 |
|------|-----|------|
| `MAX_DEPTH` | `7` | 最大递归深度 |
| `MW_TERMINAL` | `120.0` | 分子量终止阈值 |
| `MAX_REPAIR_RETRIES` | `3` | 最大修复重试次数 |
| `MAX_EXPLORATION_CALLS` | `5` | 每个决策点最大探索调用次数 |
| `MAX_GRAPH_DEPTH` | `20` | 最大图深度 |
| `MAX_TOTAL_TASKS` | `50` | 最大任务总数 |

### 模块级常量

定义于 `tools/common/constants.py`:

| 常量 | 值 | 含义 |
|------|-----|------|
| `CONFIDENCE_THRESHOLD` | `0.3` | 置信度阈值 |
| `DEFAULT_MAX_BYPRODUCTS` | `3` | 默认最大副产物数 |
| `DEFAULT_MAX_PROPOSALS` | `5` | 默认最大提议数 |

### 错误层级

定义于 `tools/common/errors.py`:

```
RetroError(code: str, message: str, severity: ErrorSeverity = ErrorSeverity.MEDIUM)
├── ValidationError    — 反应或输入验证失败
├── ExecutionError     — Skill 或工具执行失败
├── SkillError         — Skill 逻辑错误 (参数错误、前置条件不满足)
└── ChemistryError     — RDKit 或化学操作失败
```

`RetroError` 提供 `to_dict() -> Dict[str, Any]` 方法用于序列化。

---

## 核心数据类速查表

### 化学模型 (`tools/models/chem_models.py`)

#### AtomInfo

| 字段 | 类型 | 默认值 |
|------|------|--------|
| `idx` | `int` | `0` |
| `symbol` | `str` | `""` |
| `aromatic` | `bool` | `False` |
| `in_ring` | `bool` | `False` |
| `degree` | `int` | `0` |
| `num_hs` | `int` | `0` |
| `neighbors` | `List[int]` | `[]` |

#### BondInfo

| 字段 | 类型 | 默认值 |
|------|------|--------|
| `idx` | `int` | `0` |
| `atom1_idx` | `int` | `0` |
| `atom2_idx` | `int` | `0` |
| `bond_type` | `str` | `"SINGLE"` |
| `in_ring` | `bool` | `False` |
| `breakable` | `bool` | `False` |
| `chemical_label` | `str` | `""` |
| `potential_reactions` | `List[str]` | `[]` |
| `retro_hint` | `str` | `""` |
| `strategic_priority` | `int` | `3` |

#### RingInfo

| 字段 | 类型 | 默认值 |
|------|------|--------|
| `ring_id` | `int` | `0` |
| `size` | `int` | `0` |
| `atom_indices` | `List[int]` | `[]` |
| `is_aromatic` | `bool` | `False` |
| `ring_type` | `str` | `""` |
| `formation_hints` | `List[str]` | `[]` |

#### AtomBondMap

| 字段 | 类型 | 默认值 |
|------|------|--------|
| `canonical_smiles` | `str` | `""` |
| `atoms` | `List[AtomInfo]` | `[]` |
| `bonds` | `List[BondInfo]` | `[]` |
| `ring_info` | `List[RingInfo]` | `[]` |
| `retro_guidance` | `str` | `""` |

#### MoleculeInfo

| 字段 | 类型 | 默认值 |
|------|------|--------|
| `smiles` | `str` | `""` |
| `canonical_smiles` | `str` | `""` |
| `molecular_weight` | `float` | `0.0` |
| `sa_score` | `float` | `0.0` |
| `atom_bond_map` | `Optional[AtomBondMap]` | `None` |
| `functional_groups` | `Dict[str, int]` | `{}` |
| `protecting_groups` | `Dict[str, Any]` | `{}` |

### Skill 模型 (`tools/models/skill_models.py`)

#### AnalyzeArgs

| 字段 | 类型 | 默认值 |
|------|------|--------|
| `smiles` | `str` | `""` |

#### AnalysisResult

| 字段 | 类型 | 默认值 |
|------|------|--------|
| `success` | `bool` | `False` |
| `molecule_info` | `Optional[MoleculeInfo]` | `None` |
| `error` | `str` | `""` |

#### BreakBondArgs

| 字段 | 类型 | 默认值 |
|------|------|--------|
| `smiles` | `str` | `""` |
| `atom1_idx` | `int` | `0` |
| `atom2_idx` | `int` | `0` |

#### DisconnectionResult

| 字段 | 类型 | 默认值 |
|------|------|--------|
| `success` | `bool` | `False` |
| `fragments` | `List[str]` | `[]` |
| `reaction_smiles` | `str` | `""` |
| `retro_analysis_guide` | `Optional[Dict[str, Any]]` | `None` |
| `error` | `str` | `""` |

#### ValidateArgs

| 字段 | 类型 | 默认值 |
|------|------|--------|
| `reaction_smiles` | `str` | `""` |

#### ValidationResult

| 字段 | 类型 | 默认值 |
|------|------|--------|
| `success` | `bool` | `False` |
| `verdict` | `str` | `""` |
| `is_valid` | `bool` | `False` |
| `issues` | `List[Dict[str, Any]]` | `[]` |
| `byproducts` | `List[Dict[str, Any]]` | `[]` |
| `error` | `str` | `""` |

#### RepairArgs

| 字段 | 类型 | 默认值 |
|------|------|--------|
| `reaction_smiles` | `str` | `""` |
| `issues` | `List[str]` | `[]` |

#### RepairResult

| 字段 | 类型 | 默认值 |
|------|------|--------|
| `success` | `bool` | `False` |
| `repaired_smiles` | `str` | `""` |
| `repair_actions` | `List[str]` | `[]` |
| `error` | `str` | `""` |

#### AvailabilityResult

| 字段 | 类型 | 默认值 |
|------|------|--------|
| `success` | `bool` | `False` |
| `precursors` | `List[Dict[str, Any]]` | `[]` |
| `error` | `str` | `""` |

### 工作流模型 (`tools/models/workflow_models.py`)

#### RetroTask

| 字段 | 类型 | 默认值 |
|------|------|--------|
| `task_id` | `str` | `""` |
| `task_type` | `TaskType` | `TaskType.ANALYZE` |
| `status` | `TaskStatus` | `TaskStatus.PENDING` |
| `target_smiles` | `str` | `""` |
| `parent_task_id` | `Optional[str]` | `None` |
| `depth` | `int` | `0` |
| `result` | `Dict[str, Any]` | `{}` |
| `retry_count` | `int` | `0` |
| `max_retries` | `int` | `3` |
| `selected_bond` | `Optional[Dict[str, Any]]` | `None` |
| `reaction_smiles` | `str` | `""` |
| `precursors` | `List[str]` | `[]` |

#### RetroTaskList

| 字段 | 类型 | 默认值 |
|------|------|--------|
| `tasks` | `List[RetroTask]` | `[]` |

#### RetroRoute

| 字段 | 类型 | 默认值 |
|------|------|--------|
| `route_id` | `str` | `""` |
| `target_smiles` | `str` | `""` |
| `tasks` | `RetroTaskList` | `RetroTaskList()` |
| `metadata` | `Dict[str, Any]` | `{}` |

#### DecisionContext

| 字段 | 类型 | 默认值 |
|------|------|--------|
| `decision_type` | `str` | `""` |
| `task_id` | `str` | `""` |
| `context` | `Dict[str, Any]` | `{}` |
| `available_actions` | `List[Dict[str, Any]]` | `[]` |
| `decision_history` | `List[Dict[str, Any]]` | `[]` |
| `exploration_tools` | `List[Dict[str, Any]]` | `[]` |
| `exploration_budget` | `int` | `5` |

#### DecisionInstruction

| 字段 | 类型 | 默认值 |
|------|------|--------|
| `task_id` | `str` | `""` |
| `action` | `str` | `""` |
| `params` | `Dict[str, Any]` | `{}` |
| `reasoning` | `Optional[str]` | `None` |
| `exploration_log` | `Optional[List[str]]` | `None` |
| `reaction_conditions` | `Optional[Dict[str, Any]]` | `None` |

### 输出模型 (`tools/models/output_models.py`)

#### RetroGraphNode

| 字段 | 类型 | 默认值 |
|------|------|--------|
| `node_id` | `str` | `""` |
| `smiles` | `str` | `""` |
| `node_type` | `str` | `""` |
| `status` | `str` | `""` |
| `sa_score` | `Optional[float]` | `None` |
| `is_terminal` | `bool` | `False` |

#### RetroGraphEdge

| 字段 | 类型 | 默认值 |
|------|------|--------|
| `edge_id` | `str` | `""` |
| `product_id` | `str` | `""` |
| `precursor_ids` | `List[str]` | `[]` |
| `reaction_smiles` | `str` | `""` |
| `reaction_type` | `str` | `""` |
| `conditions` | `Optional[Dict[str, Any]]` | `None` |

#### RetroGraph

| 字段 | 类型 | 默认值 |
|------|------|--------|
| `target_smiles` | `str` | `""` |
| `nodes` | `Dict[str, RetroGraphNode]` | `{}` |
| `edges` | `Dict[str, RetroGraphEdge]` | `{}` |
| `metadata` | `Dict[str, Any]` | `{}` |

### SkillResult (`tools/skills/base.py`)

所有 Skill 的标准返回信封:

| 字段 | 类型 | 默认值 |
|------|------|--------|
| `success` | `bool` | *(必需)* |
| `data` | `Dict[str, Any]` | `{}` |
| `error` | `str` | `""` |

---

## 核心管理器速查表

### RetroOrchestrator (`orchestrator/orchestrator.py`)

顶层编排器，协调整个逆合成工作流:

| 方法 | 说明 |
|------|------|
| `plan(target_smiles: str) -> RetroRoute` | 创建逆合成计划，生成初始任务列表 |
| `run_all(route_id: str) -> None` | 执行路线中的所有任务，遇到决策点暂停 |
| `decide(route_id: str, instruction: DecisionInstruction) -> None` | 应用 LLM 决策并继续执行 |
| `finalize(route_id: str) -> Dict[str, Any]` | 生成最终报告和 RetroGraph |

### DecisionManager (`tools/workflow/decision_manager.py`)

管理 5 类决策点的构建、应用和回退:

| 方法 | 说明 |
|------|------|
| `check_decision_point(task, result) -> Optional[DecisionContext]` | 检查是否需要决策 |
| `check_recursion_decision_point(task) -> Optional[DecisionContext]` | 检查递归决策点 |
| `apply_decision(instruction: DecisionInstruction) -> None` | 应用 LLM 决策 |
| `get_fallback_decision(decision_type, task) -> DecisionInstruction` | 获取默认回退决策 |
| `save_pending_decision(context: DecisionContext) -> None` | 持久化待决策上下文 |
| `load_pending_decision() -> Optional[Dict]` | 加载待决策上下文 |

### RouteManager (`tools/workflow/route_manager.py`)

路线创建、对比和持久化:

| 方法 | 说明 |
|------|------|
| `create_route(target_smiles) -> RetroRoute` | 创建新路线 |
| `save_route(route) -> None` | 持久化路线 |
| `load_route(route_id) -> RetroRoute` | 加载路线 |
| `compare_routes() -> Dict` | 对比所有路线 |

### ExplorationSession (`tools/workflow/exploration.py`)

管理决策点处的探索调用:

| 方法/属性 | 说明 |
|----------|------|
| `remaining -> int` | 剩余探索预算 |
| `calls -> List[Dict]` | 已执行的探索调用记录 |
| `record_call(tool_id, params, result) -> None` | 记录一次探索调用 |
| `build_exploration_log() -> List[str]` | 构建探索日志 (写入 DecisionInstruction) |

### 其他工作流模块

| 模块 | 说明 |
|------|------|
| `tools/workflow/skill_dispatch.py` | `dispatch_skill()`: 根据 TaskType 分发到对应 Skill |
| `tools/workflow/subtask_manager.py` | `create_validate_subtask()`, `expand_precursor_subtasks()`: 子任务管理 |
| `tools/workflow/snapshot.py` | `save_snapshot()`, `load_snapshot()`: 快照持久化 |
| `tools/workflow/task_state.py` | `transition_status()`: 任务状态转换规则 |

---

## 批量执行说明

使用 `batch` 子命令可以从 JSON 文件批量执行多个任务:

```bash
python scripts/run_retro.py batch tasks.json
```

JSON 任务文件格式示例:

```json
[
  {
    "target_smiles": "CC(=O)Oc1ccccc1C(=O)O",
    "output_dir": "outputs/aspirin"
  },
  {
    "target_smiles": "CC(=O)Nc1ccc(O)cc1",
    "output_dir": "outputs/paracetamol"
  }
]
```

批量执行会依次为每个目标分子创建计划并执行，遇到决策点时暂停等待输入。

---

## 测试命令

```bash
# 运行全部测试
pytest Multi-retro/tests/ -v

# 运行特定模块测试
pytest Multi-retro/tests/chem/ -v          # 化学模块
pytest Multi-retro/tests/guidance/ -v      # 指导模块
pytest Multi-retro/tests/output/ -v        # 输出模块
pytest Multi-retro/tests/skills/ -v        # Skill 测试
pytest Multi-retro/tests/workflow/ -v      # 工作流测试
pytest Multi-retro/tests/docs/ -v          # 文档属性测试

# 运行单个测试文件
pytest Multi-retro/tests/test_cli.py -v    # CLI 测试
pytest Multi-retro/tests/test_errors.py -v # 错误处理测试
```

---

## 文档索引

| 文档 | 说明 |
|------|------|
| [README.md](README.md) | 项目入口 (本文件) |
| [SKILL.md](SKILL.md) | Skill API 完整参考 (LLM 面向) |
| [scripts/README.md](scripts/README.md) | CLI 入口说明 |
| [tests/README.md](tests/README.md) | 测试套件说明 |
| [outputs/README.md](outputs/README.md) | 输出目录说明 |
| [schemas/README.md](schemas/README.md) | JSON Schema 说明 |
| [skills/README.md](skills/README.md) | Skill 文档索引 |
| [skills/analyze_molecule/SKILL.md](skills/analyze_molecule/SKILL.md) | analyze_molecule 接口规范 |
| [skills/analyze_scaffold/SKILL.md](skills/analyze_scaffold/SKILL.md) | analyze_scaffold 接口规范 |
| [skills/analyze_selectivity/SKILL.md](skills/analyze_selectivity/SKILL.md) | analyze_selectivity 接口规范 |
| [skills/check_precursor_availability/SKILL.md](skills/check_precursor_availability/SKILL.md) | check_availability 接口规范 |
| [skills/decide_strategy/SKILL.md](skills/decide_strategy/SKILL.md) | decide_strategy 接口规范 |
| [skills/map_reaction_atoms/SKILL.md](skills/map_reaction_atoms/SKILL.md) | map_atoms 接口规范 |
| [skills/multi_turn_render/SKILL.md](skills/multi_turn_render/SKILL.md) | multi_turn_render 接口规范 |
| [skills/plan_ring_strategy/SKILL.md](skills/plan_ring_strategy/SKILL.md) | plan_ring_strategy 接口规范 |
| [skills/propose_disconnection/SKILL.md](skills/propose_disconnection/SKILL.md) | propose_disconnection 接口规范 |
| [skills/render_synthesis_report/SKILL.md](skills/render_synthesis_report/SKILL.md) | render_report 接口规范 |
| [skills/repair_reaction/SKILL.md](skills/repair_reaction/SKILL.md) | repair_reaction 接口规范 |
| [skills/validate_reaction/SKILL.md](skills/validate_reaction/SKILL.md) | validate_reaction 接口规范 |

---

## 环境要求

- Python ≥ 3.10
- 依赖包见 [requirements.txt](requirements.txt)
- 核心依赖: RDKit, rxnmapper, hypothesis (测试)
