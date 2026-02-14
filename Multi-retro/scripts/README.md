# scripts/ — 统一 CLI 入口

本目录包含 Multi-retro 项目的唯一命令行入口脚本 `run_retro.py`，通过 9 个子命令覆盖逆合成规划、执行、决策、探索、对比、报告生成和批量运行等全部功能。

> **注意**: 旧版本中的 `run_mesh.py`、`run_skill.py`、`run_taskflow.py`、`retro_batch.py`、`core_pipeline.py` 已全部删除，其功能已合并到 `run_retro.py` 的子命令中。

## 文件

| 文件 | 说明 |
|------|------|
| `run_retro.py` | Multi-retro 统一 CLI 入口，提供 9 个子命令 |

## 子命令总览

```
python scripts/run_retro.py <subcommand> [options]
```

| 子命令 | 说明 |
|--------|------|
| `skill` | 按名称调用单个 Skill (传 `list` 查看可用 Skill) |
| `plan` | 创建逆合成规划 |
| `run` | 执行会话中的任务 |
| `decide` | 为暂停的会话提交决策 |
| `explore` | 在决策点调用探索工具 |
| `compare` | 多路线对比 |
| `finalize` | 最终化并生成报告 |
| `status` | 查看会话状态 |
| `batch` | 从 JSON 任务文件批量执行 |

---

## 子命令详细说明

### `skill` — 调用单个 Skill

按名称调用任意已注册的 Skill，传 `list` 可列出所有可用 Skill。

```bash
python scripts/run_retro.py skill <skill_name> [--args '<JSON>']
```

| 参数 | 类型 | 必需 | 说明 |
|------|------|------|------|
| `skill_name` | str (位置参数) | ✅ | Skill 名称，或 `list` 查看列表 |
| `--args` | str | ❌ | Skill 参数的 JSON 字符串 |

示例:

```bash
# 列出所有可用 Skill
python scripts/run_retro.py skill list

# 分析分子
python scripts/run_retro.py skill analyze_molecule --args '{"smiles": "CCO"}'

# 断键
python scripts/run_retro.py skill break_bond --args '{"smiles": "CCO", "atom1_idx": 0, "atom2_idx": 1}'
```

### `plan` — 创建逆合成规划

为目标分子创建逆合成规划，生成会话目录和初始路线。

```bash
python scripts/run_retro.py plan --target_smiles "<SMILES>" [--output_dir <DIR>]
```

| 参数 | 类型 | 必需 | 说明 |
|------|------|------|------|
| `--target_smiles` | str | ✅ | 目标分子 SMILES |
| `--output_dir` | str | ❌ | 输出目录 (默认: `outputs/retro`) |

示例:

```bash
python scripts/run_retro.py plan --target_smiles "CC(=O)Oc1ccccc1C(=O)O" --output_dir outputs/aspirin
```

### `run` — 执行任务

执行会话中的任务，遇到决策点时暂停。

```bash
python scripts/run_retro.py run --session_dir <DIR> --route <ROUTE_ID>
```

| 参数 | 类型 | 必需 | 说明 |
|------|------|------|------|
| `--session_dir` | str | ✅ | 会话目录路径 |
| `--route` | str | ✅ | 路线 ID |

示例:

```bash
python scripts/run_retro.py run --session_dir outputs/aspirin --route route_001
```

### `decide` — 提交决策

为暂停的会话提交 LLM 决策，恢复执行。

```bash
python scripts/run_retro.py decide --session_dir <DIR> --route <ROUTE_ID> --decision '<JSON>'
```

| 参数 | 类型 | 必需 | 说明 |
|------|------|------|------|
| `--session_dir` | str | ✅ | 会话目录路径 |
| `--route` | str | ✅ | 路线 ID |
| `--decision` | str | ✅ | 决策 JSON 字符串，或 `@filepath` 从文件加载 |

示例:

```bash
python scripts/run_retro.py decide --session_dir outputs/aspirin --route route_001 \
  --decision '{"task_id": "task_003", "action": "accept", "params": {}}'
```

### `explore` — 探索工具调用

在决策点调用探索工具获取更多信息，辅助决策。

```bash
python scripts/run_retro.py explore --session_dir <DIR> --route <ROUTE_ID> --tool <TOOL_ID> [--params '<JSON>']
```

| 参数 | 类型 | 必需 | 说明 |
|------|------|------|------|
| `--session_dir` | str | ✅ | 会话目录路径 |
| `--route` | str | ✅ | 路线 ID |
| `--tool` | str | ✅ | 探索工具 ID |
| `--params` | str | ❌ | 工具参数 JSON 字符串 |

可用探索工具: `assess_feasibility`, `analyze_selectivity`, `analyze_molecule`, `validate_reaction`, `check_availability`, `recommend_protecting_groups`

示例:

```bash
python scripts/run_retro.py explore --session_dir outputs/aspirin --route route_001 \
  --tool assess_feasibility --params '{"reaction_smiles": "A>>B"}'
```

### `compare` — 多路线对比

对比同一会话中的多条合成路线。

```bash
python scripts/run_retro.py compare --session_dir <DIR>
```

| 参数 | 类型 | 必需 | 说明 |
|------|------|------|------|
| `--session_dir` | str | ✅ | 会话目录路径 |

示例:

```bash
python scripts/run_retro.py compare --session_dir outputs/aspirin
```

### `finalize` — 生成最终报告

为完成的路线生成 RetroGraph 和合成报告。

```bash
python scripts/run_retro.py finalize --session_dir <DIR> --route <ROUTE_ID>
```

| 参数 | 类型 | 必需 | 说明 |
|------|------|------|------|
| `--session_dir` | str | ✅ | 会话目录路径 |
| `--route` | str | ✅ | 路线 ID |

示例:

```bash
python scripts/run_retro.py finalize --session_dir outputs/aspirin --route route_001
```

### `status` — 查看会话状态

查看会话的当前状态、任务进度和待决策项。

```bash
python scripts/run_retro.py status --session_dir <DIR>
```

| 参数 | 类型 | 必需 | 说明 |
|------|------|------|------|
| `--session_dir` | str | ✅ | 会话目录路径 |

示例:

```bash
python scripts/run_retro.py status --session_dir outputs/aspirin
```

### `batch` — 批量执行

从 JSON 文件批量执行多个逆合成任务。

```bash
python scripts/run_retro.py batch <task_file>
```

| 参数 | 类型 | 必需 | 说明 |
|------|------|------|------|
| `task_file` | str (位置参数) | ✅ | JSON 任务文件路径 |

任务文件格式:

```json
[
  {"target_smiles": "CC(=O)Oc1ccccc1C(=O)O", "output_dir": "outputs/aspirin"},
  {"target_smiles": "CC(=O)Nc1ccc(O)cc1", "output_dir": "outputs/paracetamol"}
]
```

示例:

```bash
python scripts/run_retro.py batch tasks.json
```

---

## 完整工作流示例

```bash
# 1. 创建规划
python scripts/run_retro.py plan --target_smiles "CC(=O)Oc1ccccc1C(=O)O" --output_dir outputs/aspirin

# 2. 执行任务
python scripts/run_retro.py run --session_dir outputs/aspirin --route route_001

# 3. (可选) 探索
python scripts/run_retro.py explore --session_dir outputs/aspirin --route route_001 \
  --tool analyze_molecule --params '{"smiles": "Oc1ccccc1C(=O)O"}'

# 4. 提交决策
python scripts/run_retro.py decide --session_dir outputs/aspirin --route route_001 \
  --decision '{"task_id": "1", "action": "select_bond", "params": {"atom1_idx": 3, "atom2_idx": 4}}'

# 5. 重复步骤 2-4 直到完成

# 6. 查看状态
python scripts/run_retro.py status --session_dir outputs/aspirin

# 7. 生成报告
python scripts/run_retro.py finalize --session_dir outputs/aspirin --route route_001
```

---

## 参考

完整 API 参考见 [../SKILL.md](../SKILL.md)
