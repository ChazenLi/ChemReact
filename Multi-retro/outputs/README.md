# outputs/ - 输出目录

## 概述

`outputs/` 是默认的逆合成分析输出目录。运行 `run_retro.py` 后，结果会输出到此目录下的子文件夹中。

## 目录结构

```
outputs/
├── retro/                         # 默认会话目录
│   ├── routes/                    # 路线数据
│   │   ├── route_001.json         # 单条路线数据
│   │   └── route_002.json
│   ├── snapshots/                 # 快照文件 (用于暂停-恢复)
│   │   ├── route_001_pending.json # 待决策快照
│   │   └── route_001_final.json   # 最终快照
│   ├── retro_graph.json           # 完整 RetroGraph (DAG 格式)
│   └── SYNTHESIS_REPORT.md        # 人类可读合成报告
│
└── README.md                      # 本文档
```

## 会话目录说明

每个逆合成会话 (由 `plan` 子命令创建) 包含:

| 文件/目录 | 说明 |
|-----------|------|
| `routes/` | 存储所有 RetroRoute 数据 |
| `routes/route_XXX.json` | 单条路线的完整数据 (任务、元数据、状态) |
| `snapshots/` | 决策点快照，用于暂停-恢复 |
| `snapshots/route_XXX_pending.json` | 待决策上下文 (DecisionContext) |
| `retro_graph.json` | 最终的 RetroGraph (节点 + 边) |
| `SYNTHESIS_REPORT.md` | Markdown 格式的合成报告 |

## 使用 CLI 操作

```bash
# 创建规划 (生成会话目录)
python scripts/run_retro.py plan --target_smiles "CC(=O)Oc1ccccc1C(=O)O" --output_dir outputs/retro

# 执行任务
python scripts/run_retro.py run --session_dir outputs/retro --route route_001

# 查看状态
python scripts/run_retro.py status --session_dir outputs/retro

# 生成最终报告
python scripts/run_retro.py finalize --session_dir outputs/retro --route route_001
```

## RetroGraph 结构

`retro_graph.json` 包含:

```json
{
  "target_smiles": "CC(=O)Oc1ccccc1C(=O)O",
  "nodes": {
    "node_001": {
      "node_id": "node_001",
      "smiles": "CC(=O)Oc1ccccc1C(=O)O",
      "node_type": "target",
      "status": "completed",
      "sa_score": 1.85,
      "is_terminal": false
    }
  },
  "edges": {
    "edge_001": {
      "edge_id": "edge_001",
      "product_id": "node_001",
      "precursor_ids": ["node_002", "node_003"],
      "reaction_smiles": "Oc1ccccc1C(=O)O.CC(=O)Cl>>CC(=O)Oc1ccccc1C(=O)O",
      "reaction_type": "esterification"
    }
  },
  "metadata": {
    "route_id": "route_001",
    "created_at": "2025-01-15T10:30:00Z"
  }
}
```

## 注意事项

- 此目录内容为生成文件，不纳入版本控制
- 可通过 `--output_dir` 参数指定其他输出位置
- 快照文件用于实现 LLM 驱动的暂停-恢复决策协议
