# Cobimetinib 逆合成问题分析与优化计划

**目标分子**: Cobimetinib (MEK 抑制剂)
**SMILES**: `O=C(c1ccc(F)c(F)c1Nc1ccc(I)cc1F)N1CC(O)([C@@H]2CCCCN2)C1`
**分析日期**: 2026-02-11

---

## 一、逆合成过程中遇到的问题

### 问题 P1: 二芳胺键 (Ar-NH-Ar) 缺乏专用语义标注

**现象**: Cobimetinib 核心骨架包含一个二芳胺键 (Ar-NH-Ar)，即两个芳环通过 -NH- 连接。当前 `_annotate_bond()` 将其归类为通用的 "芳基-N 键"，与普通的 Ar-NHR (如苄胺) 混为一谈。

**影响**:
- 二芳胺的逆合成策略与普通芳基胺不同: 二芳胺通常通过 Buchwald-Hartwig 偶联或 SNAr 形成，而非还原胺化
- 当前 `potential_reactions` 列表包含 "还原胺化"，对二芳胺来说不合理
- LLM 可能选择错误的逆合成策略

**根因**: `_annotate_bond()` 中 C-N 键的分支只检查了 "是否有一端是芳香的"，没有进一步区分 "两端都连接芳环" 的情况。

**严重度**: P1 (影响逆合成策略选择)

---

### 问题 P2: break_bond 对酰胺键断裂生成醛而非酰氯

**现象**: 断裂 Cobimetinib 的酰胺键 (C(=O)-N, atom 1-19) 时，`_break_bond_and_get_fragments()` 基于拓扑断键 + 加H补价，生成的片段 A 是醛 (`O=Cc1ccc(F)c(F)c1Nc1ccc(I)cc1F`) 而非酰氯 (`O=C(Cl)c1ccc(F)c(F)c1Nc1ccc(I)cc1F`)。

**影响**:
- 需要人工/LLM 手动将醛修正为酰氯，才能构建正确的反应 SMILES
- 在 HDAC 逆合成中已为磺酰胺和 Heck 反应添加了反应感知后处理 (`_suggest_sulfonamide_retro`, `_suggest_heck_retro`)，但酰胺键缺少类似的自动修正

**根因**: `BreakBondSkill.execute()` 的反应感知后处理层只覆盖了磺酰胺和 Heck 两种反应类型，没有覆盖酰胺键断裂 (应自动建议 酰氯 + 胺 或 羧酸 + 胺)。

**严重度**: P0 (直接影响前体正确性)

---

### 问题 P3: break_bond 对子片段的二次断键需要手动传入子片段 SMILES

**现象**: 第一次断裂酰胺键后得到片段 A (二芳胺苯甲酰氯中间体)，需要对片段 A 再次调用 `break_bond` 断裂 Ar-NH-Ar 键。但 `break_bond` 的原子索引是基于输入 SMILES 的 canonical 编号，片段 A 的原子索引与原始分子完全不同。

**影响**:
- 必须先对片段 A 单独调用 `analyze_molecule` 获取新的 atom_bond_map
- 然后在新的索引体系中找到目标键
- 整个流程需要 3 次工具调用 (analyze -> 找索引 -> break_bond)，而非 1 次

**根因**: 这是拓扑断键的固有限制 — 断键后生成新分子，原子索引重新编号。没有提供 "级联断键" 或 "基于化学键类型的批量断键" 功能。

**严重度**: P1 (增加交互步骤，但不影响正确性)

---

### 问题 P4: 报告中缺少反应条件、副产物、合成注释

**现象**: `retro_graph.json` 中包含丰富的信息 (`conditions`, `byproducts`, `notes`)，但 `_generate_graph_report()` 只渲染了 `reaction_smiles` 和 `reaction_type`，丢弃了条件、副产物、合成注释等关键信息。

**影响**:
- 生成的 SYNTHESIS_REPORT.md 信息不完整，实验化学家无法直接使用
- `synthesis_summary` 中的分步说明也未渲染

**根因**: `_generate_graph_report()` 是早期版本的简单实现，只处理了最基本的字段。

**严重度**: P1 (报告质量)

---

### 问题 P5: 报告中反应状态显示 "pending" 而非 "validated"

**现象**: 两个反应 (E1 酰胺偶联, E2 SNAr) 在 `retro_graph.json` 中有 `"validated": true`，但报告中显示 `Status: pending`。

**根因**: `_generate_graph_report()` 读取 `edata.get("status", "pending")`，但 retro_graph.json 中边没有 `status` 字段，只有 `validated` 布尔字段。报告生成器不识别 `validated` 字段。

**严重度**: P2 (显示问题)

---

### 问题 P6: Windows 路径分隔符问题

**现象**: 报告中的图片路径使用反斜杠 (`outputs/cobimetinib\images\T.png`)，在 Markdown 渲染器中可能无法正确显示。

**根因**: `visualizer.py` 使用 `os.path.join()` 生成路径，在 Windows 上产生反斜杠。

**严重度**: P2 (跨平台兼容性)

---

### 问题 P7: 选择性酰化的保护基策略未自动提示

**现象**: Cobimetinib 的最后一步是酰氯与氮杂环丁烷-哌啶双胺 (P3) 的酰胺偶联。P3 同时含有氮杂环丁烷 NH 和哌啶 NH，需要选择性酰化。`retro_graph.json` 的 `notes` 中提到 "May need Boc protection on piperidine N"，但这个信息完全依赖人工判断。

**影响**:
- `analyze_selectivity` skill 可以分析化学选择性，但在交互模式中没有被自动调用
- LLM 可能忽略保护基需求

**根因**: 交互模式的 8 步流程中没有强制的选择性分析步骤。`analyze_selectivity` 只在自动化 pipeline 中被调用。

**严重度**: P1 (影响合成可行性)

---

## 二、与 HDAC 逆合成对比 — 已解决 vs 新发现

| 问题 | HDAC 中发现 | Cobimetinib 中状态 |
|------|------------|-------------------|
| GBK 编码崩溃 | 已修复 (v3.6) | ✅ 无问题 |
| 磺酰胺键标注 | 已修复 (v3.6) | ✅ 不涉及 |
| S=O 双键误标为可断 | 已修复 (v3.6) | ✅ 无问题 |
| break_bond 磺酰胺修正 | 已修复 (v3.6) | ✅ 不涉及 |
| 逗号分隔 SMILES 支持 | 已修复 (v3.6) | ✅ 无问题 |
| **二芳胺键标注** | 未涉及 | ❌ **新问题 P1** |
| **酰胺键 break_bond 修正** | 未涉及 | ❌ **新问题 P2** |
| **级联断键** | 未涉及 | ❌ **新问题 P3** |
| **报告缺少条件/副产物** | 未涉及 | ❌ **新问题 P4** |
| **validated 状态显示** | 未涉及 | ❌ **新问题 P5** |
| **Windows 路径分隔符** | 未涉及 | ❌ **新问题 P6** |
| **选择性分析未集成** | 未涉及 | ❌ **新问题 P7** |

---

## 三、优化计划

### P0 级 (必须修复 — 影响正确性)

#### P0-A: 酰胺键反应感知后处理 [对应 P2]

**目标**: `BreakBondSkill.execute()` 断裂酰胺键时，自动建议 酰氯+胺 或 羧酸+胺 前体。

**实现方案**:
```
新增方法: _suggest_amide_retro(self, mol, bond, c_idx, n_idx, fragments, canonical)
逻辑:
  1. 识别酰胺键的 C 端和 N 端
  2. C 端片段: 将 C-H 修正为 C-Cl (酰氯) 或 C-OH (羧酸)
  3. N 端片段: 保持为胺 (已正确)
  4. 返回 suggested_fragments, suggested_reaction_smiles, suggested_reaction_type
```

**修改文件**: `multiretro/core/skills/builtin.py`
- `BreakBondSkill.execute()`: 在反应感知后处理层添加酰胺键分支
- 新增 `_suggest_amide_retro()` 方法

---

### P1 级 (重要优化 — 影响策略质量)

#### P1-A: 二芳胺键专用标注 [对应 P1]

**目标**: 区分 Ar-NH-Ar (二芳胺) 和 Ar-NHR (普通芳基胺)。

**实现方案**:
```
在 _annotate_bond() 的 C-N 分支中:
  if a1_aromatic or a2_aromatic:
    # 检查 N 的另一端是否也连接芳环
    n_atom = atom1 if sym1 == "N" else atom2
    other_atom = atom2 if sym1 == "N" else atom1
    n_has_other_aromatic = any(
        nb.GetIsAromatic() and nb.GetIdx() != other_atom.GetIdx()
        for nb in n_atom.GetNeighbors()
    )
    if n_has_other_aromatic:
        label = "二芳胺键 Ar-NH-Ar"
        reactions = ["Buchwald-Hartwig偶联", "Ullmann偶联", "SNAr取代"]
        hint = "断裂 -> 芳基卤化物 + 芳基胺 (Pd催化或SNAr)"
        priority = 1
    else:
        # 保持现有的芳基-N标注
```

**修改文件**: `multiretro/core/skills/builtin.py` — `_annotate_bond()`

#### P1-B: 二芳胺 break_bond 反应感知 [对应 P1 延伸]

**目标**: 断裂 Ar-NH-Ar 键时，自动建议 芳基卤化物 + 芳基胺。

**实现方案**:
```
新增方法: _suggest_diarylamine_retro(self, mol, bond, ...)
逻辑:
  1. N 端: 保持为 ArNH2 (已正确)
  2. Ar 端: 将 Ar-H 修正为 Ar-X (X = Br/I/F，取决于芳环上已有取代基)
  3. 对于 Cobimetinib 的情况: 2,4-二氟苯甲酰氯的 ortho-F 被 SNAr 取代
```

**修改文件**: `multiretro/core/skills/builtin.py`

#### P1-C: 报告渲染增强 [对应 P4]

**目标**: 报告中包含反应条件、副产物、合成注释、分步总结。

**实现方案**:
```
_generate_graph_report() 增强:
  1. 每个反应下方添加:
     - **Conditions**: {conditions}
     - **Byproducts**: {byproducts}
     - **Notes**: {notes}
  2. 末尾添加 "## Synthesis Summary" 章节
     - 渲染 synthesis_summary 中的分步说明
  3. 每个节点添加 availability 和 sa_score 信息
```

**修改文件**: `multiretro/core/graph/visualizer.py` — `_generate_graph_report()`

#### P1-D: 选择性分析集成到交互模式 [对应 P7]

**目标**: 在交互模式中，当检测到多个竞争性官能团时，自动提示调用 `analyze_selectivity`。

**实现方案**:
```
在 retro_guidance 中添加选择性警告:
  - 检测到同一前体中有多个相同类型的亲核/亲电位点时
  - 输出: "⚠ 前体 P3 含有 2 个 NH 基团，建议调用 analyze_selectivity 评估选择性"
  - 在 sampleflow.md 中添加可选的 Step 6.5: 选择性分析
```

**修改文件**:
- `multiretro/core/skills/builtin.py` — `_build_retro_guidance()`
- `multiretro/docs/sampleflow.md` (如存在)

---

### P2 级 (改善体验)

#### P2-A: validated 状态映射 [对应 P5]

**实现**: `_generate_graph_report()` 中，如果 `edata.get("validated")` 为 True 且无 `status` 字段，则显示 "validated"。

**修改文件**: `multiretro/core/graph/visualizer.py`

#### P2-B: 路径分隔符统一 [对应 P6]

**实现**: `visualizer.py` 中所有路径使用 `pathlib.PurePosixPath` 或 `.replace("\\", "/")` 确保 Markdown 兼容。

**修改文件**: `multiretro/core/graph/visualizer.py`

#### P2-C: 级联断键辅助 [对应 P3]

**实现**: `BreakBondSkill` 返回结果中，对每个片段附带 `fragment_atom_bond_map` (调用 `_build_atom_bond_map`)，省去用户手动调用 `analyze_molecule` 的步骤。

**修改文件**: `multiretro/core/skills/builtin.py` — `BreakBondSkill.execute()`

---

## 四、优先级排序与实施顺序

```
Phase 1 (P0): 酰胺键反应感知后处理
  -> 直接影响前体正确性，最高优先级

Phase 2 (P1-A + P1-B): 二芳胺键标注 + break_bond 修正
  -> 提升策略质量，与 P0 同类问题

Phase 3 (P1-C): 报告渲染增强
  -> 提升输出质量

Phase 4 (P1-D): 选择性分析集成
  -> 提升合成可行性评估

Phase 5 (P2-A + P2-B + P2-C): 体验优化
  -> 非阻塞性改善
```

---

## 五、总结

Cobimetinib 逆合成暴露了 7 个问题，其中 1 个 P0、4 个 P1、3 个 P2。与 HDAC 逆合成相比，v3.6 的编码安全和磺酰胺修复工作良好，但新的分子结构 (二芳胺、酰胺偶联、双胺选择性) 暴露了新的覆盖盲区。核心问题集中在:

1. **反应感知后处理的覆盖面不足** — 只有磺酰胺和 Heck，缺少酰胺和二芳胺
2. **报告渲染信息丢失** — retro_graph.json 中的丰富数据未被充分利用
3. **交互模式缺少选择性分析** — 多官能团前体的保护基策略依赖人工判断
