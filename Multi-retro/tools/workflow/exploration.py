"""
ExplorationManager — enables LLM autonomous exploration at decision points.

When the LLM receives a DecisionContext (exit code 10), it can call
analysis tools *before* submitting a final DecisionInstruction. This
module provides:

1. ExplorationTool definitions — what the LLM can call at each decision type
2. execute_exploration() — runs a tool and returns results without
   mutating task state
3. ExplorationSession — tracks exploration calls per decision point,
   enforces budget limits, and accumulates an exploration_log

Migrated from ``core/taskflow/exploration.py``.

Design principles:
- Read-only: explorations never change task status or route state
- Budget-limited: max N exploration calls per decision point
  (default from ``RetroLimits.MAX_EXPLORATION_CALLS``)
- Logged: every exploration call is recorded and forwarded to the
  DecisionInstruction.exploration_log for journaling
"""

from __future__ import annotations

import json
import logging
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional

from tools.common.constants import RetroLimits

logger = logging.getLogger(__name__)


# ─────────────────────────────────────────────────────────────────────────────
# Exploration Tool Definitions
# ─────────────────────────────────────────────────────────────────────────────

EXPLORATION_TOOLS: List[Dict[str, Any]] = [
    {
        "tool_id": "assess_feasibility",
        "description": "正向合成可行性评估 — 检查官能团兼容性、副反应风险、保护基需求、机理合理性",
        "required_params": ["reaction_smiles"],
        "optional_params": ["reaction_class"],
        "applicable_decisions": [
            "disconnection_decision",
            "validation_judgment",
            "repair_judgment",
        ],
        "limitations": (
            "基于静态规则 (官能团矩阵 + 模式匹配)，非反应机理动态模拟。"
            "置信度上限 0.85。LLM 需独立验证关键决策。"
        ),
    },
    {
        "tool_id": "analyze_selectivity",
        "description": "化学选择性分析 — 识别竞争官能团、反应性冲突、保护策略",
        "required_params": ["smiles"],
        "optional_params": [],
        "applicable_decisions": [
            "strategy_selection",
            "disconnection_decision",
            "validation_judgment",
            "recursion_decision",
        ],
        "limitations": (
            "基于 SMARTS 模式匹配，不考虑具体构象和溶剂效应。"
            "冲突对列表有限 (12 种)。LLM 需评估具体化学环境。"
        ),
    },
    {
        "tool_id": "analyze_molecule",
        "description": "分子分析 — 获取 atom_bond_map、官能团、SA_Score",
        "required_params": ["smiles"],
        "optional_params": [],
        "applicable_decisions": [
            "disconnection_decision",
            "recursion_decision",
        ],
        "limitations": "结构分析工具，输出客观数据，局限性较小。",
    },
    {
        "tool_id": "validate_reaction",
        "description": "反应验证 — 原子平衡 + 化学合理性检查 + 正向可行性",
        "required_params": ["reaction_smiles"],
        "optional_params": ["reaction_class"],
        "applicable_decisions": [
            "disconnection_decision",
            "repair_judgment",
        ],
        "limitations": (
            "原子平衡检查是精确的；化学合理性判断基于规则。"
            "副产物推断覆盖 18 种模式，可能遗漏非常规副产物。"
        ),
    },
    {
        "tool_id": "check_availability",
        "description": "前体可用性检查 — 查询前体是否可购买",
        "required_params": ["smiles_list"],
        "optional_params": [],
        "applicable_decisions": [
            "recursion_decision",
            "validation_judgment",
        ],
        "limitations": "基于 SA_Score 估算，非实时数据库查询。",
    },
    {
        "tool_id": "recommend_protecting_groups",
        "description": "保护基推荐 — 基于当前活跃 PG 和官能团冲突推荐正交保护基",
        "required_params": ["functional_groups"],
        "optional_params": [],
        "applicable_decisions": [
            "disconnection_decision",
            "validation_judgment",
            "repair_judgment",
        ],
        "limitations": (
            "返回建议性推荐，不会自动在路线中插入保护/脱保护步骤。"
            "LLM 需自行评估: (1) 是否真的需要保护; "
            "(2) 保护基与后续所有步骤的兼容性; "
            "(3) 脱保护时机和条件。"
        ),
    },
]


def get_tools_for_decision(decision_type: str) -> List[Dict[str, Any]]:
    """Return exploration tools applicable to the given decision type."""
    return [
        t for t in EXPLORATION_TOOLS
        if decision_type in t["applicable_decisions"]
    ]


# ─────────────────────────────────────────────────────────────────────────────
# Exploration Execution (read-only tool dispatch)
# ─────────────────────────────────────────────────────────────────────────────

def execute_exploration(
    tool_id: str,
    params: Dict[str, Any],
) -> Dict[str, Any]:
    """Execute an exploration tool and return the result.

    This is a read-only operation — no task state is mutated.

    Args:
        tool_id: One of the EXPLORATION_TOOLS tool_ids.
        params: Parameters required by the tool.

    Returns:
        Tool result dict, or ``{"success": False, "error": "..."}`` on failure.
    """
    dispatchers = {
        "assess_feasibility": _explore_feasibility,
        "analyze_selectivity": _explore_selectivity,
        "analyze_molecule": _explore_molecule,
        "validate_reaction": _explore_validate,
        "check_availability": _explore_availability,
        "recommend_protecting_groups": _explore_pg_recommend,
    }

    dispatcher = dispatchers.get(tool_id)
    if dispatcher is None:
        return {"success": False, "error": f"Unknown exploration tool: {tool_id}"}

    try:
        return dispatcher(params)
    except Exception as exc:
        logger.error("Exploration %s failed: %s", tool_id, exc)
        return {"success": False, "error": str(exc)}


# -- Individual dispatchers (lazy imports to avoid circular deps) -----------

def _explore_feasibility(params: Dict[str, Any]) -> Dict[str, Any]:
    from tools.guidance.feasibility_assessor import assess_forward_feasibility
    return assess_forward_feasibility(
        reaction_smiles=params["reaction_smiles"],
        reaction_class=params.get("reaction_class", ""),
    )


def _explore_selectivity(params: Dict[str, Any]) -> Dict[str, Any]:
    from tools.guidance.selectivity_analyzer import analyze_selectivity
    return analyze_selectivity(params["smiles"])


def _explore_molecule(params: Dict[str, Any]) -> Dict[str, Any]:
    from tools.chem.structure_analyzer import analyze_molecule
    result = analyze_molecule(params["smiles"])
    # AnalysisResult has a to_dict() method
    if hasattr(result, "to_dict"):
        return result.to_dict()
    return result


def _explore_validate(params: Dict[str, Any]) -> Dict[str, Any]:
    from tools.chem.reaction_validator import validate_reaction
    result = validate_reaction(
        reaction_smiles=params["reaction_smiles"],
        reaction_class=params.get("reaction_class", ""),
    )
    if hasattr(result, "to_dict"):
        return result.to_dict()
    return result


def _explore_availability(params: Dict[str, Any]) -> Dict[str, Any]:
    from tools.chem.sa_scorer import classify_sa
    smiles_list = params["smiles_list"]
    if isinstance(smiles_list, str):
        smiles_list = json.loads(smiles_list)
    results = []
    for smi in smiles_list:
        results.append(classify_sa(smi))
    return {"success": True, "availability": results}


def _explore_pg_recommend(params: Dict[str, Any]) -> Dict[str, Any]:
    from tools.legacy.planner.global_strategy_manager import create_strategy_manager
    manager = create_strategy_manager()
    fg_list = params["functional_groups"]
    if isinstance(fg_list, str):
        fg_list = [fg_list]
    recommendations = manager.get_pg_recommendations(fg_list)
    return {"success": True, "recommendations": recommendations}


# ─────────────────────────────────────────────────────────────────────────────
# Exploration Session (budget tracking + log accumulation)
# ─────────────────────────────────────────────────────────────────────────────

class ExplorationSession:
    """Tracks exploration calls for a single decision point.

    Persisted to ``exploration_session.json`` alongside
    ``pending_decision.json`` so the LLM can make multiple ``explore``
    CLI calls before ``decide``.
    """

    def __init__(
        self,
        route_dir: Path,
        max_calls: int = RetroLimits.MAX_EXPLORATION_CALLS,
    ) -> None:
        self.route_dir = route_dir
        self.max_calls = max_calls
        self._path = route_dir / "exploration_session.json"
        self._calls: List[Dict[str, Any]] = []
        self._load()

    def _load(self) -> None:
        if self._path.exists():
            try:
                data = json.loads(self._path.read_text(encoding="utf-8"))
                self._calls = data.get("calls", [])
            except Exception:
                self._calls = []

    def _save(self) -> None:
        self._path.parent.mkdir(parents=True, exist_ok=True)
        self._path.write_text(
            json.dumps(
                {"calls": self._calls, "max_calls": self.max_calls},
                ensure_ascii=False, indent=2, default=str,
            ),
            encoding="utf-8",
        )

    @property
    def remaining(self) -> int:
        """Number of exploration calls remaining in the budget."""
        return max(0, self.max_calls - len(self._calls))

    @property
    def calls(self) -> List[Dict[str, Any]]:
        """Copy of all recorded exploration calls."""
        return list(self._calls)

    def record_call(
        self,
        tool_id: str,
        params: Dict[str, Any],
        result: Dict[str, Any],
    ) -> None:
        """Record an exploration call."""
        self._calls.append({
            "tool_id": tool_id,
            "params": params,
            "result_summary": _summarize_result(result),
            "timestamp": datetime.now().isoformat(),
        })
        self._save()

    def build_exploration_log(self) -> List[str]:
        """Build exploration_log entries for DecisionInstruction."""
        return [
            f"[{call['tool_id']}] {call['result_summary']}"
            for call in self._calls
        ]

    def clear(self) -> None:
        """Clear session (called after decide)."""
        self._calls = []
        try:
            if self._path.exists():
                self._path.unlink()
        except Exception:
            pass


# ─────────────────────────────────────────────────────────────────────────────
# Result summarisation (module-level helper)
# ─────────────────────────────────────────────────────────────────────────────

def _summarize_result(result: Dict[str, Any]) -> str:
    """Produce a one-line summary of a tool result."""
    if not result.get("success", True):
        return f"失败: {result.get('error', 'unknown')}"

    # Feasibility
    if "feasibility_score" in result:
        score = result["feasibility_score"]
        conf = result.get("overall_confidence", "?")
        concerns = len(result.get("concerns", []))
        pg = len(result.get("protecting_group_requirements", []))
        mech_count = result.get("mechanistic_analysis", {}).get("mechanism_count", 0)
        cond_conflicts = len(result.get("condition_conflicts", []))
        parts = [f"可行性={score}", f"置信度={conf}"]
        if concerns:
            parts.append(f"风险点={concerns}")
        if pg:
            parts.append(f"需保护基={pg}")
        if mech_count:
            parts.append(f"匹配机理={mech_count}")
        if cond_conflicts:
            parts.append(f"条件冲突={cond_conflicts}")
        parts.append("[静态规则分析,需LLM独立验证]")
        return ", ".join(parts)

    # Selectivity
    report = result.get("report", {})
    if "competing_groups" in report:
        competing = len(report["competing_groups"])
        conflicts = len(report.get("reactivity_conflicts", []))
        complexity = report.get("overall_complexity", "?")
        return (
            f"竞争基团={competing}, 冲突={conflicts}, "
            f"复杂度={complexity} [模式匹配分析,需LLM验证具体环境]"
        )

    # Validation
    if "is_valid" in result:
        valid = result["is_valid"]
        verdict = result.get("verdict", "?")
        return f"验证={verdict}, is_valid={valid}"

    # Molecule analysis
    summary = result.get("summary", {})
    if "sa_score" in summary:
        sa = summary["sa_score"]
        mw = summary.get("molecular_weight", "?")
        return f"SA={sa}, MW={mw}"

    # Availability
    if "availability" in result:
        avail = result["availability"]
        if isinstance(avail, list):
            available = sum(1 for a in avail if a.get("available"))
            return f"可购买={available}/{len(avail)}"

    # PG recommendations
    if "recommendations" in result:
        recs = result["recommendations"]
        return f"保护基建议={len(recs)}组 [建议性,需LLM评估全局兼容性]"

    # Fallback
    keys = list(result.keys())[:4]
    return f"keys={keys}"
