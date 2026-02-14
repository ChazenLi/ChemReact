"""
Strategy Advisor Module
=======================
Comprehensive strategy decision combining scaffold analysis, ring strategy,
and selectivity analysis. Recommends linear, convergent, or mixed strategy.

Migrated from core/skills/strategy.py (DecideStrategySkill) and
core/skills/global_strategy.py (GetGlobalStrategySkill).

Public API:
    decide_strategy(smiles) -> Dict[str, Any]
    analyze_scaffold(smiles) -> Dict[str, Any]
    plan_ring_strategy(smiles) -> Dict[str, Any]
"""

from __future__ import annotations

import logging
from typing import Any, Dict

logger = logging.getLogger(__name__)

_STRATEGY_THRESHOLD = 0.15


def decide_strategy(smiles: str) -> Dict[str, Any]:
    """Decide synthesis strategy based on scaffold, ring, and selectivity analyses.

    Returns:
        {
            "success": bool,
            "recommended_strategy": "linear" | "convergent" | "mixed",
            "alternative_strategies": [...],
            "rationale": str,
            "scores": {"convergent_feasibility": float, "linear_feasibility": float},
            "scaffold_summary": {...},
            "ring_strategy_summary": {...},
            "selectivity_summary": {...},
        }
    """
    try:
        from tools.legacy.analysis.chemoselectivity import analyze_chemoselectivity
        from tools.legacy.analysis.scaffold_analyzer import analyze_scaffold as _analyze_scaffold
        from tools.legacy.planner.ring_strategy import plan_ring_strategy as _plan_ring

        scaffold_result = _analyze_scaffold(smiles)
        ring_result = _plan_ring(scaffold_result)
        selectivity_result = analyze_chemoselectivity(smiles)
    except Exception as e:
        return {"success": False, "error": str(e)}

    scaffold = scaffold_result.get("report", {})
    ring_plan = ring_result.get("plan", {})
    selectivity = selectivity_result.get("report", {})

    conv_score = scaffold.get("convergent_feasibility", 0.5)
    lin_score = scaffold.get("linear_feasibility", 0.5)
    n_rings = scaffold.get("ring_count", 0)
    competing = len(selectivity.get("competing_groups", []))
    complexity = selectivity.get("overall_complexity", "medium")

    if conv_score > lin_score and n_rings >= 2:
        strategy_type = "convergent"
        rationale = "Complex ring-rich scaffold favors convergent synthesis."
    elif complexity == "high" or competing >= 2:
        strategy_type = "mixed"
        rationale = "Functional-group complexity suggests mixed strategy with selectivity controls."
    else:
        strategy_type = "linear"
        rationale = "Structure is relatively simple; linear synthesis is likely efficient."

    # Build alternatives
    candidate_pool = [("convergent", conv_score), ("linear", lin_score)]
    if complexity == "high" or competing >= 2:
        candidate_pool.append(("mixed", max(conv_score, lin_score) * 0.8))
    alternative_strategies = [
        name for name, score in sorted(candidate_pool, key=lambda x: x[1], reverse=True)
        if name != strategy_type and score >= _STRATEGY_THRESHOLD
    ]

    return {
        "success": True,
        "recommended_strategy": strategy_type,
        "alternative_strategies": alternative_strategies,
        "rationale": rationale,
        "scores": {
            "convergent_feasibility": conv_score,
            "linear_feasibility": lin_score,
        },
        "scaffold_summary": {
            "type": scaffold.get("scaffold_type"),
            "ring_count": n_rings,
            "complexity": scaffold.get("complexity_score"),
        },
        "ring_strategy_summary": {
            "total_rings": ring_plan.get("total_rings"),
            "reactions_used": ring_plan.get("ring_forming_reactions_used", []),
        },
        "selectivity_summary": {
            "competing_groups": competing,
            "overall_complexity": complexity,
            "protection_needed": len(selectivity.get("protection_requirements", [])),
        },
    }


def analyze_scaffold(smiles: str) -> Dict[str, Any]:
    """Analyze molecular scaffold for strategic disconnections.

    Returns scaffold analysis result from tools.legacy.analysis.scaffold_analyzer.
    """
    try:
        from tools.legacy.analysis.scaffold_analyzer import analyze_scaffold as _analyze
        return _analyze(smiles)
    except Exception as e:
        return {"success": False, "error": str(e)}


def plan_ring_strategy(smiles: str) -> Dict[str, Any]:
    """Plan ring formation strategy including construction order.

    Returns scaffold analysis + ring strategy plan.
    """
    try:
        from tools.legacy.analysis.scaffold_analyzer import analyze_scaffold as _analyze
        from tools.legacy.planner.ring_strategy import plan_ring_strategy as _plan

        scaffold_result = _analyze(smiles)
        ring_result = _plan(scaffold_result)
        return {
            "success": True,
            "scaffold_analysis": scaffold_result.get("report", {}),
            "ring_strategy": ring_result.get("plan", {}),
        }
    except Exception as e:
        return {"success": False, "error": str(e)}
