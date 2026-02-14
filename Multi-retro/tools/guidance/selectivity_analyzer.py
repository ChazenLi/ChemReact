"""
Selectivity Analyzer Module
============================
Chemoselectivity analysis: competing functional groups, protection requirements.
Extracted from core/skills/strategy.py (AnalyzeSelectivitySkill).

Public API:
    analyze_selectivity(smiles) -> Dict[str, Any]
"""

from __future__ import annotations

import logging
from typing import Any, Dict

logger = logging.getLogger(__name__)


def analyze_selectivity(smiles: str) -> Dict[str, Any]:
    """Analyze chemoselectivity challenges for a molecule.

    Returns:
        {
            "success": bool,
            "report": {
                "competing_groups": [...],
                "overall_complexity": str,
                "protection_requirements": [...],
                ...
            }
        }
    """
    try:
        from tools.legacy.analysis.chemoselectivity import analyze_chemoselectivity
        return analyze_chemoselectivity(smiles)
    except Exception as e:
        return {"success": False, "error": str(e)}
