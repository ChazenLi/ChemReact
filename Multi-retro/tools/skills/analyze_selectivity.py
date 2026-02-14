"""
Analyze Selectivity Skill
=========================

Thin skill wrapper around :func:`tools.guidance.selectivity_analyzer.analyze_selectivity`.
Accepts a SMILES string, delegates to the tool module, and returns
a :class:`SkillResult` envelope.
"""

from __future__ import annotations

import logging
from typing import Any, Dict

from tools.skills.base import BaseSkill, SkillResult

logger = logging.getLogger(__name__)


class AnalyzeSelectivitySkill(BaseSkill):
    """Analyze chemoselectivity, competing FGs, and protection requirements.

    Delegates to :func:`tools.guidance.selectivity_analyzer.analyze_selectivity`.
    """

    name = "analyze_selectivity"
    description = "Analyze chemoselectivity, competing FGs, and protection requirements"

    def execute(self, args: Any) -> Dict[str, Any]:
        """Run selectivity analysis and return a SkillResult dict.

        Args:
            args: dict with ``smiles`` key, or any object with a ``.smiles`` attribute.

        Returns:
            ``SkillResult.to_dict()`` with selectivity data nested under ``data``.
        """
        smiles = _extract_smiles(args)
        if smiles is None:
            return SkillResult(
                success=False, error="args must contain 'smiles'"
            ).to_dict()

        try:
            from tools.guidance.selectivity_analyzer import analyze_selectivity

            result = analyze_selectivity(smiles)
        except Exception as exc:
            logger.exception("analyze_selectivity failed")
            return SkillResult(success=False, error=str(exc)).to_dict()

        success = result.get("success", False) if isinstance(result, dict) else False
        error = result.get("error", "") if isinstance(result, dict) else ""
        return SkillResult(success=success, data=result, error=error).to_dict()


def _extract_smiles(args: Any) -> str | None:
    """Normalise various arg formats to a SMILES string."""
    if isinstance(args, dict):
        return args.get("smiles") or None
    if hasattr(args, "smiles"):
        return args.smiles or None
    return None
