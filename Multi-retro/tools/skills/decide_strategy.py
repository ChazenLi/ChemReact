"""
Decide Strategy Skill
=====================

Thin skill wrapper around :func:`tools.guidance.strategy_advisor.decide_strategy`.
Accepts a SMILES string, delegates to the tool module, and returns
a :class:`SkillResult` envelope.
"""

from __future__ import annotations

import logging
from typing import Any, Dict

from tools.skills.base import BaseSkill, SkillResult

logger = logging.getLogger(__name__)


class DecideStrategySkill(BaseSkill):
    """Comprehensive strategy decision (linear/convergent/mixed).

    Delegates to :func:`tools.guidance.strategy_advisor.decide_strategy`.
    """

    name = "decide_strategy"
    description = "Comprehensive strategy decision (linear/convergent) based on all analyses"

    def execute(self, args: Any) -> Dict[str, Any]:
        """Run strategy decision and return a SkillResult dict.

        Args:
            args: dict with ``smiles`` key, or any object with a ``.smiles`` attribute.

        Returns:
            ``SkillResult.to_dict()`` with strategy data nested under ``data``.
        """
        smiles = _extract_smiles(args)
        if smiles is None:
            return SkillResult(
                success=False, error="args must contain 'smiles'"
            ).to_dict()

        try:
            from tools.guidance.strategy_advisor import decide_strategy

            result = decide_strategy(smiles)
        except Exception as exc:
            logger.exception("decide_strategy failed")
            return SkillResult(success=False, error=str(exc)).to_dict()

        success = result.get("success", False)
        error = result.get("error", "")
        return SkillResult(success=success, data=result, error=error).to_dict()


def _extract_smiles(args: Any) -> str | None:
    """Normalise various arg formats to a SMILES string."""
    if isinstance(args, dict):
        return args.get("smiles") or None
    if hasattr(args, "smiles"):
        return args.smiles or None
    return None
