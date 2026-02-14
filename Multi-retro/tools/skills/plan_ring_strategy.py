"""
Plan Ring Strategy Skill
========================

Thin skill wrapper around :func:`tools.guidance.strategy_advisor.plan_ring_strategy`.
Accepts a SMILES string, delegates to the tool module, and returns
a :class:`SkillResult` envelope.
"""

from __future__ import annotations

import logging
from typing import Any, Dict

from tools.skills.base import BaseSkill, SkillResult

logger = logging.getLogger(__name__)


class PlanRingStrategySkill(BaseSkill):
    """Plan ring formation strategy including construction order.

    Delegates to :func:`tools.guidance.strategy_advisor.plan_ring_strategy`.
    """

    name = "plan_ring_strategy"
    description = "Plan ring formation strategy including construction order and retro-ring-opening"

    def execute(self, args: Any) -> Dict[str, Any]:
        """Run ring strategy planning and return a SkillResult dict.

        Args:
            args: dict with ``smiles`` key, or any object with a ``.smiles`` attribute.

        Returns:
            ``SkillResult.to_dict()`` with ring strategy nested under ``data``.
        """
        smiles = _extract_smiles(args)
        if smiles is None:
            return SkillResult(
                success=False, error="args must contain 'smiles'"
            ).to_dict()

        try:
            from tools.guidance.strategy_advisor import plan_ring_strategy

            result = plan_ring_strategy(smiles)
        except Exception as exc:
            logger.exception("plan_ring_strategy failed")
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
