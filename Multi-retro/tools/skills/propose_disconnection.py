"""
Propose Disconnection Skill
=============================

Thin skill wrapper around :func:`tools.guidance.disconnection_proposer.propose_disconnections`.
Accepts ``target_smiles`` (required), plus optional ``strategy_hint`` and
``max_proposals``, and returns a :class:`SkillResult` envelope.
"""

from __future__ import annotations

import logging
from typing import Any, Dict

from tools.skills.base import BaseSkill, SkillResult

logger = logging.getLogger(__name__)


class ProposeDisconnectionSkill(BaseSkill):
    """Propose retrosynthetic disconnection points for a target molecule."""

    name = "propose_disconnection"
    description = "Propose retrosynthetic disconnection points"

    def execute(self, args: Any) -> Dict[str, Any]:
        """Execute disconnection proposal.

        Args:
            args: dict with keys:
                - target_smiles (str): Required. SMILES of the target molecule.
                - strategy_hint (str): Optional strategy filter.
                - max_proposals (int): Optional cap on returned proposals.

        Returns:
            ``SkillResult.to_dict()`` envelope.
        """
        if not isinstance(args, dict) or "target_smiles" not in args:
            # Support namespace/dataclass args from skill_dispatch
            if hasattr(args, "target_smiles"):
                _d = {"target_smiles": args.target_smiles}
                for attr in ("strategy_hint", "max_proposals"):
                    if hasattr(args, attr):
                        _d[attr] = getattr(args, attr)
                args = _d
            else:
                return SkillResult(
                    success=False,
                    error="args must be a dict with 'target_smiles'",
                ).to_dict()

        target_smiles: str = args["target_smiles"]
        strategy_hint: str = args.get("strategy_hint", "")
        max_proposals: int = args.get("max_proposals", 5)

        try:
            from tools.guidance.disconnection_proposer import propose_disconnections

            result = propose_disconnections(
                target_smiles=target_smiles,
                strategy_hint=strategy_hint,
                max_proposals=max_proposals,
            )
        except Exception as exc:
            logger.exception("propose_disconnections failed")
            return SkillResult(success=False, error=str(exc)).to_dict()

        success = result.get("success", False)
        error = result.get("error", "")
        return SkillResult(success=success, data=result, error=error).to_dict()
