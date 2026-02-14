"""
Get Global Strategy Skill
=========================

Thin skill wrapper that generates global strategy context for LLM prompts.
Delegates to :func:`core.planner.global_strategy_manager.create_strategy_manager`
(via :mod:`tools.guidance.strategy_advisor`).
"""

from __future__ import annotations

import logging
from typing import Any, Dict

from tools.skills.base import BaseSkill, SkillResult

logger = logging.getLogger(__name__)


class GetGlobalStrategySkill(BaseSkill):
    """Get global strategy context for LLM prompts with PG policy and method preferences.

    Delegates to the global strategy manager to build strategy context.
    """

    name = "get_global_strategy"
    description = "Get global strategy context for LLM prompts with PG policy and method preferences"

    def execute(self, args: Any) -> Dict[str, Any]:
        """Generate global strategy context and return a SkillResult dict.

        Args:
            args: dict with keys:
                - target_smiles (str): Required. SMILES of the target molecule.
                - max_pg_operations (int): Optional. Defaults to 2.
                - max_steps (int): Optional. Defaults to 10.
                - constraints (dict): Optional. Extra constraints to apply.

        Returns:
            ``SkillResult.to_dict()`` with strategy context under ``data``.
        """
        if not isinstance(args, dict) or "target_smiles" not in args:
            # Support namespace/dataclass args from skill_dispatch
            if hasattr(args, "target_smiles"):
                _d = {"target_smiles": args.target_smiles}
                for attr in ("max_pg_operations", "max_steps", "constraints"):
                    if hasattr(args, attr):
                        _d[attr] = getattr(args, attr)
                args = _d
            else:
                return SkillResult(
                    success=False,
                    error="args must be a dict with 'target_smiles'",
                ).to_dict()

        target_smiles = args["target_smiles"]
        constraints: Dict[str, Any] = {
            "max_pg_operations": args.get("max_pg_operations", 2),
            "max_steps": args.get("max_steps", 10),
        }
        extra = args.get("constraints", {})
        if isinstance(extra, dict):
            constraints.update(extra)

        try:
            from tools.legacy.planner.global_strategy_manager import (
                PREFERRED_METHODS,
                create_strategy_manager,
            )

            manager = create_strategy_manager(constraints)
            strategy_context = manager.get_strategy_context()
        except Exception as exc:
            logger.exception("get_global_strategy failed")
            return SkillResult(success=False, error=str(exc)).to_dict()

        return SkillResult(
            success=True,
            data={
                "target_smiles": target_smiles,
                "strategy_context": strategy_context,
                "preferred_methods": PREFERRED_METHODS,
                "constraints_applied": constraints,
            },
        ).to_dict()
