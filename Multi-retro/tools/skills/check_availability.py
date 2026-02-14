"""
Check Availability Skill
========================

Thin skill wrapper around :func:`tools.chem.sa_scorer.classify_sa`.
Accepts a list of SMILES strings, delegates to the tool module, and returns
a :class:`SkillResult` envelope with availability classifications.
"""

from __future__ import annotations

import json
import logging
from typing import Any, Dict, List

from tools.skills.base import BaseSkill, SkillResult

logger = logging.getLogger(__name__)


class CheckAvailabilitySkill(BaseSkill):
    """Check precursor availability based on SA score.

    Delegates to :func:`tools.chem.sa_scorer.classify_sa`.
    """

    name = "check_availability"
    description = "Check precursor availability based on SA score"

    def execute(self, args: Any) -> Dict[str, Any]:
        """Run availability check and return a SkillResult dict.

        Args:
            args: dict with ``smiles_list`` key (list of SMILES strings,
                  JSON array string, or comma-separated string).

        Returns:
            ``SkillResult.to_dict()`` with availability results under ``data``.
        """
        if not isinstance(args, dict) or "smiles_list" not in args:
            # Support namespace/dataclass args from skill_dispatch
            if hasattr(args, "smiles_list"):
                args = {"smiles_list": args.smiles_list}
            else:
                return SkillResult(
                    success=False, error="args must be a dict with 'smiles_list'"
                ).to_dict()

        smiles_list = _parse_smiles_list(args["smiles_list"])
        if not smiles_list:
            return SkillResult(
                success=True,
                data={"results": []},
            ).to_dict()

        try:
            from tools.chem.sa_scorer import classify_sa

            results = [classify_sa(smi) for smi in smiles_list]
        except Exception as exc:
            logger.exception("check_availability failed")
            return SkillResult(success=False, error=str(exc)).to_dict()

        return SkillResult(
            success=True,
            data={"results": results},
        ).to_dict()


def _parse_smiles_list(raw: Any) -> List[str]:
    """Normalise smiles_list from list, JSON string, or comma-separated string."""
    if isinstance(raw, list):
        return [s.strip() for s in raw if isinstance(s, str) and s.strip()]
    if isinstance(raw, str):
        stripped = raw.strip()
        if stripped.startswith("["):
            try:
                parsed = json.loads(stripped)
                return [s.strip() for s in parsed if isinstance(s, str) and s.strip()]
            except (json.JSONDecodeError, TypeError):
                pass
        return [s.strip() for s in stripped.split(",") if s.strip()]
    return []
