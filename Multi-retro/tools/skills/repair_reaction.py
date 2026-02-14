"""
Repair Reaction Skill
=====================

Attempts to fix invalid reactions by consulting the repair advisor for
guidance and the reaction validator for re-validation.

Delegates to:
    - ``tools.guidance.repair_advisor.generate_repair_guidance``
    - ``tools.chem.reaction_validator.validate_reaction``

Public API:
    RepairReactionSkill.execute(args) -> dict  (SkillResult.to_dict format)
"""

from __future__ import annotations

import logging
from typing import Any, Dict, List

from tools.skills.base import BaseSkill, SkillResult

logger = logging.getLogger(__name__)


class RepairReactionSkill(BaseSkill):
    """Repair an invalid reaction using guidance and re-validation.

    Args (dict or object):
        reaction_smiles – The reaction SMILES string to repair.
        issues          – List of issue strings or {code, message} dicts.

    Returns ``SkillResult.to_dict()`` with repair guidance and
    re-validation results.
    """

    name = "repair_reaction"
    description = "Repair reaction by analysing issues and suggesting fixes"

    def execute(self, args: Any) -> Dict[str, Any]:
        # ── Parse arguments ──
        if isinstance(args, dict):
            reaction_smiles = args.get("reaction_smiles", "")
            issues_raw = args.get("issues", [])
        else:
            reaction_smiles = getattr(args, "reaction_smiles", "")
            issues_raw = getattr(args, "issues", [])

        if not reaction_smiles:
            return SkillResult(
                success=False, error="reaction_smiles is required"
            ).to_dict()

        if not issues_raw:
            return SkillResult(
                success=False, error="issues list is required"
            ).to_dict()

        # Normalise issues to [{code, message}, ...]
        structured_issues = _normalise_issues(issues_raw)

        # ── Repair guidance ──
        try:
            from tools.guidance.repair_advisor import generate_repair_guidance

            guidance_result = generate_repair_guidance(structured_issues)
        except Exception as exc:
            logger.exception("Repair guidance failed")
            return SkillResult(
                success=False, error=f"Repair guidance failed: {exc}"
            ).to_dict()

        # ── Re-validate the reaction ──
        try:
            from tools.chem.reaction_validator import validate_reaction

            validation = validate_reaction(reaction_smiles)
            validation_dict = validation.to_dict()
        except Exception as exc:
            logger.exception("Re-validation failed")
            validation_dict = {"error": str(exc)}

        return SkillResult(
            success=True,
            data={
                "reaction_smiles": reaction_smiles,
                "issues": structured_issues,
                "guidance": guidance_result.get("guidance", []),
                "repair_strategy": guidance_result.get("repair_strategy", ""),
                "categories_affected": guidance_result.get("categories_affected", []),
                "total_issues": guidance_result.get("total_issues", 0),
                "critical_issues": guidance_result.get("critical_issues", 0),
                "revalidation": validation_dict,
            },
        ).to_dict()


def _normalise_issues(issues_raw: List[Any]) -> List[Dict[str, str]]:
    """Convert a mixed list of strings / dicts into [{code, message}, ...]."""
    structured: List[Dict[str, str]] = []
    for issue in issues_raw:
        if isinstance(issue, str):
            structured.append({"code": issue, "message": issue})
        elif isinstance(issue, dict):
            structured.append({
                "code": issue.get("code", "UNKNOWN"),
                "message": issue.get("message", ""),
            })
        else:
            structured.append({"code": "UNKNOWN", "message": str(issue)})
    return structured
