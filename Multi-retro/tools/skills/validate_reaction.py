"""
Validate Reaction Skill
=======================

Thin skill wrapper around :mod:`tools.chem.reaction_validator`.
Accepts a reaction SMILES string, delegates to the tool module, and
returns a :class:`SkillResult` envelope.
"""

from __future__ import annotations

import logging
from typing import Any, Dict

from tools.skills.base import BaseSkill, SkillResult
from tools.models import ValidateArgs, ValidationResult

logger = logging.getLogger(__name__)


class ValidateReactionSkill(BaseSkill):
    """Validate reaction for atom balance and chemical plausibility.

    Delegates heavy lifting to
    :func:`tools.chem.reaction_validator.validate_reaction`.
    """

    name = "validate_reaction"
    description = "Validate reaction for atom balance and chemical plausibility"

    def execute(self, args: Any) -> Dict[str, Any]:
        """Run reaction validation and return a SkillResult dict.

        Args:
            args: Either a :class:`ValidateArgs` instance, a dict with a
                  ``reaction_smiles`` key, or any object with a
                  ``.reaction_smiles`` attribute.

        Returns:
            ``SkillResult.to_dict()`` with validation data nested under ``data``.
        """
        from tools.chem.reaction_validator import validate_reaction

        # Normalise input
        if isinstance(args, dict):
            reaction_smiles = args.get("reaction_smiles", "")
            reaction_class = args.get("reaction_class", "")
        elif isinstance(args, ValidateArgs):
            reaction_smiles = args.reaction_smiles
            reaction_class = ""
        elif hasattr(args, "reaction_smiles"):
            reaction_smiles = args.reaction_smiles
            reaction_class = getattr(args, "reaction_class", "")
        else:
            return SkillResult(
                success=False, error="Invalid args: expected 'reaction_smiles'"
            ).to_dict()

        if not reaction_smiles or not reaction_smiles.strip():
            return SkillResult(
                success=False, error="Empty reaction SMILES string"
            ).to_dict()

        try:
            result: ValidationResult = validate_reaction(
                reaction_smiles=reaction_smiles,
                reaction_class=reaction_class,
            )
        except Exception as exc:
            logger.exception("validate_reaction failed for %s", reaction_smiles)
            return SkillResult(success=False, error=str(exc)).to_dict()

        if not result.success:
            return SkillResult(
                success=False,
                data=result.to_dict(),
                error=result.error or "Validation failed",
            ).to_dict()

        return SkillResult(
            success=True,
            data=result.to_dict(),
        ).to_dict()
