"""
Map Atoms Skill
===============

Thin skill wrapper around :func:`tools.chem.atom_mapper.map_reaction`.
Accepts a reaction SMILES, delegates to the tool module, and returns
a :class:`SkillResult` envelope.
"""

from __future__ import annotations

import logging
from typing import Any, Dict

from tools.skills.base import BaseSkill, SkillResult

logger = logging.getLogger(__name__)


class MapAtomsSkill(BaseSkill):
    """Add atom-to-atom mapping using RXNMapper.

    Delegates to :func:`tools.chem.atom_mapper.map_reaction`.
    """

    name = "map_atoms"
    description = "Add atom-to-atom mapping using RXNMapper"

    def execute(self, args: Any) -> Dict[str, Any]:
        """Run atom mapping and return a SkillResult dict.

        Args:
            args: dict with ``reaction_smiles`` key, or any object with
                  a ``.reaction_smiles`` attribute.

        Returns:
            ``SkillResult.to_dict()`` with mapping data nested under ``data``.
        """
        reaction_smiles = _extract_reaction_smiles(args)
        if not reaction_smiles:
            return SkillResult(
                success=False, error="args must contain 'reaction_smiles'"
            ).to_dict()

        try:
            from tools.chem.atom_mapper import map_reaction

            result = map_reaction(reaction_smiles)
        except Exception as exc:
            logger.exception("map_reaction failed")
            return SkillResult(success=False, error=str(exc)).to_dict()

        success = result.get("success", False)
        error = result.get("error", "")
        return SkillResult(success=success, data=result, error=error).to_dict()


def _extract_reaction_smiles(args: Any) -> str:
    """Normalise various arg formats to a reaction SMILES string."""
    if isinstance(args, dict):
        return args.get("reaction_smiles", "")
    if hasattr(args, "reaction_smiles"):
        return args.reaction_smiles or ""
    return ""
