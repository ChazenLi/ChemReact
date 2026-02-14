"""
Analyze Molecule Skill
======================

Thin skill wrapper around :mod:`tools.chem.structure_analyzer`.
Accepts a SMILES string, delegates to the tool module, and returns
a :class:`SkillResult` envelope.
"""

from __future__ import annotations

import logging
from typing import Any, Dict

from tools.skills.base import BaseSkill, SkillResult
from tools.models import AnalyzeArgs, AnalysisResult

logger = logging.getLogger(__name__)


class AnalyzeMoleculeSkill(BaseSkill):
    """Analyze target molecule structure, functional groups, and properties.

    Delegates heavy lifting to :func:`tools.chem.structure_analyzer.analyze_molecule`.
    """

    name = "analyze_molecule"
    description = "Analyze target molecule structure, functional groups, properties"

    def execute(self, args: Any) -> Dict[str, Any]:
        """Run molecule analysis and return a SkillResult dict.

        Args:
            args: Either an :class:`AnalyzeArgs` instance, a dict with a
                  ``smiles`` key, or any object with a ``.smiles`` attribute.

        Returns:
            ``SkillResult.to_dict()`` with analysis data nested under ``data``.
        """
        from tools.chem.structure_analyzer import analyze_molecule

        # Normalise input to a SMILES string
        if isinstance(args, dict):
            smiles = args.get("smiles", "")
        elif isinstance(args, AnalyzeArgs):
            smiles = args.smiles
        elif hasattr(args, "smiles"):
            smiles = args.smiles
        else:
            return SkillResult(
                success=False, error="Invalid args: expected 'smiles'"
            ).to_dict()

        if not smiles or not smiles.strip():
            return SkillResult(
                success=False, error="Empty SMILES string"
            ).to_dict()

        try:
            result: AnalysisResult = analyze_molecule(smiles)
        except Exception as exc:
            logger.exception("analyze_molecule failed for %s", smiles)
            return SkillResult(success=False, error=str(exc)).to_dict()

        if not result.success:
            return SkillResult(
                success=False, error=result.error or "Analysis failed"
            ).to_dict()

        return SkillResult(
            success=True,
            data=result.to_dict(),
        ).to_dict()
