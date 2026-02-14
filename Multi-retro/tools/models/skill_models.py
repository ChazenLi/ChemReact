"""
Skill input/output data models.

Typed dataclasses for Skill arguments and results, replacing the
previous ``SimpleNamespace`` / loose-dict approach.  Every class
implements ``to_dict()`` / ``from_dict()`` with round-trip consistency
and missing-field tolerance.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional

from .chem_models import MoleculeInfo


# ---------------------------------------------------------------------------
# Analyze
# ---------------------------------------------------------------------------

@dataclass
class AnalyzeArgs:
    """Input for AnalyzeMoleculeSkill."""

    smiles: str = ""

    def to_dict(self) -> Dict[str, Any]:
        return {"smiles": self.smiles}

    @classmethod
    def from_dict(cls, d: Dict[str, Any]) -> "AnalyzeArgs":
        return cls(smiles=d.get("smiles", ""))


@dataclass
class AnalysisResult:
    """Output of AnalyzeMoleculeSkill."""

    success: bool = False
    molecule_info: Optional[MoleculeInfo] = None
    error: str = ""

    def to_dict(self) -> Dict[str, Any]:
        return {
            "success": self.success,
            "molecule_info": self.molecule_info.to_dict() if self.molecule_info else None,
            "error": self.error,
        }

    @classmethod
    def from_dict(cls, d: Dict[str, Any]) -> "AnalysisResult":
        mi_raw = d.get("molecule_info")
        return cls(
            success=d.get("success", False),
            molecule_info=MoleculeInfo.from_dict(mi_raw) if mi_raw else None,
            error=d.get("error", ""),
        )


# ---------------------------------------------------------------------------
# Break Bond / Disconnection
# ---------------------------------------------------------------------------

@dataclass
class BreakBondArgs:
    """Input for BreakBondSkill."""

    smiles: str = ""
    atom1_idx: int = 0
    atom2_idx: int = 0

    def to_dict(self) -> Dict[str, Any]:
        return {
            "smiles": self.smiles,
            "atom1_idx": self.atom1_idx,
            "atom2_idx": self.atom2_idx,
        }

    @classmethod
    def from_dict(cls, d: Dict[str, Any]) -> "BreakBondArgs":
        return cls(
            smiles=d.get("smiles", ""),
            atom1_idx=d.get("atom1_idx", 0),
            atom2_idx=d.get("atom2_idx", 0),
        )


@dataclass
class DisconnectionResult:
    """Output of BreakBondSkill / ProposeDisconnectionSkill."""

    success: bool = False
    fragments: List[str] = field(default_factory=list)
    reaction_smiles: str = ""
    retro_analysis_guide: Optional[Dict[str, Any]] = None
    error: str = ""

    def to_dict(self) -> Dict[str, Any]:
        return {
            "success": self.success,
            "fragments": list(self.fragments),
            "reaction_smiles": self.reaction_smiles,
            "retro_analysis_guide": self.retro_analysis_guide,
            "error": self.error,
        }

    @classmethod
    def from_dict(cls, d: Dict[str, Any]) -> "DisconnectionResult":
        return cls(
            success=d.get("success", False),
            fragments=list(d.get("fragments", [])),
            reaction_smiles=d.get("reaction_smiles", ""),
            retro_analysis_guide=d.get("retro_analysis_guide"),
            error=d.get("error", ""),
        )


# ---------------------------------------------------------------------------
# Validate
# ---------------------------------------------------------------------------

@dataclass
class ValidateArgs:
    """Input for ValidateReactionSkill."""

    reaction_smiles: str = ""

    def to_dict(self) -> Dict[str, Any]:
        return {"reaction_smiles": self.reaction_smiles}

    @classmethod
    def from_dict(cls, d: Dict[str, Any]) -> "ValidateArgs":
        return cls(reaction_smiles=d.get("reaction_smiles", ""))


@dataclass
class ValidationResult:
    """Output of ValidateReactionSkill."""

    success: bool = False
    verdict: str = ""  # "PASS", "CONDITIONAL", "FAIL"
    is_valid: bool = False
    issues: List[Dict[str, Any]] = field(default_factory=list)
    byproducts: List[Dict[str, Any]] = field(default_factory=list)
    error: str = ""

    def to_dict(self) -> Dict[str, Any]:
        return {
            "success": self.success,
            "verdict": self.verdict,
            "is_valid": self.is_valid,
            "issues": list(self.issues),
            "byproducts": list(self.byproducts),
            "error": self.error,
        }

    @classmethod
    def from_dict(cls, d: Dict[str, Any]) -> "ValidationResult":
        return cls(
            success=d.get("success", False),
            verdict=d.get("verdict", ""),
            is_valid=d.get("is_valid", False),
            issues=list(d.get("issues", [])),
            byproducts=list(d.get("byproducts", [])),
            error=d.get("error", ""),
        )


# ---------------------------------------------------------------------------
# Repair
# ---------------------------------------------------------------------------

@dataclass
class RepairArgs:
    """Input for RepairReactionSkill."""

    reaction_smiles: str = ""
    issues: List[str] = field(default_factory=list)

    def to_dict(self) -> Dict[str, Any]:
        return {
            "reaction_smiles": self.reaction_smiles,
            "issues": list(self.issues),
        }

    @classmethod
    def from_dict(cls, d: Dict[str, Any]) -> "RepairArgs":
        return cls(
            reaction_smiles=d.get("reaction_smiles", ""),
            issues=list(d.get("issues", [])),
        )


@dataclass
class RepairResult:
    """Output of RepairReactionSkill."""

    success: bool = False
    repaired_smiles: str = ""
    repair_actions: List[str] = field(default_factory=list)
    error: str = ""

    def to_dict(self) -> Dict[str, Any]:
        return {
            "success": self.success,
            "repaired_smiles": self.repaired_smiles,
            "repair_actions": list(self.repair_actions),
            "error": self.error,
        }

    @classmethod
    def from_dict(cls, d: Dict[str, Any]) -> "RepairResult":
        return cls(
            success=d.get("success", False),
            repaired_smiles=d.get("repaired_smiles", ""),
            repair_actions=list(d.get("repair_actions", [])),
            error=d.get("error", ""),
        )


# ---------------------------------------------------------------------------
# Availability
# ---------------------------------------------------------------------------

@dataclass
class AvailabilityResult:
    """Output of CheckAvailabilitySkill."""

    success: bool = False
    precursors: List[Dict[str, Any]] = field(default_factory=list)
    # Each precursor: {smiles, sa_score, classification, is_terminal}
    error: str = ""

    def to_dict(self) -> Dict[str, Any]:
        return {
            "success": self.success,
            "precursors": list(self.precursors),
            "error": self.error,
        }

    @classmethod
    def from_dict(cls, d: Dict[str, Any]) -> "AvailabilityResult":
        return cls(
            success=d.get("success", False),
            precursors=list(d.get("precursors", [])),
            error=d.get("error", ""),
        )
