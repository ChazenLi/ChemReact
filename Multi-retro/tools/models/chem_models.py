"""
Chemistry data models — molecular structure representations.

Provides typed dataclasses for atoms, bonds, rings, atom-bond maps,
and full molecule info.  Every class implements ``to_dict()`` /
``from_dict()`` with round-trip consistency and missing-field tolerance.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional


# ---------------------------------------------------------------------------
# AtomInfo
# ---------------------------------------------------------------------------

@dataclass
class AtomInfo:
    """Single atom information."""

    idx: int = 0
    symbol: str = ""
    aromatic: bool = False
    in_ring: bool = False
    degree: int = 0
    num_hs: int = 0
    neighbors: List[int] = field(default_factory=list)

    def to_dict(self) -> Dict[str, Any]:
        return {
            "idx": self.idx,
            "symbol": self.symbol,
            "aromatic": self.aromatic,
            "in_ring": self.in_ring,
            "degree": self.degree,
            "num_hs": self.num_hs,
            "neighbors": list(self.neighbors),
        }

    @classmethod
    def from_dict(cls, d: Dict[str, Any]) -> "AtomInfo":
        return cls(
            idx=d.get("idx", 0),
            symbol=d.get("symbol", ""),
            aromatic=d.get("aromatic", False),
            in_ring=d.get("in_ring", False),
            degree=d.get("degree", 0),
            num_hs=d.get("num_hs", 0),
            neighbors=list(d.get("neighbors", [])),
        )


# ---------------------------------------------------------------------------
# BondInfo
# ---------------------------------------------------------------------------

@dataclass
class BondInfo:
    """Single chemical bond information."""

    idx: int = 0
    atom1_idx: int = 0
    atom2_idx: int = 0
    bond_type: str = "SINGLE"  # "SINGLE", "DOUBLE", "TRIPLE", "AROMATIC"
    in_ring: bool = False
    breakable: bool = False
    chemical_label: str = ""
    potential_reactions: List[str] = field(default_factory=list)
    retro_hint: str = ""
    strategic_priority: int = 3  # 1=骨架, 2=环构建, 3=官能团

    def to_dict(self) -> Dict[str, Any]:
        return {
            "idx": self.idx,
            "atom1_idx": self.atom1_idx,
            "atom2_idx": self.atom2_idx,
            "bond_type": self.bond_type,
            "in_ring": self.in_ring,
            "breakable": self.breakable,
            "chemical_label": self.chemical_label,
            "potential_reactions": list(self.potential_reactions),
            "retro_hint": self.retro_hint,
            "strategic_priority": self.strategic_priority,
        }

    @classmethod
    def from_dict(cls, d: Dict[str, Any]) -> "BondInfo":
        return cls(
            idx=d.get("idx", 0),
            atom1_idx=d.get("atom1_idx", 0),
            atom2_idx=d.get("atom2_idx", 0),
            bond_type=d.get("bond_type", "SINGLE"),
            in_ring=d.get("in_ring", False),
            breakable=d.get("breakable", False),
            chemical_label=d.get("chemical_label", ""),
            potential_reactions=list(d.get("potential_reactions", [])),
            retro_hint=d.get("retro_hint", ""),
            strategic_priority=d.get("strategic_priority", 3),
        )


# ---------------------------------------------------------------------------
# RingInfo
# ---------------------------------------------------------------------------

@dataclass
class RingInfo:
    """Ring system information."""

    ring_id: int = 0
    size: int = 0
    atom_indices: List[int] = field(default_factory=list)
    is_aromatic: bool = False
    ring_type: str = ""
    formation_hints: List[str] = field(default_factory=list)

    def to_dict(self) -> Dict[str, Any]:
        return {
            "ring_id": self.ring_id,
            "size": self.size,
            "atom_indices": list(self.atom_indices),
            "is_aromatic": self.is_aromatic,
            "ring_type": self.ring_type,
            "formation_hints": list(self.formation_hints),
        }

    @classmethod
    def from_dict(cls, d: Dict[str, Any]) -> "RingInfo":
        return cls(
            ring_id=d.get("ring_id", 0),
            size=d.get("size", 0),
            atom_indices=list(d.get("atom_indices", [])),
            is_aromatic=d.get("is_aromatic", False),
            ring_type=d.get("ring_type", ""),
            formation_hints=list(d.get("formation_hints", [])),
        )


# ---------------------------------------------------------------------------
# AtomBondMap
# ---------------------------------------------------------------------------

@dataclass
class AtomBondMap:
    """Atom-bond map — core output of analyze_molecule."""

    canonical_smiles: str = ""
    atoms: List[AtomInfo] = field(default_factory=list)
    bonds: List[BondInfo] = field(default_factory=list)
    ring_info: List[RingInfo] = field(default_factory=list)
    retro_guidance: str = ""

    def to_dict(self) -> Dict[str, Any]:
        return {
            "canonical_smiles": self.canonical_smiles,
            "atoms": [a.to_dict() for a in self.atoms],
            "bonds": [b.to_dict() for b in self.bonds],
            "ring_info": [r.to_dict() for r in self.ring_info],
            "retro_guidance": self.retro_guidance,
        }

    @classmethod
    def from_dict(cls, d: Dict[str, Any]) -> "AtomBondMap":
        return cls(
            canonical_smiles=d.get("canonical_smiles", ""),
            atoms=[AtomInfo.from_dict(a) for a in d.get("atoms", [])],
            bonds=[BondInfo.from_dict(b) for b in d.get("bonds", [])],
            ring_info=[RingInfo.from_dict(r) for r in d.get("ring_info", [])],
            retro_guidance=d.get("retro_guidance", ""),
        )


# ---------------------------------------------------------------------------
# MoleculeInfo
# ---------------------------------------------------------------------------

@dataclass
class MoleculeInfo:
    """Comprehensive molecule information."""

    smiles: str = ""
    canonical_smiles: str = ""
    molecular_weight: float = 0.0
    sa_score: float = 0.0
    atom_bond_map: Optional[AtomBondMap] = None
    functional_groups: Dict[str, int] = field(default_factory=dict)
    protecting_groups: Dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> Dict[str, Any]:
        return {
            "smiles": self.smiles,
            "canonical_smiles": self.canonical_smiles,
            "molecular_weight": self.molecular_weight,
            "sa_score": self.sa_score,
            "atom_bond_map": self.atom_bond_map.to_dict() if self.atom_bond_map else None,
            "functional_groups": dict(self.functional_groups),
            "protecting_groups": dict(self.protecting_groups),
        }

    @classmethod
    def from_dict(cls, d: Dict[str, Any]) -> "MoleculeInfo":
        abm_raw = d.get("atom_bond_map")
        return cls(
            smiles=d.get("smiles", ""),
            canonical_smiles=d.get("canonical_smiles", ""),
            molecular_weight=d.get("molecular_weight", 0.0),
            sa_score=d.get("sa_score", 0.0),
            atom_bond_map=AtomBondMap.from_dict(abm_raw) if abm_raw else None,
            functional_groups=dict(d.get("functional_groups", {})),
            protecting_groups=dict(d.get("protecting_groups", {})),
        )
