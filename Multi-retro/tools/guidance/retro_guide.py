"""
Retro Analysis Guide Builder
=============================
Reaction-aware retro analysis guide for LLM decision-making.
Extracted from BreakBondSkill in core/skills/builtin.py.

Builds structured guide dicts describing bond semantics, fragment profiles,
and retro hints for LLM to construct correct precursor SMILES.

Public API:
    build_retro_guide(mol, bond, atom1_idx, atom2_idx, fragments,
                      bond_semantic, retro_pattern, retro_hint,
                      extra_info=None) -> Dict[str, Any]
    build_fragment_profiles(fragments) -> List[Dict[str, Any]]
"""

from __future__ import annotations

import logging
from typing import Any, Dict, List, Optional

logger = logging.getLogger(__name__)

try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False


def build_fragment_profiles(fragments: List[str]) -> List[Dict[str, Any]]:
    """Build chemical profiles for each fragment SMILES.

    Returns list of dicts with smiles, atom_count, has_aromatic,
    has_heteroatom, molecular_weight, functional_groups_hint.
    """
    profiles: List[Dict[str, Any]] = []
    for frag_smi in fragments:
        profile: Dict[str, Any] = {"smiles": frag_smi}
        if not RDKIT_AVAILABLE:
            profiles.append(profile)
            continue

        mol = Chem.MolFromSmiles(frag_smi)
        if mol is None:
            profiles.append(profile)
            continue

        profile["atom_count"] = mol.GetNumAtoms()
        profile["has_aromatic"] = any(a.GetIsAromatic() for a in mol.GetAtoms())
        profile["has_heteroatom"] = any(
            a.GetSymbol() not in ("C", "H") for a in mol.GetAtoms()
        )
        try:
            profile["molecular_weight"] = round(Descriptors.MolWt(mol), 1)
        except Exception:
            pass

        # Quick FG hints
        fg_hints: List[str] = []
        for a in mol.GetAtoms():
            sym = a.GetSymbol()
            if sym == "N":
                fg_hints.append("nitrogen-containing")
            elif sym == "O":
                fg_hints.append("oxygen-containing")
            elif sym == "S":
                fg_hints.append("sulfur-containing")
        profile["functional_groups_hint"] = list(set(fg_hints))
        profiles.append(profile)

    return profiles


def build_retro_guide(
    mol,
    bond,
    atom1_idx: int,
    atom2_idx: int,
    fragments: List[str],
    bond_semantic: str,
    retro_pattern: str,
    retro_hint: str,
    extra_info: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    """Build a retro analysis guide for LLM decision-making.

    Args:
        mol: RDKit Mol object of the target molecule.
        bond: RDKit Bond object that was broken.
        atom1_idx: Index of first atom in the broken bond.
        atom2_idx: Index of second atom in the broken bond.
        fragments: List of fragment SMILES after bond breaking.
        bond_semantic: Human-readable bond type description.
        retro_pattern: Name of the retrosynthetic pattern.
        retro_hint: Specific hint for LLM on how to handle this disconnection.
        extra_info: Optional additional info to merge into the guide.

    Returns:
        Structured guide dict with bond info, fragment profiles, and LLM task.
    """
    atom1 = mol.GetAtomWithIdx(atom1_idx)
    atom2 = mol.GetAtomWithIdx(atom2_idx)

    guide: Dict[str, Any] = {
        "bond_type": bond_semantic,
        "retro_pattern": retro_pattern,
        "broken_atoms": {
            "atom1": {
                "idx": atom1_idx,
                "symbol": atom1.GetSymbol(),
                "in_ring": atom1.IsInRing(),
                "is_aromatic": atom1.GetIsAromatic(),
            },
            "atom2": {
                "idx": atom2_idx,
                "symbol": atom2.GetSymbol(),
                "in_ring": atom2.IsInRing(),
                "is_aromatic": atom2.GetIsAromatic(),
            },
        },
        "fragment_profiles": build_fragment_profiles(fragments),
        "retro_hint": retro_hint,
        "llm_task": (
            f"Detected {retro_pattern} ({bond_semantic}). "
            f"Use fragment_profiles and retro_hint to determine each fragment's role "
            f"(nucleophile/electrophile), correct SMILES (add leaving/activating groups), "
            f"build retro reaction: precursors >> product, then validate with validate_reaction."
        ),
    }

    if extra_info:
        guide.update(extra_info)

    return guide
