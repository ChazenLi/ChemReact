"""
Disconnection Proposer Module
==============================
Retrosynthetic disconnection proposals with 35+ reaction rules.
Migrated from core/planner/disconnector.py.

Performance fix: when strategy_hint doesn't match, falls back to full
rule traversal in a SINGLE pass (no duplicate traversal).

Public API:
    propose_disconnections(target_smiles, strategy_hint="", max_proposals=5) -> Dict
"""

from __future__ import annotations

import hashlib
import logging
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Tuple

from tools.common.constants import DEFAULT_MAX_PROPOSALS

logger = logging.getLogger(__name__)

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False


@dataclass
class DisconnectionProposal:
    """A proposed retrosynthetic disconnection."""
    id: str
    strategy: str
    bond_description: str
    precursors: List[str]
    reaction_class: str
    confidence: float
    rationale: str
    reaction_smiles: str = ""
    byproducts: List[str] = field(default_factory=list)
    reagents: List[str] = field(default_factory=list)
    is_redox: bool = False
    specificity: int = 1
    priority: int = 5

    def to_dict(self) -> Dict[str, Any]:
        return {
            "id": self.id, "strategy": self.strategy,
            "bond_description": self.bond_description,
            "precursors": self.precursors,
            "reaction_class": self.reaction_class,
            "confidence": self.confidence, "rationale": self.rationale,
            "reaction_smiles": self.reaction_smiles,
            "byproducts": self.byproducts, "reagents": self.reagents,
            "is_redox": self.is_redox,
        }


# Import rules from original disconnector module
from tools.legacy.planner.disconnector import (
    DISCONNECTION_RULES,
    STRATEGY_CATEGORY_MAP,
    _apply_disconnection_rule,
    _find_pattern_matches,
    _rule_category,
    _rule_category_counts,
)


def _generate_id(prefix: str, content: str) -> str:
    return f"{prefix}_{hashlib.sha256(content.encode()).hexdigest()[:8]}"


def propose_disconnections(
    target_smiles: str,
    strategy_hint: str = "",
    max_proposals: int = DEFAULT_MAX_PROPOSALS,
) -> Dict[str, Any]:
    """Propose retrosynthetic disconnections for a target molecule.

    Performance fix (Req 9.1): when strategy_hint is provided but yields
    no matches, we collect ALL proposals in a single pass rather than
    traversing rules twice.

    Args:
        target_smiles: Target molecule SMILES
        strategy_hint: Optional strategy hint to filter rules
        max_proposals: Maximum number of proposals to return

    Returns:
        {
            "success": bool,
            "target_smiles": str,
            "disconnections": [...],
            "total_rules": int,
            "rule_categories": {...},
        }
    """
    if not RDKIT_AVAILABLE:
        return {
            "success": False, "error": "RDKit not available",
            "target_smiles": target_smiles, "disconnections": [],
        }

    mol = Chem.MolFromSmiles(target_smiles)
    if mol is None:
        return {
            "success": False, "error": f"Invalid SMILES: {target_smiles}",
            "target_smiles": target_smiles, "disconnections": [],
        }

    canonical = Chem.MolToSmiles(mol, canonical=True)
    # Re-parse canonical SMILES so atom indices are deterministic
    # and consistent with analyze_molecule / break_bond.
    mol_c = Chem.MolFromSmiles(canonical)
    if mol_c is not None:
        mol = mol_c
    hint_proposals: List[DisconnectionProposal] = []
    all_proposals: List[DisconnectionProposal] = []
    seen: set = set()

    # Single pass through all rules
    for rule in DISCONNECTION_RULES:
        matches_hint = False
        if strategy_hint:
            hint_lower = strategy_hint.lower()
            if any(hint_lower in x.lower() for x in [
                rule["id"], rule["name"], rule.get("strategy", "")
            ]):
                matches_hint = True

        matches = _find_pattern_matches(mol, rule["pattern"])
        for match in matches:
            proposal = _apply_disconnection_rule(mol, canonical, rule, match)
            if proposal:
                key = tuple(sorted(proposal.precursors))
                if key not in seen:
                    seen.add(key)
                    if matches_hint:
                        hint_proposals.append(proposal)
                    all_proposals.append(proposal)

    # Use hint-filtered proposals if any matched, else fall back to all
    proposals = hint_proposals if (strategy_hint and hint_proposals) else all_proposals

    # Sort by priority ascending, confidence descending
    proposals.sort(key=lambda p: (p.priority, -p.confidence))
    proposals = proposals[:max_proposals]

    return {
        "success": True,
        "target_smiles": canonical,
        "disconnections": [p.to_dict() for p in proposals],
        "total_rules": len(DISCONNECTION_RULES),
        "rule_categories": _rule_category_counts(),
    }
