"""
Atom Mapping Module
===================
Wrapper around RXNMapper for atom mapping chemical reactions.
Migrated from internal/atom_mappers/rxn_mapper.py.

Public API:
    RXNMapperWrapper  - class with map_reaction() method
    map_reaction()    - convenience function for single reaction mapping
"""

import logging
import re
from typing import Any, Dict, List

logger = logging.getLogger(__name__)

try:
    from rxnmapper import RXNMapper as _RXNMapper
    RXNMAPPER_AVAILABLE = True
except ImportError:
    RXNMAPPER_AVAILABLE = False
    _RXNMapper = None


class RXNMapperWrapper:
    """Lightweight wrapper for RXNMapper focused on structured data output."""

    def __init__(self):
        self._mapper = _RXNMapper() if RXNMAPPER_AVAILABLE else None

    def is_available(self) -> bool:
        return RXNMAPPER_AVAILABLE and self._mapper is not None

    def map_reaction(self, reaction_smiles: str) -> Dict[str, Any]:
        """Map atoms in a reaction and return structured data.

        Args:
            reaction_smiles: Reaction SMILES (reactants>>products)

        Returns:
            {
                "success": bool,
                "mapped_rxn": str or None,
                "confidence": float,
                "reaction_center": List[int],
                "correspondence": Dict[int, Dict],
                "original": str
            }
        """
        result: Dict[str, Any] = {
            "success": False,
            "mapped_rxn": None,
            "confidence": 0.0,
            "reaction_center": [],
            "correspondence": {},
            "original": reaction_smiles,
        }

        if not self.is_available():
            return result

        try:
            mapped_list = self._mapper.get_attention_guided_atom_maps(
                [reaction_smiles]
            )
            if not mapped_list:
                return result

            mapped = mapped_list[0]
            result["mapped_rxn"] = mapped.get("mapped_rxn")
            result["confidence"] = mapped.get("confidence", 0.0)
            result["success"] = True

            if result["mapped_rxn"]:
                correspondence = self._extract_correspondence(result["mapped_rxn"])
                result["correspondence"] = correspondence
                result["reaction_center"] = self._find_reaction_center(correspondence)
        except Exception as e:
            logger.error("Mapping error: %s", e)

        return result

    def _extract_correspondence(
        self, mapped_smiles: str
    ) -> Dict[int, Dict[str, str]]:
        """Extract atom correspondence from mapped SMILES."""
        if ">>" not in mapped_smiles:
            return {}

        reactants, products = mapped_smiles.split(">>", 1)
        pattern = re.compile(r'\[([A-Za-z][A-Za-z0-9+@-]*):(\d+)\]')

        reactant_atoms: Dict[int, str] = {}
        for match in pattern.finditer(reactants):
            atom_num = int(match.group(2))
            if atom_num not in reactant_atoms:
                reactant_atoms[atom_num] = match.group(1)

        product_atoms: Dict[int, str] = {}
        for match in pattern.finditer(products):
            atom_num = int(match.group(2))
            if atom_num not in product_atoms:
                product_atoms[atom_num] = match.group(1)

        correspondence: Dict[int, Dict[str, str]] = {}
        for atom_num, r_type in reactant_atoms.items():
            correspondence[atom_num] = {
                "reactant": r_type,
                "product": product_atoms.get(atom_num, "?"),
            }
        return correspondence

    def _find_reaction_center(
        self, correspondence: Dict[int, Dict[str, str]]
    ) -> List[int]:
        """Identify atoms that changed (reaction center)."""
        return [
            num for num, types in correspondence.items()
            if types["reactant"] != types["product"]
        ]


def map_reaction(reaction_smiles: str) -> Dict[str, Any]:
    """Convenience function for single reaction mapping."""
    return RXNMapperWrapper().map_reaction(reaction_smiles)
