"""
RXNMapper Wrapper Module

Provides a wrapper around RXNMapper for atom mapping chemical reactions.
Designed to return structured information for LLM consumption.

Key outputs for LLM:
- mapped_rxn: Reaction SMILES with atom map numbers
- confidence: Mapping confidence score
- reaction_center: Atoms that changed in the reaction
- atom_correspondence: Reactant -> Product atom type mapping
"""

import logging
from typing import Dict, List, Optional, Any

logger = logging.getLogger(__name__)

# Try to import RXNMapper
try:
    from rxnmapper import RXNMapper as _RXNMapper
    RXNMAPPER_AVAILABLE = True
except ImportError:
    RXNMAPPER_AVAILABLE = False
    _RXNMapper = None


class RXNMapperWrapper:
    """
    Lightweight wrapper for RXNMapper focused on LLM data consumption.
    """

    def __init__(self):
        if RXNMAPPER_AVAILABLE:
            self._mapper = _RXNMapper()
        else:
            self._mapper = None

    def is_available(self) -> bool:
        return RXNMAPPER_AVAILABLE and self._mapper is not None

    def map_reaction(self, reaction_smiles: str) -> Dict[str, Any]:
        """
        Map atoms in a reaction and return structured data for LLM.

        Args:
            reaction_smiles: Reaction SMILES (reactants>>products)

        Returns:
            {
                "success": bool,
                "mapped_rxn": str or None,
                "confidence": float,
                "reaction_center": List[int],  # atoms that changed
                "correspondence": Dict[int, Dict],  # atom -> {reactant, product} types
                "original": str
            }
        """
        result = {
            "success": False,
            "mapped_rxn": None,
            "confidence": 0.0,
            "reaction_center": [],
            "correspondence": {},
            "original": reaction_smiles
        }

        if not self.is_available():
            return result

        try:
            # Get mapping from RXNMapper
            mapped_list = self._mapper.get_attention_guided_atom_maps([reaction_smiles])

            if not mapped_list:
                return result

            mapped = mapped_list[0]
            result["mapped_rxn"] = mapped.get("mapped_rxn")
            result["confidence"] = mapped.get("confidence", 0.0)
            result["success"] = True

            # Extract structured information for LLM
            if result["mapped_rxn"]:
                correspondence = self._extract_correspondence(result["mapped_rxn"])
                result["correspondence"] = correspondence
                result["reaction_center"] = self._find_reaction_center(correspondence)

        except Exception as e:
            logger.error(f"Mapping error: {e}")

        return result

    def _extract_correspondence(self, mapped_smiles: str) -> Dict[int, Dict[str, str]]:
        """Extract atom correspondence from mapped SMILES."""
        import re

        if ">>" not in mapped_smiles:
            return {}

        reactants, products = mapped_smiles.split(">>", 1)
        pattern = re.compile(r'\[([A-Za-z][A-Za-z0-9+@-]*):(\d+)\]')

        reactant_atoms = {}
        for match in pattern.finditer(reactants):
            atom_type = match.group(1)
            atom_num = int(match.group(2))
            if atom_num not in reactant_atoms:
                reactant_atoms[atom_num] = atom_type

        product_atoms = {}
        for match in pattern.finditer(products):
            atom_type = match.group(1)
            atom_num = int(match.group(2))
            if atom_num not in product_atoms:
                product_atoms[atom_num] = atom_type

        # Build correspondence
        correspondence = {}
        for atom_num, r_type in reactant_atoms.items():
            correspondence[atom_num] = {
                "reactant": r_type,
                "product": product_atoms.get(atom_num, "?")
            }

        return correspondence

    def _find_reaction_center(self, correspondence: Dict[int, Dict[str, str]]) -> List[int]:
        """Identify atoms that changed (reaction center)."""
        reaction_center = []
        for atom_num, types in correspondence.items():
            if types["reactant"] != types["product"]:
                reaction_center.append(atom_num)
        return reaction_center


def map_reaction(reaction_smiles: str) -> Dict[str, Any]:
    """Convenience function for single reaction mapping."""
    mapper = RXNMapperWrapper()
    return mapper.map_reaction(reaction_smiles)
