"""
Atom Mappers Module - LLM Focused

Provides atom mapping functionality for chemical reactions using RXNMapper.
Returns structured information for LLM analysis.

Main exports:
- RXNMapperWrapper: Main wrapper for atom mapping
- map_reaction: Convenience function for single reactions
- extract_atom_correspondence: Get atom-to-atom mapping
- get_reaction_center: Identify atoms that changed
- get_atom_changes: Detailed change analysis
- summarize_for_llm: Text summary for LLM consumption
"""

from .rxn_mapper import RXNMapperWrapper, map_reaction
from .utils import (
    extract_atom_correspondence,
    get_reaction_center,
    get_atom_changes,
    get_bond_changes,
    summarize_for_llm,
    remove_atom_numbers,
    validate_mapped_reaction
)

__all__ = [
    "RXNMapperWrapper",
    "map_reaction",
    "extract_atom_correspondence",
    "get_reaction_center",
    "get_atom_changes",
    "get_bond_changes",
    "summarize_for_llm",
    "remove_atom_numbers",
    "validate_mapped_reaction",
]

__version__ = "0.1.0"
