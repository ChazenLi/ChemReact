"""
Molecular Parser & SMILES Utilities
====================================
Consolidates internal/common/mol_cache.py and internal/normalize/smiles_normalizer.py
into a single, cached, type-safe module.

Public API:
    parse_smiles(smiles)          → Optional[Mol]   (LRU-cached)
    normalize_smiles(smiles)      → str              (canonical SMILES)
    normalize_smiles_full(smiles) → Dict[str, str]   (original/canonical/display)
    validate_smiles(smiles)       → Tuple[bool, str] (safety checks)
    remove_atom_numbers(smiles)   → str
    normalize_reaction_smiles(reaction_smiles) → Dict[str, str]
    normalize_smiles_list(smiles_list)          → List[Dict[str, str]]
    clear_cache()                 → None
"""

import re
from functools import lru_cache
from typing import Dict, List, Optional, Tuple

from rdkit import Chem

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

_MAX_SMILES_LENGTH = 5000
_MAX_NESTING_DEPTH = 40  # max parenthesis / bracket nesting

# Characters that may legitimately appear in a SMILES string.
# Covers atoms, bonds, branches, rings, charges, chirality, isotopes, etc.
_VALID_SMILES_CHARS = re.compile(
    r'^[A-Za-z0-9@+\-\[\]\(\)\.\#\=\:\%\/\\]+$'
)


# ---------------------------------------------------------------------------
# Core cached helpers (from mol_cache.py)
# ---------------------------------------------------------------------------

@lru_cache(maxsize=512)
def parse_smiles(smiles: str) -> Optional[Chem.Mol]:
    """Parse a SMILES string into an RDKit Mol object with sanitization.

    Results are LRU-cached (maxsize=512) so repeated calls with the same
    SMILES return the *same* object without re-parsing.

    Returns ``None`` for invalid or empty SMILES.
    """
    if not smiles:
        return None
    mol = Chem.MolFromSmiles(smiles)
    if mol is not None:
        try:
            Chem.SanitizeMol(mol)
        except Exception:
            return None
    return mol


# ---------------------------------------------------------------------------
# Normalization (from smiles_normalizer.py — simplified + full variants)
# ---------------------------------------------------------------------------

def normalize_smiles(smiles: str) -> str:
    """Return the canonical SMILES string for *smiles*.

    This is the simplified interface: it returns a plain ``str`` rather than
    the three-format dict produced by :func:`normalize_smiles_full`.
    Invalid or empty input is returned unchanged.
    """
    if not smiles:
        return smiles
    mol = parse_smiles(smiles)
    if mol is None:
        return smiles
    return Chem.MolToSmiles(mol)


def normalize_smiles_full(smiles: str) -> Dict[str, str]:
    """Normalize *smiles* into three formats (backward-compatible dict).

    Returns::

        {
            "original":  <input as-is>,
            "canonical": <RDKit canonical, mapping preserved>,
            "display":   <atom-mapping stripped, for visualisation>,
        }
    """
    result: Dict[str, str] = {
        "original": smiles,
        "canonical": smiles,
        "display": smiles,
    }
    if not smiles:
        return result

    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return result

        result["canonical"] = Chem.MolToSmiles(mol)

        display_smiles = remove_atom_numbers(smiles)
        display_mol = Chem.MolFromSmiles(display_smiles)
        if display_mol:
            result["display"] = Chem.MolToSmiles(display_mol)
        else:
            result["display"] = display_smiles
    except Exception:
        pass

    return result


# ---------------------------------------------------------------------------
# Validation
# ---------------------------------------------------------------------------

def validate_smiles(smiles: str) -> Tuple[bool, str]:
    """Validate a SMILES string for safety and correctness.

    Checks performed (in order):
    1. Non-empty
    2. Length ≤ ``_MAX_SMILES_LENGTH``
    3. Only valid SMILES characters
    4. Nesting depth (parentheses + brackets) ≤ ``_MAX_NESTING_DEPTH``
    5. RDKit parsability

    Returns:
        ``(True, "")`` on success, or ``(False, <reason>)`` on failure.
    """
    if not smiles or not smiles.strip():
        return False, "SMILES string is empty"

    if len(smiles) > _MAX_SMILES_LENGTH:
        return False, f"SMILES exceeds maximum length ({len(smiles)} > {_MAX_SMILES_LENGTH})"

    if not _VALID_SMILES_CHARS.match(smiles):
        return False, "SMILES contains invalid characters"

    # Nesting depth check
    depth = 0
    max_depth = 0
    for ch in smiles:
        if ch in ("(", "["):
            depth += 1
            max_depth = max(max_depth, depth)
        elif ch in (")", "]"):
            depth -= 1
    if max_depth > _MAX_NESTING_DEPTH:
        return False, f"SMILES nesting depth too deep ({max_depth} > {_MAX_NESTING_DEPTH})"

    # RDKit parsability
    mol = parse_smiles(smiles)
    if mol is None:
        return False, "RDKit cannot parse this SMILES"

    return True, ""


# ---------------------------------------------------------------------------
# Atom-mapping helpers (from smiles_normalizer.py)
# ---------------------------------------------------------------------------

def remove_atom_numbers(smiles: str) -> str:
    """Strip atom-mapping numbers from a SMILES string.

    Example::

        >>> remove_atom_numbers("[CH3:1][OH:2]")
        '[CH3][OH]'
    """
    return re.sub(r'\[([A-Za-z][A-Za-z0-9+@-]*):\d+\]', r'[\1]', smiles)


def normalize_reaction_smiles(reaction_smiles: str) -> Dict[str, str]:
    """Normalize a reaction SMILES (``Reactants>>Products``) into three formats.

    Returns the same ``{original, canonical, display}`` dict as
    :func:`normalize_smiles_full`, but applied component-wise across the
    ``>>`` separator.
    """
    result: Dict[str, str] = {
        "original": reaction_smiles,
        "canonical": reaction_smiles,
        "display": reaction_smiles,
    }
    if not reaction_smiles or ">>" not in reaction_smiles:
        return result

    try:
        parts = reaction_smiles.split(">>")
        if len(parts) != 2:
            return result

        reactants_raw, products_raw = parts

        reactants_norm, reactants_disp = [], []
        products_norm, products_disp = [], []

        for r in reactants_raw.split("."):
            r = r.strip()
            if r:
                n = normalize_smiles_full(r)
                reactants_norm.append(n["canonical"])
                reactants_disp.append(n["display"])

        for p in products_raw.split("."):
            p = p.strip()
            if p:
                n = normalize_smiles_full(p)
                products_norm.append(n["canonical"])
                products_disp.append(n["display"])

        result["canonical"] = ".".join(reactants_norm) + ">>" + ".".join(products_norm)
        result["display"] = ".".join(reactants_disp) + ">>" + ".".join(products_disp)
    except Exception:
        pass

    return result


def normalize_smiles_list(smiles_list: List[str]) -> List[Dict[str, str]]:
    """Batch-normalize a list of SMILES strings.

    Returns a list of ``{original, canonical, display}`` dicts.
    """
    return [normalize_smiles_full(smi) for smi in smiles_list]


# ---------------------------------------------------------------------------
# Cache management
# ---------------------------------------------------------------------------

def clear_cache() -> None:
    """Clear the LRU caches (useful for testing or memory management)."""
    parse_smiles.cache_clear()
