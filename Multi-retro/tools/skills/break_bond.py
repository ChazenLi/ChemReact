"""
Break-Bond Skill
=================
Breaks a specified bond (by atom indices) in a molecule and generates
precursor fragments with retrosynthetic guidance for LLM decision-making.

v4: Three-layer fix
  - Layer 1: Hydrogen补齐 — RemoveBond 后根据键类型给断裂端补 H
  - Layer 2: 环键断裂支持 — 环内键断裂返回开环产物而非报错
  - Layer 3: 化学语义断键 — 识别酯/酰胺/醚等键型，生成化学正确的前体

Public API::

    skill = BreakBondSkill()
    result = skill.execute({"smiles": "...", "atom1_idx": 0, "atom2_idx": 3})
"""

from __future__ import annotations

import logging
from typing import Any, Dict, List, Optional, Tuple

from tools.skills.base import BaseSkill, SkillResult

logger = logging.getLogger(__name__)

try:
    from rdkit import Chem
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False


# ---------------------------------------------------------------------------
# Layer 1: Bond-breaking with hydrogen补齐
# ---------------------------------------------------------------------------

def _sanitize_fragment(frag) -> Optional[str]:
    """Sanitize a molecular fragment with multiple fallback strategies.

    Strategy 1: Full SanitizeMol (standard path)
    Strategy 2: Partial sanitize — skip KEKULIZE, then round-trip validate
    Strategy 3: Add explicit Hs to aromatic N in rings (fixes heterocycles)
    Strategy 4: Last resort — raw SMILES round-trip
    """
    import re as _re

    # Strategy 1: standard full sanitize
    try:
        frag_copy = Chem.RWMol(frag)
        Chem.SanitizeMol(frag_copy)
        smiles = Chem.MolToSmiles(frag_copy, canonical=True)
        if smiles:
            return smiles
    except Exception:
        pass

    # Strategy 2: partial sanitize — skip kekulize step
    try:
        frag_copy = Chem.RWMol(frag)
        Chem.SanitizeMol(
            frag_copy,
            sanitizeOps=Chem.SanitizeFlags.SANITIZE_ALL
            ^ Chem.SanitizeFlags.SANITIZE_KEKULIZE,
        )
        smiles = Chem.MolToSmiles(frag_copy, canonical=True)
        if smiles:
            check = Chem.MolFromSmiles(smiles)
            if check is not None:
                return Chem.MolToSmiles(check, canonical=True)
            fixed = _re.sub(r'(?<!\[)n(\d)', r'[nH]\1', smiles)
            check2 = Chem.MolFromSmiles(fixed)
            if check2 is not None:
                return Chem.MolToSmiles(check2, canonical=True)
    except Exception:
        pass

    # Strategy 3: add explicit Hs to aromatic N in rings (from legacy)
    try:
        frag_copy = Chem.RWMol(frag)
        for atom in frag_copy.GetAtoms():
            if (atom.GetSymbol() == "N"
                    and atom.GetIsAromatic()
                    and atom.GetTotalNumHs() == 0
                    and atom.GetFormalCharge() == 0
                    and atom.GetDegree() == 2):
                atom.SetNumExplicitHs(1)
                atom.SetNoImplicit(True)
        Chem.SanitizeMol(frag_copy)
        smiles = Chem.MolToSmiles(frag_copy, canonical=True)
        if smiles:
            return smiles
    except Exception:
        pass

    # Strategy 4: last resort — raw SMILES round-trip
    try:
        smiles = Chem.MolToSmiles(frag, canonical=True)
        if smiles:
            check = Chem.MolFromSmiles(smiles)
            if check is not None:
                return Chem.MolToSmiles(check, canonical=True)
    except Exception:
        pass

    return None


def _break_bond_and_get_fragments(
    mol, atom_idx1: int, atom_idx2: int
) -> List[str]:
    """Remove a bond and return sanitized fragment SMILES.

    Layer 1 fix: after RemoveBond, atoms with explicit H that had their
    valence reduced by the bond removal get their explicit H count bumped
    to maintain correct total valence.  For atoms using only implicit H
    (the common case), RDKit's sanitize automatically recalculates the
    correct implicit H count — no manual intervention needed.

    For double/triple bond breaks, the broken ends need extra explicit H
    to compensate for the higher bond order lost.

    Returns a list of canonical SMILES strings (may be 1 for ring bonds,
    >=2 for chain bonds, 0 on failure).
    """
    try:
        emol = Chem.RWMol(mol)
        bond = emol.GetBondBetweenAtoms(atom_idx1, atom_idx2)
        if bond is None:
            return []

        bt = bond.GetBondType()
        # Bond order: how many valence units the bond consumed
        bond_order = {
            Chem.BondType.SINGLE: 1,
            Chem.BondType.DOUBLE: 2,
            Chem.BondType.TRIPLE: 3,
        }.get(bt, 1)

        emol.RemoveBond(atom_idx1, atom_idx2)

        # ── Layer 1: explicit H compensation ──
        # For atoms that already had explicit H set, we must bump the count
        # because RDKit won't auto-recalculate for them.
        # For atoms with only implicit H AND bond_order > 1 (double/triple),
        # we also need to add explicit H because the implicit H recalculation
        # only accounts for single-bond valence changes.
        for idx in (atom_idx1, atom_idx2):
            atom = emol.GetAtomWithIdx(idx)
            if atom.GetIsAromatic():
                continue
            cur_explicit = atom.GetNumExplicitHs()
            if cur_explicit > 0 or atom.GetNoImplicit():
                # Atom uses explicit H tracking — must compensate manually
                atom.SetNumExplicitHs(cur_explicit + bond_order)
            elif bond_order > 1:
                # Implicit-H atom losing a multi-order bond: sanitize will
                # recalculate implicit H for single-bond loss, but not for
                # the extra valence from double/triple.  Add the difference.
                atom.SetNumExplicitHs(bond_order - 1)

        frags = Chem.GetMolFrags(emol.GetMol(), asMols=True, sanitizeFrags=False)
        result: List[str] = []
        for frag in frags:
            smiles = _sanitize_fragment(frag)
            if smiles:
                result.append(smiles)
        return result
    except Exception:
        return []


# ---------------------------------------------------------------------------
# Bond-type detection helpers
# ---------------------------------------------------------------------------

def _is_ester_bond(mol, bond) -> bool:
    """Return True if *bond* is the C-O single bond in an ester C(=O)-O-C."""
    a1 = bond.GetBeginAtom()
    a2 = bond.GetEndAtom()

    def _check(carbon, oxygen):
        if carbon.GetSymbol() != "C" or oxygen.GetSymbol() != "O":
            return False
        has_carbonyl = False
        for nb in carbon.GetNeighbors():
            if nb.GetIdx() == oxygen.GetIdx():
                continue
            if nb.GetSymbol() == "O":
                b = mol.GetBondBetweenAtoms(carbon.GetIdx(), nb.GetIdx())
                if b and str(b.GetBondType()).endswith("DOUBLE"):
                    has_carbonyl = True
                    break
        if not has_carbonyl:
            return False
        for nb in oxygen.GetNeighbors():
            if nb.GetIdx() == carbon.GetIdx():
                continue
            if nb.GetSymbol() == "C":
                return True
        return False

    return _check(a1, a2) or _check(a2, a1)


def _is_amide_bond(mol, bond) -> bool:
    """Return True if *bond* is the C-N single bond in an amide C(=O)-N."""
    a1 = bond.GetBeginAtom()
    a2 = bond.GetEndAtom()

    def _check(carbon, nitrogen):
        if carbon.GetSymbol() != "C" or nitrogen.GetSymbol() != "N":
            return False
        for nb in carbon.GetNeighbors():
            if nb.GetIdx() == nitrogen.GetIdx():
                continue
            if nb.GetSymbol() == "O":
                b = mol.GetBondBetweenAtoms(carbon.GetIdx(), nb.GetIdx())
                if b and str(b.GetBondType()).endswith("DOUBLE"):
                    return True
        return False

    return _check(a1, a2) or _check(a2, a1)


def _is_alpha_carbonyl_cc(mol, bond) -> bool:
    """Return True if *bond* is a C-C bond alpha to a carbonyl (C=O)."""
    a1 = bond.GetBeginAtom()
    a2 = bond.GetEndAtom()
    if a1.GetSymbol() != "C" or a2.GetSymbol() != "C":
        return False
    for atom in (a1, a2):
        for nb in atom.GetNeighbors():
            if nb.GetSymbol() == "O":
                b = mol.GetBondBetweenAtoms(atom.GetIdx(), nb.GetIdx())
                if b and str(b.GetBondType()).endswith("DOUBLE"):
                    return True
    return False


def _is_sulfonamide_bond(mol, bond) -> bool:
    """Return True if *bond* is the S-N bond in a sulfonamide -SO2-N-."""
    a1 = bond.GetBeginAtom()
    a2 = bond.GetEndAtom()

    def _check(s_atom, n_atom):
        if s_atom.GetSymbol() != "S" or n_atom.GetSymbol() != "N":
            return False
        double_o = 0
        for nb in s_atom.GetNeighbors():
            if nb.GetIdx() == n_atom.GetIdx():
                continue
            if nb.GetSymbol() == "O":
                b = mol.GetBondBetweenAtoms(s_atom.GetIdx(), nb.GetIdx())
                if b and str(b.GetBondType()).endswith("DOUBLE"):
                    double_o += 1
        return double_o >= 2

    return _check(a1, a2) or _check(a2, a1)


def _is_sulfonyl_ester_bond(mol, bond) -> bool:
    """Return True if *bond* is the S-O single bond in a sulfonyl ester -SO2-O-R."""
    a1 = bond.GetBeginAtom()
    a2 = bond.GetEndAtom()
    if str(bond.GetBondType()).endswith("DOUBLE"):
        return False

    def _check(s_atom, o_atom):
        if s_atom.GetSymbol() != "S" or o_atom.GetSymbol() != "O":
            return False
        double_o = 0
        for nb in s_atom.GetNeighbors():
            if nb.GetIdx() == o_atom.GetIdx():
                continue
            if nb.GetSymbol() == "O":
                b = mol.GetBondBetweenAtoms(s_atom.GetIdx(), nb.GetIdx())
                if b and str(b.GetBondType()).endswith("DOUBLE"):
                    double_o += 1
        if double_o < 2:
            return False
        for nb in o_atom.GetNeighbors():
            if nb.GetIdx() == s_atom.GetIdx():
                continue
            if nb.GetSymbol() != "H":
                return True
        return o_atom.GetDegree() > 1

    return _check(a1, a2) or _check(a2, a1)


def _is_ether_bond(mol, bond) -> bool:
    """Return True if *bond* is a C-O single bond in an ether C-O-C (not ester)."""
    if _is_ester_bond(mol, bond):
        return False
    a1 = bond.GetBeginAtom()
    a2 = bond.GetEndAtom()
    if str(bond.GetBondType()).endswith("DOUBLE"):
        return False

    def _check(c_atom, o_atom):
        if c_atom.GetSymbol() != "C" or o_atom.GetSymbol() != "O":
            return False
        # O must connect to another C (not just H)
        for nb in o_atom.GetNeighbors():
            if nb.GetIdx() == c_atom.GetIdx():
                continue
            if nb.GetSymbol() == "C":
                return True
        return False

    return _check(a1, a2) or _check(a2, a1)


def _is_diarylamine_bond(mol, bond) -> bool:
    """Return True if *bond* is a C-N bond in a diarylamine Ar-NH-Ar (non-amide)."""
    if _is_amide_bond(mol, bond):
        return False
    a1 = bond.GetBeginAtom()
    a2 = bond.GetEndAtom()

    def _check(c_atom, n_atom):
        if c_atom.GetSymbol() != "C" or n_atom.GetSymbol() != "N":
            return False
        if not c_atom.GetIsAromatic():
            return False
        for nb in n_atom.GetNeighbors():
            if nb.GetIdx() == c_atom.GetIdx():
                continue
            if nb.GetSymbol() == "C" and nb.GetIsAromatic():
                return True
        return False

    return _check(a1, a2) or _check(a2, a1)


def _is_sulfonyl_double_o(mol, bond) -> bool:
    """Return True if *bond* is an internal S=O double bond in a sulfonyl group."""
    a1 = bond.GetBeginAtom()
    a2 = bond.GetEndAtom()
    if not str(bond.GetBondType()).endswith("DOUBLE"):
        return False
    s_atom = a1 if a1.GetSymbol() == "S" else (a2 if a2.GetSymbol() == "S" else None)
    o_atom = a1 if a1.GetSymbol() == "O" else (a2 if a2.GetSymbol() == "O" else None)
    if s_atom is None or o_atom is None:
        return False
    return s_atom.GetDegree() >= 3


def _is_phosphoryl_double_o(mol, bond) -> bool:
    """Return True if *bond* is an internal P=O double bond in a phosphoryl group."""
    a1 = bond.GetBeginAtom()
    a2 = bond.GetEndAtom()
    if not str(bond.GetBondType()).endswith("DOUBLE"):
        return False
    p_atom = a1 if a1.GetSymbol() == "P" else (a2 if a2.GetSymbol() == "P" else None)
    o_atom = a1 if a1.GetSymbol() == "O" else (a2 if a2.GetSymbol() == "O" else None)
    return p_atom is not None and o_atom is not None


# ---------------------------------------------------------------------------
# Layer 3: SMARTS-driven semantic bond-breaking
# ---------------------------------------------------------------------------
# Reuses the 35+ DISCONNECTION_RULES from the disconnector module.
# Each rule carries a SMARTS pattern, bond_to_break indices, and metadata.
# We match the broken bond against these rules and apply fragment corrections
# based on a SEMANTIC_CORRECTIONS table keyed by rule strategy.

try:
    from tools.legacy.planner.disconnector import DISCONNECTION_RULES as _DISC_RULES
except ImportError:
    try:
        from multiretro.core.planner.disconnector import DISCONNECTION_RULES as _DISC_RULES
    except ImportError:
        _DISC_RULES = []


def _replace_h_with_atom_on_fragment(
    frag_smiles: str,
    target_symbol: str,
    new_atom_num: int,
    new_bond_type=None,
    prefer_carbonyl_neighbor: bool = False,
) -> Optional[str]:
    """In *frag_smiles*, find an atom with symbol *target_symbol* that has
    at least one H, remove one H and attach *new_atom_num* via *new_bond_type*.

    Used to convert e.g. R-C(=O)H → R-C(=O)OH  or  R-C(=O)H → R-C(=O)Cl.
    """
    if new_bond_type is None:
        new_bond_type = Chem.BondType.SINGLE
    frag_mol = Chem.MolFromSmiles(frag_smiles)
    if frag_mol is None:
        return None

    rw = Chem.RWMol(frag_mol)
    target_idx = None
    for a in rw.GetAtoms():
        if a.GetSymbol() != target_symbol or a.GetTotalNumHs() == 0:
            continue
        if prefer_carbonyl_neighbor:
            has_co = any(
                nb.GetSymbol() == "O"
                and str(rw.GetBondBetweenAtoms(a.GetIdx(), nb.GetIdx()).GetBondType()).endswith("DOUBLE")
                for nb in a.GetNeighbors() if nb.GetSymbol() == "O"
            )
            if has_co:
                target_idx = a.GetIdx()
                break
        if target_idx is None:
            target_idx = a.GetIdx()

    if target_idx is None:
        return None

    ta = rw.GetAtomWithIdx(target_idx)
    explicit_h = ta.GetNumExplicitHs()
    if explicit_h > 0:
        ta.SetNumExplicitHs(explicit_h - 1)
    new_idx = rw.AddAtom(Chem.Atom(new_atom_num))
    rw.AddBond(target_idx, new_idx, new_bond_type)
    try:
        Chem.SanitizeMol(rw)
        return Chem.MolToSmiles(rw, canonical=True)
    except Exception:
        return None


def _find_fragment_with_feature(
    fragments: List[str],
    symbol: str,
    require_double_o: bool = False,
    require_carbonyl: bool = False,
) -> Optional[int]:
    """Find fragment index containing an atom matching criteria."""
    for i, smi in enumerate(fragments):
        m = Chem.MolFromSmiles(smi)
        if m is None:
            continue
        for a in m.GetAtoms():
            if a.GetSymbol() != symbol:
                continue
            if require_double_o:
                dbl_o = sum(
                    1 for nb in a.GetNeighbors()
                    if nb.GetSymbol() == "O"
                    and str(m.GetBondBetweenAtoms(a.GetIdx(), nb.GetIdx()).GetBondType()).endswith("DOUBLE")
                )
                if dbl_o >= 2:
                    return i
                continue
            if require_carbonyl:
                for nb in a.GetNeighbors():
                    if nb.GetSymbol() == "O":
                        b = m.GetBondBetweenAtoms(a.GetIdx(), nb.GetIdx())
                        if b and str(b.GetBondType()).endswith("DOUBLE"):
                            if a.GetTotalNumHs() > 0:
                                return i
                continue
            if a.GetTotalNumHs() > 0:
                return i
    return None


def _find_aromatic_fragment(fragments: List[str]) -> Optional[int]:
    """Return index of the fragment containing aromatic atoms."""
    for i, smi in enumerate(fragments):
        m = Chem.MolFromSmiles(smi)
        if m and any(a.GetIsAromatic() for a in m.GetAtoms()):
            return i
    return None


def _smaller_fragment_idx(fragments: List[str]) -> int:
    """Return index of the smaller fragment by heavy atom count."""
    sizes = []
    for smi in fragments:
        m = Chem.MolFromSmiles(smi)
        sizes.append(m.GetNumHeavyAtoms() if m else 999)
    return sizes.index(min(sizes))


# ── Correction functions keyed by strategy ──
# Each returns (corrected_fragments, reaction_type, note) or None.

def _correct_acylation(rule, fragments):
    """Ester/amide/acyl: C=O end → add OH (acid) or Cl (acyl chloride)."""
    c_idx = _find_fragment_with_feature(fragments, "C", require_carbonyl=True)
    if c_idx is None:
        return None
    new_frags = list(fragments)
    corrected = _replace_h_with_atom_on_fragment(
        new_frags[c_idx], "C", 8, prefer_carbonyl_neighbor=True,  # O → acid
    )
    if corrected is None:
        return None
    new_frags[c_idx] = corrected
    return (new_frags, rule["reaction_class"], f"{rule['rationale']}; C=O end → C(=O)OH")


def _correct_cross_coupling_halide(rule, fragments):
    """Cross-coupling: aromatic end → add Br (aryl halide)."""
    ar_idx = _find_aromatic_fragment(fragments)
    if ar_idx is None:
        return None
    new_frags = list(fragments)
    corrected = _replace_h_with_atom_on_fragment(new_frags[ar_idx], "C", 35)  # Br
    if corrected is None:
        return None
    new_frags[ar_idx] = corrected
    func = rule.get("functionalize", "")
    note = f"{rule['rationale']}; aryl end → Ar-Br"
    if func == "suzuki":
        # Non-aryl end → boronic acid (informational, hard to do via atom add)
        note += "; other end ideally → Ar-B(OH)2"
    return (new_frags, rule["reaction_class"], note)


def _correct_sulfonyl(rule, fragments):
    """Sulfonamide/sulfonyl ester: S end → add Cl (sulfonyl chloride)."""
    s_idx = _find_fragment_with_feature(fragments, "S", require_double_o=True)
    if s_idx is None:
        return None
    new_frags = list(fragments)
    corrected = _replace_h_with_atom_on_fragment(new_frags[s_idx], "S", 17)  # Cl
    if corrected is None:
        return None
    new_frags[s_idx] = corrected
    return (new_frags, rule["reaction_class"], f"{rule['rationale']}; S end → SO2Cl")


def _correct_ether_williamson(rule, fragments):
    """Ether: smaller end → add Br (alkyl halide)."""
    if len(fragments) < 2:
        return None
    halide_idx = _smaller_fragment_idx(fragments)
    new_frags = list(fragments)
    corrected = _replace_h_with_atom_on_fragment(new_frags[halide_idx], "C", 35)  # Br
    if corrected is not None:
        new_frags[halide_idx] = corrected
    return (new_frags, rule["reaction_class"], f"{rule['rationale']}; smaller end → R-Br")


def _correct_grignard(rule, fragments):
    """Grignard: C-C(OH) break → one end becomes aldehyde/ketone (already
    has C=O after break), other end → RMgBr (informational)."""
    # The fragment with OH is the product side; the other is the Grignard
    # After topological break, both are valid — just add note
    return (list(fragments), rule["reaction_class"], rule["rationale"])


def _correct_friedel_crafts(rule, fragments):
    """Friedel-Crafts: non-aromatic end → add Cl (alkyl/acyl chloride)."""
    if len(fragments) < 2:
        return None
    ar_idx = _find_aromatic_fragment(fragments)
    if ar_idx is None:
        return None
    other_idx = 1 - ar_idx if len(fragments) == 2 else None
    if other_idx is None:
        return None
    new_frags = list(fragments)
    corrected = _replace_h_with_atom_on_fragment(new_frags[other_idx], "C", 17)  # Cl
    if corrected is None:
        return None
    new_frags[other_idx] = corrected
    return (new_frags, rule["reaction_class"], f"{rule['rationale']}; alkyl/acyl end → R-Cl")


def _correct_protection(rule, fragments):
    """Protecting group removal: just return fragments as-is with metadata."""
    return (list(fragments), rule["reaction_class"], rule["rationale"])


def _correct_silyl(rule, fragments):
    """Silyl ether: O-Si break → alcohol + silyl halide."""
    return (list(fragments), rule["reaction_class"], rule["rationale"])


# Strategy → correction function mapping
_STRATEGY_CORRECTIONS: Dict[str, Any] = {
    "ester_hydrolysis": _correct_acylation,
    "amide_hydrolysis": _correct_acylation,
    "acyl_substitution": _correct_acylation,
    "aldol": _correct_acylation,
    "claisen": _correct_acylation,
    "cross_coupling": _correct_cross_coupling_halide,
    "grignard": _correct_grignard,
    "wittig": _correct_grignard,  # similar: fragments are already reasonable
    "electrophilic_aromatic": _correct_friedel_crafts,
    "ether_cleavage": _correct_ether_williamson,
    "substitution": _correct_ether_williamson,
    "reductive_amination": _correct_grignard,
    "protection": _correct_protection,
    "ring_opening": _correct_protection,
    "conjugate_addition": _correct_grignard,
    "metathesis": _correct_grignard,
    "cycloaddition": _correct_grignard,
    "heterocycle_formation": _correct_protection,
}


def _match_smarts_rule(
    mol, atom_idx1: int, atom_idx2: int,
) -> Optional[Dict[str, Any]]:
    """Find the highest-specificity SMARTS rule whose bond_to_break matches
    the given atom pair.  Returns the matched rule dict or None."""
    if not RDKIT_AVAILABLE or not _DISC_RULES:
        return None

    best_rule = None
    best_score = -1  # (specificity * 10 + confidence)

    for rule in _DISC_RULES:
        if rule.get("is_redox"):
            continue  # redox rules don't break bonds
        bond_indices = rule.get("bond_to_break")
        if bond_indices is None:
            continue
        pattern_str = rule.get("pattern", "")
        try:
            pattern = Chem.MolFromSmarts(pattern_str)
        except Exception:
            continue
        if pattern is None:
            continue

        matches = mol.GetSubstructMatches(pattern)
        for match in matches:
            try:
                rule_a1 = match[bond_indices[0]]
                rule_a2 = match[bond_indices[1]]
            except IndexError:
                continue
            if {rule_a1, rule_a2} == {atom_idx1, atom_idx2}:
                score = rule.get("specificity", 1) * 10 + rule.get("confidence", 0.5)
                if score > best_score:
                    best_score = score
                    best_rule = rule

    return best_rule


def _semantic_transform_fragments(
    mol,
    bond,
    atom_idx1: int,
    atom_idx2: int,
    raw_fragments: List[str],
) -> Optional[Tuple[List[str], str, str]]:
    """Apply retrosynthetic semantic transforms to raw fragments.

    Uses the 35+ SMARTS DISCONNECTION_RULES to identify the reaction type,
    then applies the appropriate fragment correction (add leaving groups,
    convert to proper precursors).

    Falls back to the hardcoded _is_*_bond detectors if no SMARTS rule
    matches (e.g. when the disconnector module is unavailable).

    Returns:
        (corrected_fragments, reaction_type, correction_note) or None.
    """
    if len(raw_fragments) < 2:
        return None

    # ── Primary path: SMARTS rule matching ──
    try:
        rule = _match_smarts_rule(mol, atom_idx1, atom_idx2)
        if rule is not None:
            strategy = rule.get("strategy", "")
            correction_fn = _STRATEGY_CORRECTIONS.get(strategy)
            if correction_fn is not None:
                result = correction_fn(rule, raw_fragments)
                if result is not None:
                    return result
    except Exception:
        pass

    # ── Fallback: hardcoded bond-type detectors ──
    try:
        if _is_ester_bond(mol, bond):
            c_idx = _find_fragment_with_feature(raw_fragments, "C", require_carbonyl=True)
            if c_idx is not None:
                new_frags = list(raw_fragments)
                corrected = _replace_h_with_atom_on_fragment(
                    new_frags[c_idx], "C", 8, prefer_carbonyl_neighbor=True,
                )
                if corrected:
                    new_frags[c_idx] = corrected
                    return (new_frags, "Fischer Esterification (reverse)",
                            "Ester retro: C=O end → C(=O)OH; O-end → R-OH")

        if _is_amide_bond(mol, bond):
            c_idx = _find_fragment_with_feature(raw_fragments, "C", require_carbonyl=True)
            if c_idx is not None:
                new_frags = list(raw_fragments)
                corrected = _replace_h_with_atom_on_fragment(
                    new_frags[c_idx], "C", 8, prefer_carbonyl_neighbor=True,
                )
                if corrected:
                    new_frags[c_idx] = corrected
                    return (new_frags, "Amide Coupling (reverse)",
                            "Amide retro: C=O end → C(=O)OH; N-end → amine")

        if _is_sulfonamide_bond(mol, bond):
            s_idx = _find_fragment_with_feature(raw_fragments, "S", require_double_o=True)
            if s_idx is not None:
                new_frags = list(raw_fragments)
                corrected = _replace_h_with_atom_on_fragment(new_frags[s_idx], "S", 17)
                if corrected:
                    new_frags[s_idx] = corrected
                    return (new_frags, "Sulfonamide Formation (reverse)",
                            "Sulfonamide retro: S end → SO2Cl; N-end → amine")

        if _is_sulfonyl_ester_bond(mol, bond):
            s_idx = _find_fragment_with_feature(raw_fragments, "S", require_double_o=True)
            if s_idx is not None:
                new_frags = list(raw_fragments)
                corrected = _replace_h_with_atom_on_fragment(new_frags[s_idx], "S", 17)
                if corrected:
                    new_frags[s_idx] = corrected
                    return (new_frags, "Sulfonyl Ester Formation (reverse)",
                            "Sulfonyl ester retro: S end → SO2Cl; O-end → R-OH")

        if _is_ether_bond(mol, bond):
            halide_idx = _smaller_fragment_idx(raw_fragments)
            new_frags = list(raw_fragments)
            corrected = _replace_h_with_atom_on_fragment(new_frags[halide_idx], "C", 35)
            if corrected:
                new_frags[halide_idx] = corrected
            return (new_frags, "Williamson Ether Synthesis (reverse)",
                    "Ether retro: smaller end → R-Br; larger end → R-OH")
    except Exception:
        pass

    return None


# ---------------------------------------------------------------------------
# Fragment profiling (detailed)
# ---------------------------------------------------------------------------

def _build_fragment_profiles(fragments: List[str]) -> List[Dict[str, Any]]:
    """Build detailed chemical profiles for each fragment SMILES."""
    if not RDKIT_AVAILABLE:
        return [{"index": i, "smiles": f, "valid": False}
                for i, f in enumerate(fragments)]

    profiles: List[Dict[str, Any]] = []
    for i, frag in enumerate(fragments):
        frag_mol = Chem.MolFromSmiles(frag)
        profile: Dict[str, Any] = {
            "index": i, "smiles": frag, "valid": frag_mol is not None,
        }
        if frag_mol is None:
            profiles.append(profile)
            continue

        ring_n = sum(1 for a in frag_mol.GetAtoms()
                     if a.GetSymbol() == "N" and a.IsInRing())
        total_n = sum(1 for a in frag_mol.GetAtoms()
                      if a.GetSymbol() == "N")
        has_nh = any(a.GetSymbol() == "N" and a.GetTotalNumHs() > 0
                     for a in frag_mol.GetAtoms())
        has_aromatic = any(a.GetIsAromatic() for a in frag_mol.GetAtoms())
        num_heavy = frag_mol.GetNumHeavyAtoms()

        fg_tags: List[str] = []
        for a in frag_mol.GetAtoms():
            sym = a.GetSymbol()
            if sym == "O" and a.GetTotalNumHs() > 0 and not a.IsInRing():
                fg_tags.append("OH")
            if sym == "N" and a.GetTotalNumHs() >= 2 and not a.IsInRing():
                fg_tags.append("NH2")
            if sym == "N" and a.GetTotalNumHs() == 1 and not a.IsInRing():
                fg_tags.append("NH")
            if sym == "S":
                dbl_o = sum(
                    1 for nb in a.GetNeighbors()
                    if nb.GetSymbol() == "O"
                    and str(frag_mol.GetBondBetweenAtoms(
                        a.GetIdx(), nb.GetIdx()).GetBondType()).endswith("DOUBLE")
                )
                if dbl_o >= 2:
                    fg_tags.append("SO2")
            if sym in ("F", "Cl", "Br", "I"):
                fg_tags.append(sym)
        for b in frag_mol.GetBonds():
            ba, bb = b.GetBeginAtom(), b.GetEndAtom()
            syms = {ba.GetSymbol(), bb.GetSymbol()}
            bt = str(b.GetBondType())
            if syms == {"C", "O"} and bt.endswith("DOUBLE"):
                fg_tags.append("C=O")
            if syms == {"C", "N"} and bt.endswith("TRIPLE"):
                fg_tags.append("CN_triple")

        profile.update({
            "ring_N_count": ring_n,
            "total_N_count": total_n,
            "has_NH": has_nh,
            "has_aromatic_ring": has_aromatic,
            "heavy_atom_count": num_heavy,
            "functional_groups": sorted(set(fg_tags)),
        })
        profiles.append(profile)
    return profiles


# ---------------------------------------------------------------------------
# Retro-guide builder
# ---------------------------------------------------------------------------

def _build_retro_guide(
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
    """Build a retro analysis guide dict for LLM decision-making."""
    guide: Dict[str, Any] = {
        "bond_type": bond_semantic,
        "retro_pattern": retro_pattern,
        "broken_atoms": {
            "atom1": {
                "idx": atom1_idx,
                "symbol": mol.GetAtomWithIdx(atom1_idx).GetSymbol(),
                "in_ring": mol.GetAtomWithIdx(atom1_idx).IsInRing(),
                "is_aromatic": mol.GetAtomWithIdx(atom1_idx).GetIsAromatic(),
            },
            "atom2": {
                "idx": atom2_idx,
                "symbol": mol.GetAtomWithIdx(atom2_idx).GetSymbol(),
                "in_ring": mol.GetAtomWithIdx(atom2_idx).IsInRing(),
                "is_aromatic": mol.GetAtomWithIdx(atom2_idx).GetIsAromatic(),
            },
        },
        "fragment_profiles": _build_fragment_profiles(fragments),
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


# ---------------------------------------------------------------------------
# BreakBondSkill
# ---------------------------------------------------------------------------

class BreakBondSkill(BaseSkill):
    """Break a bond by atom indices to generate precursor SMILES automatically.

    Args (dict):
        smiles     – target molecule SMILES
        atom1_idx  – index of the first atom in the bond to break
        atom2_idx  – index of the second atom in the bond to break

    Returns ``SkillResult.to_dict()`` with fragments, reaction SMILES,
    and an optional ``retro_analysis_guide`` when a known reaction pattern
    is detected.

    v4 enhancements:
      - Layer 1: H补齐 on broken ends
      - Layer 2: ring bond breaking returns opened-ring product
      - Layer 3: semantic transforms produce correct precursors
    """

    name = "break_bond"
    description = "Break a bond by atom indices to generate precursor SMILES automatically"

    def execute(self, args: Any) -> Dict[str, Any]:  # noqa: C901
        if not RDKIT_AVAILABLE:
            return SkillResult(success=False, error="RDKit not available").to_dict()

        # ── Parse arguments ──
        if isinstance(args, dict):
            smiles = args.get("smiles", "")
            atom1_idx = int(args.get("atom1_idx", 0))
            atom2_idx = int(args.get("atom2_idx", 0))
        else:
            smiles = getattr(args, "smiles", "")
            atom1_idx = int(getattr(args, "atom1_idx", 0))
            atom2_idx = int(getattr(args, "atom2_idx", 0))

        # ── Parse & canonicalize ──
        from tools.chem.mol_parser import parse_smiles, normalize_smiles

        mol = parse_smiles(smiles)
        if mol is None:
            return SkillResult(
                success=False, error=f"Invalid SMILES: {smiles}"
            ).to_dict()

        # canonical 仅用于输出显示，不重新 parse mol。
        # 重新 parse 会重排原子索引，导致用户传入的 atom_idx 对不上。
        canonical = normalize_smiles(smiles)

        # ── Validate indices ──
        n_atoms = mol.GetNumAtoms()
        if atom1_idx < 0 or atom1_idx >= n_atoms or atom2_idx < 0 or atom2_idx >= n_atoms:
            return SkillResult(
                success=False,
                error=(
                    f"Atom index out of range. Molecule has {n_atoms} atoms "
                    f"(0-{n_atoms - 1}). Got atom1_idx={atom1_idx}, atom2_idx={atom2_idx}"
                ),
            ).to_dict()

        bond = mol.GetBondBetweenAtoms(atom1_idx, atom2_idx)
        if bond is None:
            return SkillResult(
                success=False,
                error=(
                    f"No bond between atom {atom1_idx} "
                    f"({mol.GetAtomWithIdx(atom1_idx).GetSymbol()}) and atom "
                    f"{atom2_idx} ({mol.GetAtomWithIdx(atom2_idx).GetSymbol()})"
                ),
            ).to_dict()

        # ── Bond metadata ──
        atom1_sym = mol.GetAtomWithIdx(atom1_idx).GetSymbol()
        atom2_sym = mol.GetAtomWithIdx(atom2_idx).GetSymbol()
        bond_type = str(bond.GetBondType()).replace(
            "rdkit.Chem.rdchem.BondType.", ""
        )
        is_aromatic = bond.GetIsAromatic()
        in_ring = bond.IsInRing()

        if is_aromatic:
            return SkillResult(
                success=False,
                error=(
                    f"Bond {atom1_idx}-{atom2_idx} is aromatic. "
                    "Breaking aromatic bonds is chemically unreasonable."
                ),
            ).to_dict()

        # ── Break the bond (Layer 1: with H补齐) ──
        fragments = _break_bond_and_get_fragments(mol, atom1_idx, atom2_idx)

        if len(fragments) == 0:
            return SkillResult(
                success=False,
                error="Bond break produced no valid fragments.",
                data={"fragments": []},
            ).to_dict()

        # ── Layer 2: ring bond support ──
        ring_opened = False
        if len(fragments) == 1 and in_ring:
            # Ring bond break → molecule stays connected (ring opened).
            # This is chemically valid (e.g. lactone hydrolysis, epoxide
            # opening).  Return the single opened-ring product.
            ring_opened = True
        elif len(fragments) < 2:
            return SkillResult(
                success=False,
                error=(
                    "Bond break did not produce valid fragments. "
                    "This may happen with ring bonds that don't split the molecule."
                ),
                data={"fragments": fragments},
            ).to_dict()

        # ── Layer 3: semantic transform ──
        semantic_result = None
        if len(fragments) >= 2:
            semantic_result = _semantic_transform_fragments(
                mol, bond, atom1_idx, atom2_idx, fragments,
            )

        # Determine which fragments to use for the reaction SMILES
        if semantic_result is not None:
            corrected_fragments, reaction_type, correction_note = semantic_result
            display_fragments = corrected_fragments
        else:
            display_fragments = fragments
            reaction_type = None
            correction_note = None

        if ring_opened:
            reaction_smiles = display_fragments[0] + ">>" + canonical
        else:
            reaction_smiles = ".".join(display_fragments) + ">>" + canonical

        # ── Retro-analysis guidance ──
        retro_guide = self._detect_retro_pattern(
            mol, bond, atom1_idx, atom2_idx, atom1_sym, atom2_sym,
            bond_type, display_fragments,
        )

        # ── Build result ──
        data: Dict[str, Any] = {
            "target_smiles": canonical,
            "bond_broken": {
                "atom1_idx": atom1_idx,
                "atom1_symbol": atom1_sym,
                "atom2_idx": atom2_idx,
                "atom2_symbol": atom2_sym,
                "bond_type": bond_type,
                "in_ring": in_ring,
            },
            "fragments": display_fragments,
            "num_fragments": len(display_fragments),
            "reaction_smiles": reaction_smiles,
            "next_step": (
                "Use validate_reaction to verify this reaction: "
                f'--reaction_smiles "{reaction_smiles}"'
            ),
        }

        if ring_opened:
            data["ring_opened"] = True
            data["raw_fragments"] = fragments
            data["next_step"] = (
                "Ring bond broken — molecule opened (not split). "
                "The single fragment is the ring-opened precursor. "
                f'Validate: --reaction_smiles "{reaction_smiles}"'
            )

        # Attach raw (pre-semantic) fragments when semantic correction applied
        if semantic_result is not None:
            data["raw_fragments"] = fragments
            data["semantic_correction"] = {
                "reaction_type": reaction_type,
                "note": correction_note,
                "corrected_fragments": corrected_fragments,
            }
            data["next_step"] = (
                f"[Semantic] {correction_note} "
                f'Validate: --reaction_smiles "{reaction_smiles}"'
            )

        if retro_guide:
            data["retro_analysis_guide"] = retro_guide
            if semantic_result is None:
                data["next_step"] = (
                    f"[Guided] Detected {retro_guide['retro_pattern']}. "
                    f"Refer to retro_analysis_guide (fragment_profiles & retro_hint) "
                    f"to build correct precursor SMILES, then validate with validate_reaction."
                )

        return SkillResult(success=True, data=data).to_dict()

    # ------------------------------------------------------------------
    # Internal: pattern detection dispatcher
    # ------------------------------------------------------------------

    def _detect_retro_pattern(
        self,
        mol,
        bond,
        atom1_idx: int,
        atom2_idx: int,
        atom1_sym: str,
        atom2_sym: str,
        bond_type: str,
        fragments: List[str],
    ) -> Optional[Dict[str, Any]]:
        """Detect known retrosynthetic patterns and return a guide dict."""
        try:
            return self._detect_retro_pattern_inner(
                mol, bond, atom1_idx, atom2_idx, atom1_sym, atom2_sym,
                bond_type, fragments,
            )
        except Exception:
            return None

    def _detect_retro_pattern_inner(  # noqa: C901
        self,
        mol,
        bond,
        atom1_idx: int,
        atom2_idx: int,
        atom1_sym: str,
        atom2_sym: str,
        bond_type: str,
        fragments: List[str],
    ) -> Optional[Dict[str, Any]]:
        """Core pattern matching — separated so the outer method can catch errors."""

        if _is_ester_bond(mol, bond):
            return _build_retro_guide(
                mol, bond, atom1_idx, atom2_idx, fragments,
                bond_semantic="C(=O)-O (ester bond)",
                retro_pattern="ester hydrolysis",
                retro_hint=(
                    "Ester retrosynthesis: C=O end → carboxylic acid C(=O)OH; "
                    "O-C end → alcohol R-OH. Fischer esterification reverse. "
                    "Reagents: H2SO4 (cat.), reflux. Byproduct: H2O."
                ),
            )

        if _is_amide_bond(mol, bond):
            return _build_retro_guide(
                mol, bond, atom1_idx, atom2_idx, fragments,
                bond_semantic="C(=O)-N (amide bond)",
                retro_pattern="amide coupling",
                retro_hint=(
                    "Amide retrosynthesis: after C(=O)-N cleavage, the C=O end should "
                    "become a carboxylic acid C(=O)OH or acyl chloride C(=O)Cl; the N end "
                    "stays as an amine/NH. Common methods: EDC/DCC coupling, "
                    "acyl chloride + amine."
                ),
            )

        if _is_sulfonamide_bond(mol, bond):
            return _build_retro_guide(
                mol, bond, atom1_idx, atom2_idx, fragments,
                bond_semantic="S(=O)2-N (sulfonamide bond)",
                retro_pattern="sulfonamide formation",
                retro_hint=(
                    "Sulfonamide retrosynthesis: after S-N cleavage, the S end should "
                    "become a sulfonyl chloride -SO2Cl; the N end stays as an amine/NH. "
                    "Basic conditions (Et3N/pyridine)."
                ),
            )

        if _is_sulfonyl_ester_bond(mol, bond):
            return _build_retro_guide(
                mol, bond, atom1_idx, atom2_idx, fragments,
                bond_semantic="S(=O)2-O (sulfonyl ester bond)",
                retro_pattern="sulfonyl ester formation",
                retro_hint=(
                    "Sulfonyl ester retrosynthesis: after S-O cleavage, the S end becomes "
                    "a sulfonyl chloride -SO2Cl; the O end stays as an alcohol -OH. "
                    "Basic conditions."
                ),
            )

        if _is_ether_bond(mol, bond):
            return _build_retro_guide(
                mol, bond, atom1_idx, atom2_idx, fragments,
                bond_semantic="C-O-C (ether bond)",
                retro_pattern="ether cleavage (Williamson)",
                retro_hint=(
                    "Ether retrosynthesis (Williamson): one end → alcohol R-OH "
                    "(nucleophile), other end → alkyl halide R-X (electrophile). "
                    "Base (NaH/K2CO3), solvent (DMF/THF)."
                ),
            )

        if _is_diarylamine_bond(mol, bond):
            return _build_retro_guide(
                mol, bond, atom1_idx, atom2_idx, fragments,
                bond_semantic="Ar-N-Ar (diarylamine bond)",
                retro_pattern="diarylamine formation",
                retro_hint=(
                    "Diarylamine retrosynthesis: after Ar-N cleavage, one end becomes "
                    "an aryl amine Ar-NH2, the other an aryl halide Ar-X. Possible "
                    "methods: Buchwald-Hartwig (Pd-catalyzed, Ar-Br/I + Ar-NH2), "
                    "SNAr (activated aryl halide F/Cl + EWG)."
                ),
            )

        # Heck-type: Ar-C(vinyl) single bond
        if (
            {atom1_sym, atom2_sym} == {"C", "C"}
            and bond_type == "SINGLE"
            and (
                mol.GetAtomWithIdx(atom1_idx).GetIsAromatic()
                or mol.GetAtomWithIdx(atom2_idx).GetIsAromatic()
            )
            and not (
                mol.GetAtomWithIdx(atom1_idx).GetIsAromatic()
                and mol.GetAtomWithIdx(atom2_idx).GetIsAromatic()
            )
        ):
            aryl_idx = (
                atom1_idx
                if mol.GetAtomWithIdx(atom1_idx).GetIsAromatic()
                else atom2_idx
            )
            other_idx = atom2_idx if aryl_idx == atom1_idx else atom1_idx
            other_atom = mol.GetAtomWithIdx(other_idx)
            has_vinyl = any(
                nb.GetSymbol() == "C"
                and nb.GetIdx() != aryl_idx
                and mol.GetBondBetweenAtoms(other_idx, nb.GetIdx()) is not None
                and str(
                    mol.GetBondBetweenAtoms(other_idx, nb.GetIdx()).GetBondType()
                ).endswith("DOUBLE")
                for nb in other_atom.GetNeighbors()
                if nb.GetSymbol() == "C" and nb.GetIdx() != aryl_idx
            )
            if has_vinyl:
                return _build_retro_guide(
                    mol, bond, atom1_idx, atom2_idx, fragments,
                    bond_semantic="Ar-C (aryl-vinyl, Heck type)",
                    retro_pattern="Heck reaction",
                    retro_hint=(
                        "Heck retrosynthesis: after aryl-vinyl C-C cleavage, the aryl "
                        "end should become Ar-Br or Ar-I (halide precursor); the vinyl "
                        "end stays as an alkene. Pd-catalyzed, base (Et3N/K2CO3), "
                        "solvent (DMF/NMP)."
                    ),
                )

        # N-alkylation: N-C single bond (non-amide, non-sulfonamide)
        if (
            {atom1_sym, atom2_sym} == {"N", "C"}
            and bond_type == "SINGLE"
            and not _is_amide_bond(mol, bond)
            and not _is_sulfonamide_bond(mol, bond)
        ):
            n_idx = atom1_idx if atom1_sym == "N" else atom2_idx
            c_idx = atom2_idx if atom1_sym == "N" else atom1_idx
            c_atom = mol.GetAtomWithIdx(c_idx)
            is_benzylic = any(
                nb.GetIsAromatic()
                for nb in c_atom.GetNeighbors()
                if nb.GetIdx() != n_idx
            )
            halide_hint = (
                "benzylic halide (BnCl/BnBr)"
                if is_benzylic
                else "alkyl halide (R-Cl/R-Br/R-OTs)"
            )
            return _build_retro_guide(
                mol, bond, atom1_idx, atom2_idx, fragments,
                bond_semantic="N-C (non-amide/sulfonamide)",
                retro_pattern="N-alkylation",
                retro_hint=(
                    f"N-alkylation retrosynthesis: after N-C cleavage, the N end stays "
                    f"as N-H (nucleophile); the C end should become a {halide_hint} "
                    f"(electrophile). "
                    f"{'Benzylic C: suggest BnCl/BnBr. ' if is_benzylic else ''}"
                    f"Basic conditions (K2CO3/NaH), solvent (DMF/THF)."
                ),
                extra_info={
                    "broken_N_in_ring": mol.GetAtomWithIdx(n_idx).IsInRing(),
                    "broken_C_is_benzylic": is_benzylic,
                },
            )

        return None
