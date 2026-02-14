"""
Structure Analyzer
==================
Consolidates internal/rdkit_utils/analyze_structure.py and calculate_features.py
into a single module with typed dataclass outputs.

Public API:
    analyze_molecule(smiles) → AnalysisResult
    build_atom_bond_map(smiles) → AtomBondMap
"""

from __future__ import annotations

import logging
from typing import Any, Dict, List, Optional, Tuple

from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors

from tools.chem.mol_parser import parse_smiles, normalize_smiles
from tools.models import (
    AnalysisResult,
    AtomBondMap,
    AtomInfo,
    BondInfo,
    MoleculeInfo,
    RingInfo,
)

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def analyze_molecule(smiles: str) -> AnalysisResult:
    """Full molecule analysis returning a typed AnalysisResult.

    Steps:
        1. Parse & validate SMILES
        2. Compute molecular weight and SA score
        3. Build AtomBondMap (atoms, bonds with labels, rings)
        4. Identify functional groups (lightweight)
        5. Return AnalysisResult wrapping MoleculeInfo
    """
    if not smiles or not smiles.strip():
        return AnalysisResult(success=False, error="Empty SMILES")

    mol = parse_smiles(smiles)
    if mol is None:
        return AnalysisResult(success=False, error="Invalid SMILES: RDKit cannot parse")

    try:
        canonical = normalize_smiles(smiles)
        # Always re-parse canonical SMILES so atom indices are deterministic
        # and consistent with break_bond / propose_disconnection.
        mol_c = parse_smiles(canonical)
        if mol_c is None:
            mol_c = mol  # fallback: should never happen
        mol_wt = round(Descriptors.MolWt(mol_c), 3)
        sa_score = _calculate_sa_score(mol_c)
        atom_bond_map = _build_atom_bond_map_from_mol(mol_c, canonical)
        functional_groups = _identify_functional_groups(mol_c)

        mol_info = MoleculeInfo(
            smiles=smiles,
            canonical_smiles=canonical,
            molecular_weight=mol_wt,
            sa_score=sa_score,
            atom_bond_map=atom_bond_map,
            functional_groups=functional_groups,
        )
        return AnalysisResult(success=True, molecule_info=mol_info)
    except Exception as exc:
        logger.exception("analyze_molecule failed for %s", smiles)
        return AnalysisResult(success=False, error=str(exc))



def build_atom_bond_map(smiles: str) -> AtomBondMap:
    """Build an AtomBondMap for *smiles*.

    Returns an empty AtomBondMap (with error-like retro_guidance) on failure.
    """
    mol = parse_smiles(smiles)
    if mol is None:
        return AtomBondMap(retro_guidance="无法解析 SMILES")

    canonical = normalize_smiles(smiles)
    # Re-parse canonical to ensure atom indices are deterministic
    mol_c = parse_smiles(canonical)
    if mol_c is None:
        mol_c = mol  # fallback
    return _build_atom_bond_map_from_mol(mol_c, canonical)


# ---------------------------------------------------------------------------
# SA Score (from calculate_features.py)
# ---------------------------------------------------------------------------


def _calculate_sa_score(mol) -> float:
    """SA Score (1-10, lower = easier). Tries RDKit contrib, falls back to heuristic."""
    try:
        from rdkit.Chem import RDConfig
        import os
        import importlib.util

        sascore_path = os.path.join(RDConfig.RDContribDir, "SA_Score", "sascorer.py")
        if os.path.exists(sascore_path):
            spec = importlib.util.spec_from_file_location("sascorer", sascore_path)
            sascorer = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(sascorer)
            return round(sascorer.calculateScore(mol), 3)
    except Exception:
        pass

    # Improved fallback heuristic (considers rings, heteroatoms, chirality, MW)
    num_rings = rdMolDescriptors.CalcNumRings(mol)
    num_stereo = len(Chem.FindMolChiralCenters(mol, includeUnassigned=True))
    mol_wt = Descriptors.MolWt(mol)
    num_hetero = rdMolDescriptors.CalcNumHeteroatoms(mol)
    sa = 2.0 + 0.3 * num_rings + 0.5 * num_stereo + 0.01 * max(0, mol_wt - 150) + 0.1 * num_hetero
    return round(min(10.0, max(1.0, sa)), 3)


# ---------------------------------------------------------------------------
# Functional group identification (lightweight)
# ---------------------------------------------------------------------------

_FG_SMARTS: Dict[str, str] = {
    "hydroxyl": "[OX2H]",
    "carboxyl": "[CX3](=O)[OX2H1]",
    "amine_primary": "[NX3H2]",
    "amine_secondary": "[NX3H1]([#6])[#6]",
    "amide": "[NX3][CX3](=[OX1])",
    "ester": "[#6][CX3](=O)[OX2][#6]",
    "ketone": "[#6][CX3](=O)[#6]",
    "aldehyde": "[CX3H1](=O)",
    "nitro": "[$([NX3](=O)=O),$([NX3+](=O)[O-])]",
    "sulfonamide": "[SX4](=[OX1])(=[OX1])[NX3]",
    "nitrile": "[CX2]#[NX1]",
    "halide": "[F,Cl,Br,I]",
    "tetrazole": "c1nnn[nH]1",
}


def _identify_functional_groups(mol) -> Dict[str, int]:
    """Return {group_name: count} for common functional groups."""
    groups: Dict[str, int] = {}
    for name, smarts in _FG_SMARTS.items():
        pat = Chem.MolFromSmarts(smarts)
        if pat is None:
            continue
        matches = mol.GetSubstructMatches(pat)
        if matches:
            groups[name] = len(matches)
    return groups


# ---------------------------------------------------------------------------
# Bond classification helpers (from builtin.py)
# ---------------------------------------------------------------------------


def _is_ester_bond(mol, bond) -> bool:
    """C(=O)-O-C ester bond."""
    a1, a2 = bond.GetBeginAtom(), bond.GetEndAtom()

    def _check(c, o):
        if c.GetSymbol() != "C" or o.GetSymbol() != "O":
            return False
        for nb in c.GetNeighbors():
            if nb.GetIdx() == o.GetIdx():
                continue
            if nb.GetSymbol() == "O":
                b = mol.GetBondBetweenAtoms(c.GetIdx(), nb.GetIdx())
                if b and str(b.GetBondType()).endswith("DOUBLE"):
                    # O must connect to another C
                    for nb2 in o.GetNeighbors():
                        if nb2.GetIdx() != c.GetIdx() and nb2.GetSymbol() == "C":
                            return True
        return False

    return _check(a1, a2) or _check(a2, a1)


def _is_amide_bond(mol, bond) -> bool:
    """C(=O)-N amide bond."""
    a1, a2 = bond.GetBeginAtom(), bond.GetEndAtom()

    def _check(c, n):
        if c.GetSymbol() != "C" or n.GetSymbol() != "N":
            return False
        for nb in c.GetNeighbors():
            if nb.GetIdx() == n.GetIdx():
                continue
            if nb.GetSymbol() == "O":
                b = mol.GetBondBetweenAtoms(c.GetIdx(), nb.GetIdx())
                if b and str(b.GetBondType()).endswith("DOUBLE"):
                    return True
        return False

    return _check(a1, a2) or _check(a2, a1)


def _is_alpha_carbonyl_cc(mol, bond) -> bool:
    """C-C bond where one C is alpha to a carbonyl."""
    a1, a2 = bond.GetBeginAtom(), bond.GetEndAtom()
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
    """-SO2-N- sulfonamide bond."""
    a1, a2 = bond.GetBeginAtom(), bond.GetEndAtom()

    def _check(s, n):
        if s.GetSymbol() != "S" or n.GetSymbol() != "N":
            return False
        double_o = 0
        for nb in s.GetNeighbors():
            if nb.GetIdx() == n.GetIdx():
                continue
            if nb.GetSymbol() == "O":
                b = mol.GetBondBetweenAtoms(s.GetIdx(), nb.GetIdx())
                if b and str(b.GetBondType()).endswith("DOUBLE"):
                    double_o += 1
        return double_o >= 2

    return _check(a1, a2) or _check(a2, a1)


def _is_sulfonyl_ester_bond(mol, bond) -> bool:
    """-SO2-O-R sulfonyl ester bond."""
    if str(bond.GetBondType()).endswith("DOUBLE"):
        return False
    a1, a2 = bond.GetBeginAtom(), bond.GetEndAtom()

    def _check(s, o):
        if s.GetSymbol() != "S" or o.GetSymbol() != "O":
            return False
        double_o = 0
        for nb in s.GetNeighbors():
            if nb.GetIdx() == o.GetIdx():
                continue
            if nb.GetSymbol() == "O":
                b = mol.GetBondBetweenAtoms(s.GetIdx(), nb.GetIdx())
                if b and str(b.GetBondType()).endswith("DOUBLE"):
                    double_o += 1
        if double_o < 2:
            return False
        for nb in o.GetNeighbors():
            if nb.GetIdx() != s.GetIdx() and nb.GetSymbol() != "H":
                return True
        return o.GetDegree() > 1

    return _check(a1, a2) or _check(a2, a1)


# ---------------------------------------------------------------------------
# Bond annotation (from builtin.py AnalyzeMoleculeSkill._annotate_bond)
# ---------------------------------------------------------------------------


def _annotate_bond(
    mol, bond, atom1, atom2,
    sym1: str, sym2: str, bond_type: str, in_ring: bool,
) -> Tuple[str, List[str], str, int]:
    """Return (chemical_label, potential_reactions, retro_hint, strategic_priority)."""
    reactions: List[str] = []
    label = f"{sym1}-{sym2} ({bond_type})"
    hint = ""
    priority = 3

    a1_ar = atom1.GetIsAromatic()
    a2_ar = atom2.GetIsAromatic()

    if _is_ester_bond(mol, bond):
        label, reactions = "酯键 C(=O)-O", ["酯水解/Fischer酯化", "转酯化"]
        hint, priority = "断裂 -> 羧酸 + 醇 (经典逆合成断裂)", 1

    elif _is_amide_bond(mol, bond):
        label, reactions = "酰胺键 C(=O)-N", ["酰胺偶联 (EDC/HOBt)", "Schotten-Baumann"]
        hint, priority = "断裂 -> 羧酸 + 胺", 1

    elif _is_sulfonamide_bond(mol, bond):
        label, reactions = "磺酰胺键 SO2-N", ["磺酰胺形成 (磺酰氯+胺)", "Buchwald磺酰胺化"]
        hint = "断裂 -> 磺酰氯 (ArSO2Cl) + 胺 (RNH2)"
        priority = 1

    elif _is_sulfonyl_ester_bond(mol, bond):
        label, reactions = "磺酸酯键 SO2-O", ["磺酰化 (磺酰氯+醇)", "甲磺酰化/甲苯磺酰化"]
        hint, priority = "断裂 -> 磺酰氯 + 醇", 1

    elif {sym1, sym2} == {"C", "N"} and not _is_amide_bond(mol, bond):
        n_atom = atom1 if sym1 == "N" else atom2
        ar_atom = atom2 if sym1 == "N" else atom1
        if a1_ar or a2_ar:
            n_has_other_ar = any(
                nb.GetIsAromatic() and nb.GetIdx() != ar_atom.GetIdx()
                for nb in n_atom.GetNeighbors()
            )
            if n_has_other_ar:
                label, reactions = "二芳胺键 Ar-NH-Ar", ["Buchwald-Hartwig偶联", "Ullmann偶联", "SNAr取代"]
                hint, priority = "断裂 -> 芳基卤化物 + 芳基胺 (Pd催化或SNAr)", 1
            else:
                label, reactions = "芳基-N 键", ["Buchwald-Hartwig胺化", "Ullmann偶联", "Chan-Lam偶联"]
                hint, priority = "断裂 -> 芳基卤化物 + 胺 (Pd/Cu催化)", 2
        else:
            label, reactions = "烷基-N 键", ["还原胺化", "SN2取代", "Mannich反应"]
            hint, priority = "断裂 -> 醛/酮 + 胺 或 卤代物 + 胺", 2

    elif {sym1, sym2} == {"C", "O"} and bond_type == "SINGLE" and not _is_ester_bond(mol, bond):
        if a1_ar or a2_ar:
            label, reactions = "芳基-O 醚键", ["Williamson醚合成", "SNAr取代", "Mitsunobu反应"]
            hint, priority = "断裂 -> 酚/醇 + 烷基卤化物", 2
        else:
            label, reactions = "烷基-O 醚键", ["Williamson醚合成", "SN2取代"]
            hint, priority = "断裂 -> 醇 + 卤代物", 2

    elif {sym1, sym2} == {"N"} or (sym1 == "N" and sym2 == "N"):
        if bond_type == "SINGLE":
            label, reactions = "N-N 单键 (肼键)", ["肼形成", "偶氮偶联还原", "Curtius重排"]
            hint, priority = "断裂 -> 肼 + 羰基 或 两个胺片段", 2
        elif bond_type == "DOUBLE":
            label, reactions = "N=N 双键 (偶氮)", ["偶氮偶联", "重氮化"]
            hint, priority = "断裂 -> 芳基重氮盐 + 酚/胺", 2

    elif {sym1, sym2} == {"C", "S"}:
        label, reactions = "C-S 键", ["硫醚形成", "SN2取代"]
        hint, priority = "断裂 -> 硫醇 + 卤代物", 2

    elif {sym1, sym2} == {"S", "N"} and not _is_sulfonamide_bond(mol, bond):
        label, reactions = "S-N 键", ["S-N 偶联"]
        hint, priority = "断裂 -> 含硫亲电体 + 胺", 2

    elif {sym1, sym2} == {"N", "O"} and bond_type == "SINGLE":
        label, reactions = "N-O 键", ["N-氧化物还原", "Cope消除"]
        hint, priority = "断裂 -> 胺 + 含氧片段", 3

    elif {sym1, sym2} == {"P", "O"} and bond_type == "SINGLE":
        label, reactions = "P-O 键 (磷酸酯)", ["磷酸酯化", "Mitsunobu反应"]
        hint, priority = "断裂 -> 磷酰氯 + 醇", 2

    elif {sym1, sym2} == {"C", "B"}:
        label, reactions = "C-B 键 (硼化)", ["Miyaura硼化", "Suzuki偶联前体"]
        hint, priority = "断裂 -> 芳基/烷基 + 硼酸/硼酸酯", 2

    elif {sym1, sym2} == {"B", "O"}:
        label, reactions = "B-O 键 (硼酸酯)", ["硼酸酯化"]
        hint, priority = "断裂 -> 硼酸 + 醇", 2

    elif sym1 == "C" and sym2 == "C" and bond_type == "DOUBLE":
        label, reactions = "C=C 双键", ["Wittig反应", "HWE反应", "烯烃复分解", "Heck反应"]
        hint, priority = "断裂 -> 醛/酮 + 磷叶立德 或 两个烯烃片段", 2

    elif sym1 == "C" and sym2 == "C" and bond_type == "SINGLE":
        if a1_ar and a2_ar:
            label, reactions = "联芳 C-C 键", ["Suzuki偶联", "Negishi偶联", "Kumada偶联"]
            hint, priority = "断裂 -> 芳基卤化物 + 芳基金属试剂 (Pd催化)", 1
        elif a1_ar or a2_ar:
            label, reactions = "芳基-烷基 C-C 键", ["Friedel-Crafts", "Heck反应", "Negishi偶联"]
            hint, priority = "断裂 -> 芳烃 + 烷基卤化物/酰氯", 2
        elif _is_alpha_carbonyl_cc(mol, bond):
            label, reactions = "α-羰基 C-C 键", ["Aldol反应", "Michael加成", "Claisen缩合"]
            hint, priority = "断裂 -> 两个羰基片段 (烯醇 + 亲电体)", 1
        else:
            label, reactions = "烷基 C-C 键", ["Grignard加成", "烷基化"]
            hint, priority = "断裂 -> 格氏试剂 + 羰基 或 两个自由基片段", 3

    elif {sym1, sym2} == {"Si", "O"}:
        label, reactions = "Si-O 键 (硅醚保护基)", ["TBS/TIPS/TMS保护", "TBAF脱保护"]
        hint, priority = "断裂 -> 自由醇 + 硅基试剂", 2

    # Ring bond fallback
    if in_ring and not hint:
        label = f"环内 {sym1}-{sym2} 键"
        hint = "环键断裂只能开环。考虑环形成的逆反应。"
        priority = 3

    return label, reactions, hint, priority


# ---------------------------------------------------------------------------
# Core: build AtomBondMap from an RDKit Mol
# ---------------------------------------------------------------------------


def _classify_ring(mol, ring_atoms: tuple, ring_id: int) -> RingInfo:
    """Classify a single ring and produce formation hints."""
    size = len(ring_atoms)
    is_aromatic = all(mol.GetAtomWithIdx(a).GetIsAromatic() for a in ring_atoms)
    heteroatoms = [
        mol.GetAtomWithIdx(a).GetSymbol()
        for a in ring_atoms
        if mol.GetAtomWithIdx(a).GetSymbol() != "C"
    ]

    # Ring type label
    if is_aromatic and size == 6 and not heteroatoms:
        ring_type = "苯环"
    elif is_aromatic and heteroatoms:
        ring_type = f"芳香杂环 ({size}元, 含 {', '.join(heteroatoms)})"
    elif heteroatoms:
        ring_type = f"杂环 ({size}元, 含 {', '.join(heteroatoms)})"
    else:
        ring_type = f"碳环 ({size}元)"

    # Formation hints
    hints: List[str] = []
    het_set = set(heteroatoms)
    if size == 3 and "O" in het_set:
        hints.append("环氧化物 -> 烯烃 + 氧化剂 (mCPBA/Sharpless)")
    elif size == 5 and heteroatoms:
        if "N" in het_set:
            hints.append("含氮五元环 -> Paal-Knorr / Fischer吲哚 / 1,3-偶极环加成")
        if "O" in het_set:
            hints.append("含氧五元环 -> 分子内环化/缩合/内酯化")
        if "S" in het_set:
            hints.append("含硫五元环 -> 噻吩合成 / Gewald反应")
        hints.append(f"{size}元杂环 -> 考虑分子内环化或[3+2]环加成")
    elif size == 6 and is_aromatic and heteroatoms:
        hints.append("芳香杂环 -> Hantzsch / Chichibabin / 环化缩合")
    elif size == 6 and not is_aromatic:
        hints.append("六元环 -> Diels-Alder [4+2] / RCM / 分子内环化")
    elif size >= 7:
        hints.append(f"大环 ({size}元) -> RCM / 大环内酯化 / 模板效应")

    return RingInfo(
        ring_id=ring_id,
        size=size,
        atom_indices=list(ring_atoms),
        is_aromatic=is_aromatic,
        ring_type=ring_type,
        formation_hints=hints,
    )


def _build_atom_bond_map_from_mol(mol, canonical: str) -> AtomBondMap:
    """Build a typed AtomBondMap from an already-parsed RDKit Mol."""
    ring_info_obj = mol.GetRingInfo()
    atom_rings = ring_info_obj.AtomRings()

    # Ring classification
    rings: List[RingInfo] = [
        _classify_ring(mol, ring_atoms, rid)
        for rid, ring_atoms in enumerate(atom_rings)
    ]

    # Atoms
    atoms: List[AtomInfo] = []
    for atom in mol.GetAtoms():
        atoms.append(AtomInfo(
            idx=atom.GetIdx(),
            symbol=atom.GetSymbol(),
            aromatic=atom.GetIsAromatic(),
            in_ring=atom.IsInRing(),
            degree=atom.GetDegree(),
            num_hs=atom.GetTotalNumHs(),
            neighbors=[n.GetIdx() for n in atom.GetNeighbors()],
        ))

    # Bonds with semantic annotation
    bonds: List[BondInfo] = []
    strategic_bonds: List[Dict[str, Any]] = []
    breakable_count = 0

    for bond in mol.GetBonds():
        a1 = bond.GetBeginAtomIdx()
        a2 = bond.GetEndAtomIdx()
        in_ring = bond.IsInRing()
        is_aromatic = bond.GetIsAromatic()
        bt_str = str(bond.GetBondType()).replace("rdkit.Chem.rdchem.BondType.", "")
        breakable = not is_aromatic

        # Exclude bonds that should not be broken
        if breakable:
            s1 = mol.GetAtomWithIdx(a1).GetSymbol()
            s2 = mol.GetAtomWithIdx(a2).GetSymbol()
            syms = {s1, s2}
            if (syms == {"S", "O"} and bt_str == "DOUBLE") or \
               (syms == {"P", "O"} and bt_str == "DOUBLE") or \
               (syms == {"C", "O"} and bt_str == "DOUBLE") or \
               (syms == {"N", "O"} and bt_str == "DOUBLE"):
                breakable = False

        chem_label = ""
        potential_rxns: List[str] = []
        retro_hint = ""
        strat_priority = 3

        if breakable:
            breakable_count += 1
            atom1 = mol.GetAtomWithIdx(a1)
            atom2 = mol.GetAtomWithIdx(a2)
            chem_label, potential_rxns, retro_hint, strat_priority = _annotate_bond(
                mol, bond, atom1, atom2,
                atom1.GetSymbol(), atom2.GetSymbol(), bt_str, in_ring,
            )
            if strat_priority <= 2:
                strategic_bonds.append({
                    "atom1_idx": a1, "atom2_idx": a2,
                    "chemical_label": chem_label,
                    "potential_reactions": potential_rxns,
                    "retro_hint": retro_hint,
                    "strategic_priority": strat_priority,
                    "in_ring": in_ring,
                })

        bonds.append(BondInfo(
            idx=bond.GetIdx(),
            atom1_idx=a1,
            atom2_idx=a2,
            bond_type=bt_str,
            in_ring=in_ring,
            breakable=breakable,
            chemical_label=chem_label,
            potential_reactions=potential_rxns,
            retro_hint=retro_hint,
            strategic_priority=strat_priority,
        ))

    retro_guidance = _build_retro_guidance(mol, rings, strategic_bonds, breakable_count)

    return AtomBondMap(
        canonical_smiles=canonical,
        atoms=atoms,
        bonds=bonds,
        ring_info=rings,
        retro_guidance=retro_guidance,
    )


# ---------------------------------------------------------------------------
# Retro guidance summary (from builtin.py _build_retro_guidance)
# ---------------------------------------------------------------------------


def _build_retro_guidance(
    mol, rings: List[RingInfo],
    strategic_bonds: List[Dict[str, Any]], breakable_count: int,
) -> str:
    """Generate a retro-synthesis guidance summary string."""
    lines: List[str] = []

    if strategic_bonds:
        sorted_bonds = sorted(strategic_bonds, key=lambda x: x["strategic_priority"])
        lines.append("【推荐断键位置 (按优先级排序)】")
        for i, b in enumerate(sorted_bonds, 1):
            a1, a2 = b["atom1_idx"], b["atom2_idx"]
            rxns = ", ".join(b.get("potential_reactions", [])[:3])
            lines.append(f"  {i}. atom {a1}-{a2} ({b['chemical_label']}), priority={b['strategic_priority']}")
            lines.append(f"     反应: {rxns}")
            if b.get("retro_hint"):
                lines.append(f"     提示: {b['retro_hint']}")

        lines.append("")
        lines.append("【建议操作顺序】")
        non_ring = [b for b in sorted_bonds if not b.get("in_ring")]
        ring_bonds = [b for b in sorted_bonds if b.get("in_ring")]
        step = 1
        for b in non_ring:
            a1, a2 = b["atom1_idx"], b["atom2_idx"]
            lines.append(f"  Step {step}: break_bond atom1_idx={a1} atom2_idx={a2} ({b['chemical_label']})")
            step += 1
        for b in ring_bonds:
            a1, a2 = b["atom1_idx"], b["atom2_idx"]
            lines.append(f"  Step {step}: break_bond atom1_idx={a1} atom2_idx={a2} ({b['chemical_label']}) [环键]")
            step += 1
    else:
        lines.append("【注意】未发现高优先级断键位置。")
        lines.append("  建议: 使用 propose_disconnection 获取规则匹配提案。")

    # Heterocyclic ring hints
    hetero_rings = [r for r in rings if not r.is_aromatic and r.formation_hints]
    if hetero_rings:
        lines.append("")
        lines.append("【杂环系统 -- 可能需要环形成逆反应】")
        for r in hetero_rings:
            lines.append(f"  - Ring {r.ring_id}: {r.ring_type}")
            for fh in r.formation_hints:
                lines.append(f"    -> {fh}")

    # Fused ring hints
    atom_ring_count: Dict[int, int] = {}
    for r in rings:
        for aidx in r.atom_indices:
            atom_ring_count[aidx] = atom_ring_count.get(aidx, 0) + 1
    fused_atoms = [a for a, c in atom_ring_count.items() if c >= 2]
    if fused_atoms:
        lines.append("")
        lines.append("【稠合/桥接环系统】")
        lines.append(f"  共 {len(fused_atoms)} 个原子属于多个环 (稠合/桥接)。")
        lines.append("  建议: 考虑环形成的逆反应 (分子内环化、缩合、环加成的逆反应)。")

    lines.append("")
    lines.append("【工具使用提示】")
    lines.append("  - break_bond 做拓扑断键 (加H补价)，某些键类型需要手动修正前体")
    lines.append("  - 每次 break_bond 后务必 validate_reaction 验证原子平衡")

    return "\n".join(lines)
