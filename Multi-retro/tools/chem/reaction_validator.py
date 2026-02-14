"""
Reaction Validation Module
==========================
Comprehensive reaction validation with structured output.
Migrated from internal/analysis/reaction_validator.py.

18 byproduct patterns, 12 selectivity hints, three-tier verdict (PASS/CONDITIONAL/FAIL).
Returns typed ValidationResult from tools.models.skill_models.

Public API:
    validate_reaction(reaction_smiles, ...) -> ValidationResult
"""

from __future__ import annotations

import logging
from collections import Counter
from typing import Any, Dict, List, Optional

from rdkit import Chem

from tools.chem.mol_parser import normalize_reaction_smiles
from tools.models.skill_models import ValidationResult

logger = logging.getLogger(__name__)


def _element_counter_from_smiles(smiles: str) -> Counter:
    """Count elements (including implicit H) from a dot-separated SMILES."""
    counter: Counter = Counter()
    for frag in smiles.split("."):
        frag = frag.strip()
        if not frag:
            continue
        mol = Chem.MolFromSmiles(frag)
        if mol is None:
            continue
        mol = Chem.AddHs(mol)
        for atom in mol.GetAtoms():
            counter[atom.GetSymbol()] += 1
    return counter


# 18 byproduct patterns (original 11 + 7 new)
_BYPRODUCT_PATTERNS = [
    {"formula": "H2", "elements": {"H": 2}, "name": "Hydrogen", "smiles": "[H][H]"},
    {"formula": "H2O", "elements": {"H": 2, "O": 1}, "name": "Water", "smiles": "O"},
    {"formula": "HCl", "elements": {"H": 1, "Cl": 1}, "name": "Hydrochloric acid", "smiles": "Cl"},
    {"formula": "HBr", "elements": {"H": 1, "Br": 1}, "name": "Hydrobromic acid", "smiles": "Br"},
    {"formula": "HI", "elements": {"H": 1, "I": 1}, "name": "Hydroiodic acid", "smiles": "I"},
    {"formula": "HF", "elements": {"H": 1, "F": 1}, "name": "Hydrofluoric acid", "smiles": "F"},
    {"formula": "NH3", "elements": {"N": 1, "H": 3}, "name": "Ammonia", "smiles": "N"},
    {"formula": "CO2", "elements": {"C": 1, "O": 2}, "name": "Carbon dioxide", "smiles": "O=C=O"},
    {"formula": "CO", "elements": {"C": 1, "O": 1}, "name": "Carbon monoxide", "smiles": "[C-]#[O+]"},
    {"formula": "SO2", "elements": {"S": 1, "O": 2}, "name": "Sulfur dioxide", "smiles": "O=S=O"},
    {"formula": "N2", "elements": {"N": 2}, "name": "Nitrogen", "smiles": "N#N"},
    {"formula": "NaCl", "elements": {"Na": 1, "Cl": 1}, "name": "Sodium chloride", "smiles": "[Na]Cl"},
    {"formula": "NaBr", "elements": {"Na": 1, "Br": 1}, "name": "Sodium bromide", "smiles": "[Na]Br"},
    {"formula": "KBr", "elements": {"K": 1, "Br": 1}, "name": "Potassium bromide", "smiles": "[K]Br"},
    {"formula": "LiCl", "elements": {"Li": 1, "Cl": 1}, "name": "Lithium chloride", "smiles": "[Li]Cl"},
    {"formula": "CH3OH", "elements": {"C": 1, "H": 4, "O": 1}, "name": "Methanol", "smiles": "CO"},
    {"formula": "C2H5OH", "elements": {"C": 2, "H": 6, "O": 1}, "name": "Ethanol", "smiles": "CCO"},
    {"formula": "AcOH", "elements": {"C": 2, "H": 4, "O": 2}, "name": "Acetic acid", "smiles": "CC(=O)O"},
]

# 12 selectivity hint sets
_SELECTIVITY_HINTS: Dict[str, List[str]] = {
    "esterification": [
        "Carboxylic acid is more reactive than alcohol",
        "Consider acid catalyst for faster reaction",
    ],
    "amide_coupling": [
        "Primary amines more reactive than secondary",
        "Watch for racemization at alpha-carbon",
    ],
    "nucleophilic_substitution": [
        "SN2: steric hindrance affects rate",
        "SN1: carbocation stability determines pathway",
    ],
    "oxidation": [
        "Primary alcohols -> aldehydes -> carboxylic acids",
        "Control conditions to stop at aldehyde stage",
    ],
    "reduction": [
        "Nitriles reduce to amines",
        "Aldehydes reduce to alcohols",
    ],
    "cross_coupling": [
        "Palladium-catalyzed: halide leaving group order Br > I > Cl",
        "Aryl halides more reactive than alkyl",
    ],
    "diels_alder": [
        "Endo product kinetically favored, exo thermodynamically favored",
        "Electron-rich diene + electron-poor dienophile preferred",
    ],
    "click": [
        "CuAAC gives 1,4-disubstituted triazole selectively",
        "RuAAC gives 1,5-disubstituted triazole",
    ],
    "grignard": [
        "1,2-addition preferred for aldehydes/ketones",
        "1,4-addition possible with alpha,beta-unsaturated carbonyls",
    ],
    "protection": [
        "Boc: acid-labile (TFA/HCl)",
        "Cbz: hydrogenolysis (H2/Pd-C)",
        "Fmoc: base-labile (piperidine)",
    ],
    "epoxidation": [
        "mCPBA: cis-selective for cyclic alkenes",
        "Sharpless: enantioselective for allylic alcohols",
    ],
    "metathesis": [
        "Ring-closing metathesis: 5-7 membered rings preferred",
        "Cross metathesis: selectivity depends on olefin type",
    ],
}


def _infer_byproducts(delta: Dict[str, int]) -> List[Dict[str, Any]]:
    """Infer possible byproducts from atom imbalance delta.

    Each returned dict includes ``_remaining_imbalance`` for verdict logic.
    """
    byproducts: List[Dict[str, Any]] = []
    excess_reactant = {k: v for k, v in delta.items() if v > 0}
    excess_product = {k: -v for k, v in delta.items() if v < 0}

    remaining = dict(excess_reactant)
    for pattern in _BYPRODUCT_PATTERNS:
        while all(remaining.get(e, 0) >= c for e, c in pattern["elements"].items()):
            existing = next((bp for bp in byproducts if bp["formula"] == pattern["formula"]), None)
            if existing:
                existing["count"] = existing.get("count", 1) + 1
            else:
                byproducts.append({
                    "formula": pattern["formula"],
                    "name": pattern["name"],
                    "smiles": pattern["smiles"],
                    "count": 1,
                    "type": "byproduct",
                    "confidence": 0.8,
                })
            for elem, count in pattern["elements"].items():
                remaining[elem] = remaining.get(elem, 0) - count
                if remaining[elem] == 0:
                    del remaining[elem]

    if excess_product:
        for pattern in _BYPRODUCT_PATTERNS:
            while all(excess_product.get(e, 0) >= c for e, c in pattern["elements"].items()):
                byproducts.append({
                    "formula": pattern["formula"],
                    "name": pattern["name"],
                    "smiles": pattern["smiles"],
                    "count": 1,
                    "type": "missing_reactant",
                    "confidence": 0.6,
                })
                for elem, count in pattern["elements"].items():
                    excess_product[elem] = excess_product.get(elem, 0) - count
                    if excess_product[elem] == 0:
                        del excess_product[elem]

    unexplained = sum(abs(v) for v in remaining.values()) + sum(excess_product.values())
    if unexplained > 0 and byproducts:
        penalty = min(0.4, unexplained * 0.1)
        for bp in byproducts:
            bp["confidence"] = round(bp["confidence"] - penalty, 2)

    for bp in byproducts:
        bp["_remaining_imbalance"] = unexplained

    return byproducts


def _generate_selectivity_notes(
    reaction_type: str, bond_changes: Dict[str, List]
) -> List[str]:
    """Generate selectivity notes based on reaction type and bond changes."""
    notes: List[str] = []
    for key, hints in _SELECTIVITY_HINTS.items():
        if key.lower() in reaction_type.lower():
            notes.extend(hints)
            break

    if bond_changes.get("formed") and len(bond_changes["formed"]) > 1:
        notes.append("Multiple bonds forming - check regioselectivity")
    if bond_changes.get("broken") and len(bond_changes["broken"]) > 1:
        notes.append("Multiple bonds breaking - consider reaction order")

    return notes


def validate_reaction(
    reaction_smiles: str,
    reaction_class: str = "",
    include_bond_changes: bool = True,
    include_selectivity: bool = True,
) -> ValidationResult:
    """Validate a reaction and return a typed ValidationResult.

    Three-tier verdict:
        PASS        - atom balance OK, no issues
        CONDITIONAL - imbalance explainable by inferred byproducts
        FAIL        - unexplainable imbalance or format error
    """
    issues: List[Dict[str, Any]] = []
    byproducts: List[Dict[str, Any]] = []

    if ">>" not in reaction_smiles:
        issues.append({"code": "REACTION_FORMAT_INVALID", "message": "Missing '>>' separator"})
        return ValidationResult(
            success=False, verdict="FAIL", is_valid=False,
            issues=issues, byproducts=[], error="",
        )

    try:
        parts = reaction_smiles.split(">>")
        if len(parts) != 2:
            issues.append({"code": "REACTION_FORMAT_INVALID", "message": "Invalid reaction format"})
            return ValidationResult(
                success=False, verdict="FAIL", is_valid=False,
                issues=issues, byproducts=[], error="",
            )

        reactants, products = parts

        # 1. Atom conservation check
        reactant_counter = _element_counter_from_smiles(reactants)
        product_counter = _element_counter_from_smiles(products)
        all_elements = set(reactant_counter.keys()) | set(product_counter.keys())
        delta = {
            elem: reactant_counter.get(elem, 0) - product_counter.get(elem, 0)
            for elem in all_elements
        }
        delta = {k: v for k, v in delta.items() if v != 0}

        atom_balanced = True
        if delta:
            atom_balanced = False
            issues.append({
                "code": "ATOM_CONSERVATION_FAIL",
                "message": f"Atom imbalance: {dict(delta)}",
            })
            byproducts = _infer_byproducts(delta)

        # 2. Bond changes (optional, requires atom_mappers.utils)
        bond_changes: Dict[str, List] = {"formed": [], "broken": []}
        if include_bond_changes:
            try:
                from tools.legacy.atom_mappers import utils
                bond_changes = utils.get_bond_changes(reaction_smiles)
            except Exception as e:
                issues.append({"code": "BOND_ANALYSIS_FAILED", "message": str(e)})

        # 3. Reaction type classification (optional)
        predicted_type = reaction_class or "unknown"
        try:
            from tools.chem.reaction_classifier import classify_reaction
            type_info = classify_reaction(reactants, products, reaction_class)
            predicted_type = type_info.get("name", predicted_type)
        except Exception:
            pass

        # 4. Selectivity notes
        selectivity_notes: List[str] = []
        if include_selectivity:
            selectivity_notes = _generate_selectivity_notes(predicted_type, bond_changes)

        # 5. Determine verdict
        if not atom_balanced:
            if byproducts:
                remaining = byproducts[0].get("_remaining_imbalance", -1)
                verdict = "CONDITIONAL"
            else:
                verdict = "FAIL"
        elif issues:
            verdict = "CONDITIONAL"
        else:
            verdict = "PASS"

        # Strip internal keys from byproducts before returning
        clean_byproducts = [
            {k: v for k, v in bp.items() if not k.startswith("_")}
            for bp in byproducts
        ]

        return ValidationResult(
            success=True,
            verdict=verdict,
            is_valid=(verdict == "PASS"),
            issues=issues,
            byproducts=clean_byproducts,
            error="",
        )

    except Exception as e:
        issues.append({"code": "VALIDATION_ERROR", "message": str(e)})
        return ValidationResult(
            success=False, verdict="FAIL", is_valid=False,
            issues=issues, byproducts=[], error=str(e),
        )
