"""
Feasibility Assessor Module
============================
Forward synthesis feasibility assessment: FG compatibility, side reaction risks,
protecting group requirements, mechanistic plausibility.

Migrated from internal/analysis/feasibility_assessor.py.
Uses tools.chem.functional_groups for FG analysis instead of internal module.

Public API:
    assess_forward_feasibility(reaction_smiles, precursors=None, reaction_class="") -> Dict
"""

from __future__ import annotations

import logging
from typing import Any, Dict, List, Optional

logger = logging.getLogger(__name__)

try:
    from rdkit import Chem
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False


# ---- Incompatibility Rules Matrix ----
_INCOMPATIBILITY_RULES = [
    ("alcohol", "acyl_halide", "high", "Alcohol + acyl halide: competing esterification"),
    ("amine", "acyl_halide", "high", "Amine + acyl halide: competing amidation"),
    ("amine", "aldehyde", "medium", "Amine + aldehyde: possible imine formation"),
    ("alcohol", "carboxylic_acid", "medium", "Multiple OH may cause esterification selectivity issues"),
    ("thiol", "halide", "medium", "Thiol + halide: nucleophilic substitution risk"),
    ("amine", "epoxide", "high", "Amine opens epoxide ring"),
    ("alcohol", "aldehyde", "medium", "Alcohol + aldehyde: hemiacetal formation"),
]

# ---- Side Reaction Patterns ----
_SIDE_REACTION_PATTERNS = [
    ({"aldehyde": 1, "amine": 1}, "imine_formation", "medium",
     "Aldehyde + amine condensation to imine"),
    ({"epoxide": 1, "thiol": 1}, "thiol_ring_opening", "high",
     "Thiol rapidly opens epoxide"),
    ({"nitro": 1, "thiol": 1}, "redox_side_reaction", "high",
     "Nitro + thiol redox incompatibility"),
    ({"azide": 1}, "azide_decomposition", "high",
     "Azide may decompose under acidic/thermal conditions"),
    ({"diazo": 1}, "diazo_decomposition", "high",
     "Diazo compound unstable, may decompose"),
    ({"acyl_halide": 1}, "hydrolysis", "medium",
     "Acyl halide moisture-sensitive, may hydrolyze"),
    ({"anhydride": 1, "amine": 1}, "competing_acylation", "medium",
     "Anhydride + amine: competing acylation"),
    ({"aldehyde": 1, "alcohol": 1}, "hemiacetal_formation", "low",
     "Aldehyde + alcohol: reversible hemiacetal"),
    ({"ketone": 1, "amine": 1}, "enamine_formation", "medium",
     "Ketone + secondary amine: enamine byproduct"),
    ({"epoxide": 1, "amine": 1}, "amine_ring_opening", "high",
     "Amine nucleophilic ring opening of epoxide"),
    ({"epoxide": 1, "carboxylic_acid": 1}, "acid_ring_opening", "medium",
     "Acid-catalyzed epoxide ring opening"),
    ({"vinyl": 1, "amine": 1}, "michael_addition", "medium",
     "Activated alkene + amine: Michael addition"),
    ({"isocyanate": 1, "alcohol": 1}, "carbamate_formation", "high",
     "Isocyanate + alcohol: rapid carbamate formation"),
    ({"isocyanate": 1, "amine": 1}, "urea_formation", "high",
     "Isocyanate + amine: rapid urea formation"),
    ({"aldehyde": 2}, "aldol_self_condensation", "medium",
     "Multiple aldehydes: self-condensation (Aldol)"),
    ({"boronic_acid": 1, "alcohol": 1}, "boronate_ester", "low",
     "Boronic acid + diol: boronate ester may affect Suzuki coupling"),
    ({"halide": 1, "amine": 1}, "over_alkylation", "medium",
     "Halide + amine: over-alkylation risk"),
    ({"carboxylic_acid": 1, "amine": 1}, "salt_formation", "low",
     "Acid + amine: salt formation reduces nucleophilicity"),
]

# ---- Protecting Group Suggestions ----
_PG_SUGGESTIONS = [
    ("alcohol", "acyl_halide", "Protect OH to avoid esterification", "TBS or THP"),
    ("amine", "acyl_halide", "Protect amine to avoid amidation", "Boc or Fmoc"),
    ("amine", "aldehyde", "Protect amine to avoid imine formation", "Boc or Cbz"),
    ("amine", "epoxide", "Protect amine to avoid ring opening", "Boc or Fmoc"),
    ("alcohol", "aldehyde", "Protect OH to avoid hemiacetal", "TBS or Acetyl"),
    ("thiol", "halide", "Protect thiol to avoid substitution", "Trityl or Acm"),
    ("phenol", "acyl_halide", "Protect phenol to avoid O-acylation", "Bn or MOM"),
    ("alcohol", "isocyanate", "Protect OH to avoid carbamate", "TBS"),
    ("amine", "anhydride", "Protect amine to avoid competing acylation", "Boc or Fmoc"),
]


def _collect_functional_groups(precursor_list: List[str]) -> Dict[str, Dict[str, int]]:
    """Analyze FGs for each precursor using tools.chem.functional_groups."""
    from tools.chem.functional_groups import analyze_functional_groups

    result: Dict[str, Dict[str, int]] = {}
    for smi in precursor_list:
        fg_result = analyze_functional_groups(smi)
        result[smi] = fg_result.get("functional_groups", {}) if fg_result.get("success") else {}
    return result


def _aggregate_fgs(precursor_groups: Dict[str, Dict[str, int]]) -> Dict[str, int]:
    merged: Dict[str, int] = {}
    for fg_dict in precursor_groups.values():
        for fg_name, count in fg_dict.items():
            if count > 0:
                merged[fg_name] = merged.get(fg_name, 0) + count
    return merged


def _check_compatibility(all_fg: Dict[str, int]) -> List[Dict[str, Any]]:
    issues: List[Dict[str, Any]] = []
    for ga, gb, severity, desc in _INCOMPATIBILITY_RULES:
        if all_fg.get(ga, 0) > 0 and all_fg.get(gb, 0) > 0:
            issues.append({
                "type": f"{ga}_{gb}_incompatibility",
                "severity": severity,
                "description": desc,
                "affected_groups": [ga, gb],
            })
    return issues


def _assess_side_reactions(all_fg: Dict[str, int]) -> List[Dict[str, Any]]:
    risks: List[Dict[str, Any]] = []
    for trigger, risk_name, severity, desc in _SIDE_REACTION_PATTERNS:
        if all(all_fg.get(fg, 0) >= cnt for fg, cnt in trigger.items()):
            risks.append({"risk": risk_name, "severity": severity, "description": desc})
    return risks


def _detect_pg_requirements(all_fg: Dict[str, int]) -> List[Dict[str, Any]]:
    reqs: List[Dict[str, Any]] = []
    for sensitive, interfering, reason, suggested in _PG_SUGGESTIONS:
        if all_fg.get(sensitive, 0) > 0 and all_fg.get(interfering, 0) > 0:
            reqs.append({"group": sensitive, "reason": reason, "suggested_pg": suggested})
    return reqs


def _compute_score(
    compat_issues: List[Dict[str, Any]],
    side_risks: List[Dict[str, Any]],
) -> tuple:
    high = sum(1 for x in compat_issues + side_risks if x.get("severity") == "high")
    medium = sum(1 for x in compat_issues + side_risks if x.get("severity") == "medium")

    if high == 0:
        score = "high"
    elif high <= 1:
        score = "medium"
    else:
        score = "low"

    confidence = max(0.0, min(0.85, 1.0 - high * 0.3 - medium * 0.1))
    return score, round(confidence, 2)


def assess_forward_feasibility(
    reaction_smiles: str,
    precursors: Optional[List[str]] = None,
    reaction_class: str = "",
) -> Dict[str, Any]:
    """Assess forward synthesis feasibility of a retrosynthetic proposal.

    Args:
        reaction_smiles: Reaction SMILES (reactants>>products)
        precursors: Optional explicit precursor list (else parsed from reaction_smiles)
        reaction_class: Optional reaction class hint

    Returns:
        {
            "feasibility_score": "high" | "medium" | "low",
            "overall_confidence": float,
            "functional_group_analysis": {...},
            "side_reaction_risks": [...],
            "protecting_group_requirements": [...],
            "concerns": [...],
            "recommendations": [...],
            "success": bool,
            "error": str | None,
        }
    """
    empty = {
        "feasibility_score": "low", "overall_confidence": 0.0,
        "functional_group_analysis": {"precursor_groups": {}, "compatibility_issues": []},
        "side_reaction_risks": [], "protecting_group_requirements": [],
        "concerns": [], "recommendations": [], "success": False, "error": None,
    }

    if not reaction_smiles or ">>" not in reaction_smiles:
        empty["error"] = "Invalid reaction SMILES: missing '>>' separator"
        return empty

    try:
        # Parse precursors
        if precursors:
            precursor_list = [s.strip() for s in precursors if s and s.strip()]
        else:
            reactant_side = reaction_smiles.split(">>")[0].strip()
            precursor_list = [s.strip() for s in reactant_side.split(".") if s.strip()]

        if not precursor_list:
            empty["error"] = "No precursors found"
            return empty

        precursor_groups = _collect_functional_groups(precursor_list)
        all_fg = _aggregate_fgs(precursor_groups)
        compat_issues = _check_compatibility(all_fg)
        side_risks = _assess_side_reactions(all_fg)
        pg_reqs = _detect_pg_requirements(all_fg)
        score, confidence = _compute_score(compat_issues, side_risks)

        concerns = [x["description"] for x in compat_issues] + [x["description"] for x in side_risks]
        recommendations = [
            f"Consider protecting {pg['group']} with {pg['suggested_pg']}: {pg['reason']}"
            for pg in pg_reqs
        ]
        if not concerns:
            recommendations.append("No major compatibility concerns detected")

        return {
            "feasibility_score": score,
            "overall_confidence": confidence,
            "functional_group_analysis": {
                "precursor_groups": precursor_groups,
                "compatibility_issues": compat_issues,
            },
            "side_reaction_risks": side_risks,
            "protecting_group_requirements": pg_reqs,
            "concerns": concerns,
            "recommendations": recommendations,
            "success": True,
            "error": None,
        }
    except Exception as e:
        empty["error"] = f"Feasibility assessment failed: {e}"
        return empty
