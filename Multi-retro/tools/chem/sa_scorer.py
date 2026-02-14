"""
SA Score Estimation Module
==========================
Consolidates SA scoring logic from executor.py and availability.py
into a single module using SAThresholds from the configuration registry.

Public API:
    estimate_sa_score(smiles) -> float
    classify_sa(smiles)       -> Dict[str, Any]
"""

from __future__ import annotations

import logging
from typing import Any, Dict

from tools.chem.mol_parser import parse_smiles
from tools.common.constants import SAThresholds

logger = logging.getLogger(__name__)


def estimate_sa_score(smiles: str) -> float:
    """Estimate SA Score for a SMILES string.

    Uses RDKit SA_Score (sascorer) when available, otherwise falls back
    to an improved heuristic considering ring count, heteroatom count,
    chirality, molecular weight, and SMILES complexity.

    Returns:
        Float in [1.0, 10.0]. Lower = easier to synthesize.
    """
    try:
        from rdkit.Contrib.SA_Score import sascorer
        mol = parse_smiles(smiles)
        if mol is not None:
            return float(sascorer.calculateScore(mol))
    except Exception:
        pass
    return _heuristic_sa_score(smiles)


def _heuristic_sa_score(smiles: str) -> float:
    """Fallback heuristic SA score when RDKit sascorer is unavailable.

    Considers ring count, heteroatom count, chirality indicators,
    molecular weight, and SMILES length as complexity proxy.
    """
    # Ring count from SMILES ring-closure digits and %nn notation
    ring_count = sum(1 for ch in smiles if ch.isdigit() and ch != "0") // 2
    ring_count += smiles.count("%")

    # Heteroatom count
    upper = smiles.upper()
    heteroatom_count = sum(upper.count(a) for a in ("N", "O", "S", "P", "F"))
    heteroatom_count += smiles.count("Cl") + smiles.count("Br") + upper.count("I")

    # Stereocenters / stereochemistry
    stereo_count = smiles.count("@") + smiles.count("/") + smiles.count("\\")

    # Molecular weight contribution
    mw_penalty = 0.0
    mol = parse_smiles(smiles)
    if mol is not None:
        try:
            from rdkit.Chem import Descriptors
            mw = Descriptors.MolWt(mol)
            if mw > 500:
                mw_penalty = 0.8
            elif mw > 300:
                mw_penalty = 0.3
        except Exception:
            pass

    score = 2.5
    score += ring_count * 0.5
    score += heteroatom_count * 0.3
    score += stereo_count * 0.8
    score += mw_penalty
    if len(smiles) > 40:
        score += 0.5

    return max(1.0, min(score, 10.0))


def classify_sa(smiles: str) -> Dict[str, Any]:
    """Estimate SA score and classify using SAThresholds.

    Returns:
        {
            "smiles": str,
            "sa_score": float,
            "classification": str,  # "purchasable" | "easily_synthesizable" | "complex"
            "is_terminal": bool,    # True if purchasable
        }
    """
    sa_score = estimate_sa_score(smiles)
    classification = SAThresholds.classify(sa_score)
    return {
        "smiles": smiles,
        "sa_score": round(sa_score, 3),
        "classification": classification,
        "is_terminal": classification == "purchasable",
    }
