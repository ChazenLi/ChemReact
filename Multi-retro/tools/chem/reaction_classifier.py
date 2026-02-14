"""
Reaction Classification Module
==============================
SMARTS pattern matching + keyword fallback for reaction type classification.
Migrated from internal/reaction_classifier/ (classifier.py + rules.py + rules.json).

The rules.json file remains at its original location and is loaded via
``internal.reaction_classifier.rules.load_rules``.

Public API:
    classify_reaction(left_side, right_side, reaction_class, rules_path=None) -> Dict
    infer_missing_byproducts_for_delta(reaction_class, delta, ...) -> List[Dict]
"""

from __future__ import annotations

from collections import Counter
from typing import Any, Dict, List, Optional, Tuple

from rdkit import Chem

from tools.common.constants import CONFIDENCE_THRESHOLD

# Re-use the existing rules loader (reads rules.json)
from tools.legacy.reaction_classifier.rules import load_rules


def _mols_from_side(side: str) -> Tuple[List[Any], bool]:
    mols = []
    for frag in side.split("."):
        frag = frag.strip()
        if not frag:
            continue
        mol = Chem.MolFromSmiles(frag)
        if mol is None:
            return [], False
        mols.append(mol)
    return mols, True


def _has_smarts(mols: List[Any], smarts: str) -> bool:
    try:
        patt = Chem.MolFromSmarts(smarts)
    except Exception:
        return False
    if patt is None:
        return False
    return any(mol.HasSubstructMatch(patt) for mol in mols)


def _keyword_fallback_candidates(reaction_class: str) -> List[Dict[str, Any]]:
    name = (reaction_class or "").lower()
    fallback: List[Dict[str, Any]] = []

    if any(k in name for k in ["ester", "amid", "dehydr", "condensation"]):
        fallback.append({"smiles": "O", "confidence": 0.45, "source": "keyword"})
    if any(k in name for k in ["acid chloride", "acyl chloride", "schotten", "amid"]):
        fallback.append({"smiles": "Cl", "confidence": 0.4, "source": "keyword"})
    if any(k in name for k in ["aminolysis", "transester"]):
        fallback.extend([
            {"smiles": "CO", "confidence": 0.35, "source": "keyword"},
            {"smiles": "CCO", "confidence": 0.35, "source": "keyword"},
        ])
    if any(k in name for k in ["acetyl", "acylation"]):
        fallback.append({"smiles": "CC(=O)O", "confidence": 0.3, "source": "keyword"})

    dedup: Dict[str, Dict[str, Any]] = {}
    for item in fallback:
        s = item["smiles"]
        if s not in dedup or item["confidence"] > dedup[s]["confidence"]:
            dedup[s] = item
    return list(dedup.values())


def _merge_byproduct_candidates(*candidate_lists: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    merged: Dict[str, Dict[str, Any]] = {}
    for candidates in candidate_lists:
        for item in candidates:
            smiles = item.get("smiles")
            confidence = item.get("confidence")
            source = item.get("source", "unknown")
            if not isinstance(smiles, str) or not smiles:
                continue
            if not isinstance(confidence, (int, float)):
                continue
            conf = float(max(0.0, min(1.0, confidence)))
            if smiles not in merged or conf > merged[smiles]["confidence"]:
                merged[smiles] = {"smiles": smiles, "confidence": round(conf, 3), "source": source}
    return sorted(merged.values(), key=lambda x: x["confidence"], reverse=True)


def _formula_counter_from_smiles(smiles: str) -> Optional[Counter]:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    mol = Chem.AddHs(mol)
    counts: Counter = Counter()
    for atom in mol.GetAtoms():
        counts[atom.GetSymbol()] += 1
    return counts


def classify_reaction(
    left_side: str,
    right_side: str,
    reaction_class: str,
    rules_path: Optional[str] = None,
) -> Dict[str, Any]:
    """Classify a reaction using SMARTS pattern matching + keyword fallback.

    Scoring: 55% reactant SMARTS + 35% product SMARTS + 10% keyword match.

    Returns dict with id, name, category, confidence, matched_features,
    candidate_byproducts, byproduct_candidates.
    """
    left_mols, left_ok = _mols_from_side(left_side)
    right_mols, right_ok = _mols_from_side(right_side)
    reaction_name = (reaction_class or "").lower()

    fallback_candidates = _keyword_fallback_candidates(reaction_class)
    best: Dict[str, Any] = {
        "id": "unknown",
        "name": "Unknown",
        "category": "uncategorized",
        "confidence": 0.0,
        "source": "none",
        "matched_features": [],
        "candidate_byproducts": [c["smiles"] for c in fallback_candidates],
        "byproduct_candidates": fallback_candidates,
    }

    if not left_ok or not right_ok:
        return best

    rules = load_rules(rules_path)
    for rule in rules:
        reactant_smarts = rule.get("reactant_smarts", [])
        product_smarts = rule.get("product_smarts", [])

        reactant_hits = [p for p in reactant_smarts if _has_smarts(left_mols, p)]
        product_hits = [p for p in product_smarts if _has_smarts(right_mols, p)]
        keyword_hit = any(k in reaction_name for k in rule.get("keywords", []))

        r_score = (len(reactant_hits) / max(1, len(reactant_smarts))) if reactant_smarts else 0.0
        p_score = (len(product_hits) / max(1, len(product_smarts))) if product_smarts else 0.0
        k_score = 1.0 if keyword_hit else 0.0
        score = 0.55 * r_score + 0.35 * p_score + 0.1 * k_score

        if score <= best["confidence"]:
            continue

        base_byproducts: List[Dict[str, Any]] = []
        for item in rule.get("byproducts", []):
            smiles = item.get("smiles")
            weight = item.get("weight", 1.0)
            if not isinstance(smiles, str) or not smiles.strip():
                continue
            if not isinstance(weight, (int, float)):
                weight = 1.0
            conf = float(max(0.0, min(1.0, score * float(weight))))
            base_byproducts.append({
                "smiles": smiles.strip(),
                "confidence": round(conf, 3),
                "source": f"rule:{rule['id']}",
            })

        byproduct_candidates = _merge_byproduct_candidates(base_byproducts, fallback_candidates)

        matched_features: List[str] = []
        matched_features.extend([f"reactant:{p}" for p in reactant_hits])
        matched_features.extend([f"product:{p}" for p in product_hits])
        if keyword_hit:
            matched_features.append("keyword")

        best = {
            "id": rule["id"],
            "name": rule["name"],
            "category": str(rule.get("category") or "uncategorized"),
            "confidence": round(min(1.0, score), 3),
            "source": "smarts+keyword" if keyword_hit else "smarts",
            "matched_features": matched_features,
            "candidate_byproducts": [c["smiles"] for c in byproduct_candidates],
            "byproduct_candidates": byproduct_candidates,
        }

    if best["id"] == "unknown" and fallback_candidates:
        best["source"] = "keyword_fallback"

    return best


def infer_missing_byproducts_for_delta(
    reaction_class: str,
    delta: Counter,
    reaction_type_info: Optional[Dict[str, Any]] = None,
    max_byproducts: int = 2,
) -> List[Dict[str, Any]]:
    """Infer missing byproducts that would explain an atom-count delta.

    Args:
        reaction_class: Reaction class name for keyword fallback.
        delta: Counter of excess atoms (all values must be >= 0).
        reaction_type_info: Optional classify_reaction() output for typed candidates.
        max_byproducts: Max number of byproducts to combine (1 or 2).

    Returns:
        List of {smiles, confidence, source} dicts.
    """
    if not delta or not all(v >= 0 for v in delta.values()):
        return []

    typed_candidates: List[Dict[str, Any]] = []
    if isinstance(reaction_type_info, dict):
        for item in reaction_type_info.get("byproduct_candidates", []):
            if not isinstance(item, dict):
                continue
            smiles = item.get("smiles")
            confidence = item.get("confidence")
            source = item.get("source", "reaction_type")
            if isinstance(smiles, str) and isinstance(confidence, (int, float)):
                typed_candidates.append({
                    "smiles": smiles, "confidence": float(confidence), "source": source,
                })

    fallback_candidates = _keyword_fallback_candidates(reaction_class)
    candidates = _merge_byproduct_candidates(typed_candidates, fallback_candidates)
    if not candidates:
        return []

    candidate_counters: List[Tuple[Dict[str, Any], Counter]] = []
    for item in candidates:
        cnt = _formula_counter_from_smiles(item["smiles"])
        if cnt is not None:
            candidate_counters.append((item, cnt))
    if not candidate_counters:
        return []

    best_combo: List[Dict[str, Any]] = []
    best_score = -1.0

    for item, cnt in candidate_counters:
        if cnt == delta and item["confidence"] > best_score:
            best_combo = [item]
            best_score = item["confidence"]

    if max_byproducts >= 2:
        for item1, cnt1 in candidate_counters:
            for item2, cnt2 in candidate_counters:
                if cnt1 + cnt2 != delta:
                    continue
                score = float(item1["confidence"]) + float(item2["confidence"])
                if score > best_score:
                    best_combo = [item1, item2]
                    best_score = score

    return [
        {
            "smiles": item["smiles"],
            "confidence": round(float(item["confidence"]), 3),
            "source": item.get("source", "unknown"),
        }
        for item in best_combo
    ]
