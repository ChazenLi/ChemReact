from collections import Counter
import importlib.util
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

from rdkit import Chem

try:
    from .rules import load_rules
except ImportError:
    _RULES_PATH = Path(__file__).with_name("rules.py")
    _RULES_SPEC = importlib.util.spec_from_file_location("reaction_classifier_rules", str(_RULES_PATH))
    if _RULES_SPEC is None or _RULES_SPEC.loader is None:
        raise RuntimeError(f"Cannot load reaction rules from {_RULES_PATH}")
    _RULES_MOD = importlib.util.module_from_spec(_RULES_SPEC)
    _RULES_SPEC.loader.exec_module(_RULES_MOD)
    load_rules = _RULES_MOD.load_rules


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


def _compile_smarts(pattern: str):
    try:
        return Chem.MolFromSmarts(pattern)
    except Exception:
        return None


def _has_smarts(mols: List[Any], smarts: str) -> bool:
    patt = _compile_smarts(smarts)
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
        fallback.extend(
            [
                {"smiles": "CO", "confidence": 0.35, "source": "keyword"},
                {"smiles": "CCO", "confidence": 0.35, "source": "keyword"},
            ]
        )
    if any(k in name for k in ["acetyl", "acylation"]):
        fallback.append({"smiles": "CC(=O)O", "confidence": 0.3, "source": "keyword"})

    dedup: Dict[str, Dict[str, Any]] = {}
    for item in fallback:
        smiles = item["smiles"]
        if smiles not in dedup or item["confidence"] > dedup[smiles]["confidence"]:
            dedup[smiles] = item
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
    left_mols, left_ok = _mols_from_side(left_side)
    right_mols, right_ok = _mols_from_side(right_side)
    reaction_name = (reaction_class or "").lower()

    fallback_candidates = _keyword_fallback_candidates(reaction_class)
    best: Dict[str, Any] = {
        "id": "unknown",
        "name": "Unknown",
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

        reactant_score = (len(reactant_hits) / max(1, len(reactant_smarts))) if reactant_smarts else 0.0
        product_score = (len(product_hits) / max(1, len(product_smarts))) if product_smarts else 0.0
        keyword_score = 1.0 if keyword_hit else 0.0
        score = (0.55 * reactant_score) + (0.35 * product_score) + (0.1 * keyword_score)

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
            base_byproducts.append(
                {"smiles": smiles.strip(), "confidence": round(conf, 3), "source": f"rule:{rule['id']}"}
            )

        byproduct_candidates = _merge_byproduct_candidates(base_byproducts, fallback_candidates)

        matched_features: List[str] = []
        matched_features.extend([f"reactant:{p}" for p in reactant_hits])
        matched_features.extend([f"product:{p}" for p in product_hits])
        if keyword_hit:
            matched_features.append("keyword")

        best = {
            "id": rule["id"],
            "name": rule["name"],
            "confidence": round(min(1.0, score), 3),
            "source": "smarts+keyword" if keyword_hit else "smarts",
            "matched_features": matched_features,
            "candidate_byproducts": [c["smiles"] for c in byproduct_candidates],
            "byproduct_candidates": byproduct_candidates,
        }

    if best["id"] == "unknown" and fallback_candidates:
        best["source"] = "keyword_fallback"

    return best


def _counter_is_nonnegative(delta: Counter) -> bool:
    return all(v >= 0 for v in delta.values())


def infer_missing_byproducts_for_delta(
    reaction_class: str,
    delta: Counter,
    reaction_type_info: Optional[Dict[str, Any]] = None,
    max_byproducts: int = 2,
) -> List[Dict[str, Any]]:
    if not delta or not _counter_is_nonnegative(delta):
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
                typed_candidates.append({"smiles": smiles, "confidence": float(confidence), "source": source})

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
