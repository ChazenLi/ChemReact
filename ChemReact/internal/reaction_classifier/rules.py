import json
import os
from pathlib import Path
from typing import Any, Dict, List, Optional


_DEFAULT_RULES_PATH = Path(__file__).with_name("rules.json")


def _normalize_byproducts(values: Any) -> List[Dict[str, Any]]:
    if not isinstance(values, list):
        return []
    normalized: List[Dict[str, Any]] = []
    for entry in values:
        if isinstance(entry, str):
            normalized.append({"smiles": entry, "weight": 1.0})
            continue
        if not isinstance(entry, dict):
            continue
        smiles = entry.get("smiles")
        if not isinstance(smiles, str) or not smiles.strip():
            continue
        weight = entry.get("weight", 1.0)
        if not isinstance(weight, (int, float)):
            weight = 1.0
        normalized.append({"smiles": smiles.strip(), "weight": float(max(0.0, min(1.0, weight)))})
    return normalized


def _normalize_rule(rule: Any) -> Optional[Dict[str, Any]]:
    if not isinstance(rule, dict):
        return None
    rule_id = rule.get("id")
    name = rule.get("name")
    if not isinstance(rule_id, str) or not rule_id.strip():
        return None
    if not isinstance(name, str) or not name.strip():
        return None

    def _string_list(key: str) -> List[str]:
        values = rule.get(key, [])
        if not isinstance(values, list):
            return []
        return [str(v).strip() for v in values if isinstance(v, str) and str(v).strip()]

    return {
        "id": rule_id.strip(),
        "name": name.strip(),
        "reactant_smarts": _string_list("reactant_smarts"),
        "product_smarts": _string_list("product_smarts"),
        "keywords": [k.lower() for k in _string_list("keywords")],
        "byproducts": _normalize_byproducts(rule.get("byproducts", [])),
    }


def load_rules(path: Optional[str] = None) -> List[Dict[str, Any]]:
    rule_path = path or os.getenv("CHEMREACT_REACTION_RULES_PATH")
    resolved_path = Path(rule_path) if rule_path else _DEFAULT_RULES_PATH

    if not resolved_path.exists():
        return []

    try:
        payload = json.loads(resolved_path.read_text(encoding="utf-8-sig"))
    except Exception:
        return []

    if not isinstance(payload, list):
        return []

    rules: List[Dict[str, Any]] = []
    for raw in payload:
        normalized = _normalize_rule(raw)
        if normalized is not None:
            rules.append(normalized)
    return rules
