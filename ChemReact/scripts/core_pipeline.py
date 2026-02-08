import argparse
import json
import sys
import time
from pathlib import Path
import importlib.util
from typing import Any, Dict, List, Tuple, Optional
from render_report import render_and_report
from planner_audit import (
    audit_routes_with_gate,
    build_planner_prompt,
    build_planner_request,
    invoke_host_planner,
    planner_request_payload,
    render_host_command,
    run_host_planner_repair_loop,
    resolve_planner_routes_candidate,
)


VIS_PLAN_TEMPLATE: Dict[str, Any] = {
    "selected_route_ids": [1],
    "target_image": {
        "legend": "Target with key motif",
        "highlight_atoms": []
    },
    "routes": [
        {
            "route_id": 1,
            "why_selected": "Best score/risk balance",
            "overview_grid": {
                "enabled": True,
                "precursor_legends": ["Fragment A", "Fragment B"]
            },
            "tree_view": {
                "enabled": True,
                "caption": "Main strategic disconnection"
            },
            "step_views": [
                {
                    "step_index": 1,
                    "caption": "Key bond-forming step"
                }
            ],
            "focus_molecules": []
        }
    ],
    "global_notes": [
        "Use concise captions.",
        "Mark high-risk steps clearly."
    ]
}

PLANNER_ROUTES_TEMPLATE: List[Dict[str, Any]] = [
    {
        "route_id": 1,
        "synthesis_style": "convergent",
        "score": 5.0,
        "audit_verdict": "CONDITIONAL",
        "critical_issues": ["Planner output requires validation."],
        "precursors": ["<SMILES_1>", "<SMILES_2>"],
        "steps": [
            {
                "reaction_class": "Proposed Reaction Class",
                "reagents": ["Reagent A"],
                "conditions": "TBD",
                "reaction_smiles": "<REACTION_SMILES>"
            }
        ]
    }
]

# Route portfolio constraint: minimum 3 routes, no hard upper bound in gate logic.
MIN_REQUIRED_ROUTES = 3
MAX_AUTO_PROPOSE_ROUTES = 5
MAX_PLANNER_ITERATIONS = 5


def _load_module(module_name: str, file_path: Path):
    spec = importlib.util.spec_from_file_location(module_name, str(file_path))
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot load module from {file_path}")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _is_non_production_host_command(cmd: str) -> bool:
    text = (cmd or "").strip().lower()
    if not text:
        return True
    blocked_markers = [
        "your_host_llm_runner.py",
        "host_planner_fallback.py",
        "__set_",
        "__replace_",
        "<your_",
    ]
    return any(marker in text for marker in blocked_markers)


_NORMALIZE_MOD = None
_ANALYSIS_MOD = None
_PROMPTS_PERSONAS_MOD = None
_REACTION_RULES_CACHE: Optional[List[Dict[str, Any]]] = None


def _get_normalize_mod():
    """Lazy load the SMILES normalization module."""
    global _NORMALIZE_MOD
    if _NORMALIZE_MOD is None:
        skill_root = Path(__file__).resolve().parents[1]
        normalize_path = skill_root / "internal" / "normalize" / "smiles_normalizer.py"
        _NORMALIZE_MOD = _load_module("smiles_normalizer_loop", normalize_path)
    return _NORMALIZE_MOD


def _get_analysis_mod():
    """Lazy load the analysis module (precursor, reaction, repair guidance)."""
    global _ANALYSIS_MOD
    if _ANALYSIS_MOD is None:
        skill_root = Path(__file__).resolve().parents[1]
        repair_guidance_path = skill_root / "internal" / "analysis" / "repair_guidance.py"
        _ANALYSIS_MOD = _load_module("repair_guidance_loop", repair_guidance_path)
    return _ANALYSIS_MOD


def _get_prompts_personas_mod():
    """Lazy load persona prompt builders for repair-specialist integration."""
    global _PROMPTS_PERSONAS_MOD
    if _PROMPTS_PERSONAS_MOD is None:
        skill_root = Path(__file__).resolve().parents[1]
        prompts_path = skill_root / "internal" / "retroskill" / "prompts_personas.py"
        _PROMPTS_PERSONAS_MOD = _load_module("prompts_personas_loop", prompts_path)
    return _PROMPTS_PERSONAS_MOD


def _read_json(path: str, default: Any) -> Any:
    if not path:
        return default
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)


def _write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2, ensure_ascii=False)


def _validate_strategy(strategy: Any) -> List[str]:
    errs: List[str] = []
    if not isinstance(strategy, dict):
        return ["strategy must be a JSON object"]
    analysis = strategy.get("analysis", {})
    if analysis and not isinstance(analysis, dict):
        errs.append("strategy.analysis must be an object")
    if isinstance(analysis, dict):
        if "complexity_features" in analysis and not isinstance(analysis["complexity_features"], list):
            errs.append("strategy.analysis.complexity_features must be an array")
        if "key_disconnections" in analysis and not isinstance(analysis["key_disconnections"], list):
            errs.append("strategy.analysis.key_disconnections must be an array")
        if "global_constraints" in analysis and not isinstance(analysis["global_constraints"], list):
            errs.append("strategy.analysis.global_constraints must be an array")
        if "risk_budget" in analysis and not isinstance(analysis["risk_budget"], dict):
            errs.append("strategy.analysis.risk_budget must be an object")
        risk_budget = analysis.get("risk_budget", {})
        if isinstance(risk_budget, dict):
            for key in ["max_hard_issues", "max_soft_issues"]:
                if key in risk_budget and not isinstance(risk_budget[key], int):
                    errs.append(f"strategy.analysis.risk_budget.{key} must be an integer")
            if "min_quality_score" in risk_budget and not isinstance(risk_budget["min_quality_score"], (int, float)):
                errs.append("strategy.analysis.risk_budget.min_quality_score must be a number")
    return errs


def _validate_steps(steps: Any, route_idx: int) -> List[str]:
    errs: List[str] = []
    if not isinstance(steps, list):
        return [f"routes[{route_idx}].steps must be an array"]
    for j, step in enumerate(steps):
        if not isinstance(step, dict):
            errs.append(f"routes[{route_idx}].steps[{j}] must be an object")
            continue
        if "reaction_class" in step and not isinstance(step["reaction_class"], str):
            errs.append(f"routes[{route_idx}].steps[{j}].reaction_class must be a string")
        if "reagents" in step and not isinstance(step["reagents"], list):
            errs.append(f"routes[{route_idx}].steps[{j}].reagents must be an array")
    return errs


def _validate_routes(routes: Any) -> List[str]:
    errs: List[str] = []
    if not isinstance(routes, list):
        return ["routes must be a JSON array"]
    for i, route in enumerate(routes):
        if not isinstance(route, dict):
            errs.append(f"routes[{i}] must be an object")
            continue
        for req in ["route_id", "score", "audit_verdict", "critical_issues", "steps"]:
            if req not in route:
                errs.append(f"routes[{i}].{req} is required")
        if "route_id" in route and not isinstance(route["route_id"], int):
            errs.append(f"routes[{i}].route_id must be an integer")
        if "score" in route and not isinstance(route["score"], (int, float)):
            errs.append(f"routes[{i}].score must be a number")
        if "audit_verdict" in route and not isinstance(route["audit_verdict"], str):
            errs.append(f"routes[{i}].audit_verdict must be a string")
        if "synthesis_style" in route:
            if not isinstance(route["synthesis_style"], str):
                errs.append(f"routes[{i}].synthesis_style must be a string")
            elif route["synthesis_style"].lower() not in {"linear", "convergent"}:
                errs.append(f"routes[{i}].synthesis_style must be linear or convergent")
        if "critical_issues" in route and not isinstance(route["critical_issues"], list):
            errs.append(f"routes[{i}].critical_issues must be an array")
        if "precursors" in route and not isinstance(route["precursors"], list):
            errs.append(f"routes[{i}].precursors must be an array")
        if "steps" in route:
            errs.extend(_validate_steps(route["steps"], i))
    return errs


def _coerce_routes_payload(payload: Any) -> Any:
    if isinstance(payload, list):
        return payload
    if isinstance(payload, dict) and isinstance(payload.get("routes"), list):
        return payload["routes"]
    return payload


def _validate_vis_plan(vis_plan: Any, valid_route_ids: List[int]) -> List[str]:
    errs: List[str] = []
    if vis_plan in ({}, None):
        return errs
    if not isinstance(vis_plan, dict):
        return ["vis_plan must be a JSON object"]

    selected = vis_plan.get("selected_route_ids", [])
    if not isinstance(selected, list):
        errs.append("vis_plan.selected_route_ids must be an array")
    else:
        for i, rid in enumerate(selected):
            if not isinstance(rid, int):
                errs.append(f"vis_plan.selected_route_ids[{i}] must be an integer")
            elif valid_route_ids and rid not in valid_route_ids:
                errs.append(f"vis_plan.selected_route_ids[{i}] not found in routes route_id set")

    target_image = vis_plan.get("target_image", {})
    if target_image and not isinstance(target_image, dict):
        errs.append("vis_plan.target_image must be an object")
    if isinstance(target_image, dict):
        if "legend" in target_image and not isinstance(target_image["legend"], str):
            errs.append("vis_plan.target_image.legend must be a string")
        if "highlight_atoms" in target_image and not isinstance(target_image["highlight_atoms"], list):
            errs.append("vis_plan.target_image.highlight_atoms must be an array")

    routes = vis_plan.get("routes", [])
    if routes and not isinstance(routes, list):
        errs.append("vis_plan.routes must be an array")
    if isinstance(routes, list):
        for i, route in enumerate(routes):
            if not isinstance(route, dict):
                errs.append(f"vis_plan.routes[{i}] must be an object")
                continue
            rid = route.get("route_id")
            if not isinstance(rid, int):
                errs.append(f"vis_plan.routes[{i}].route_id must be an integer")
            elif valid_route_ids and rid not in valid_route_ids:
                errs.append(f"vis_plan.routes[{i}].route_id not found in routes route_id set")
            step_views = route.get("step_views", [])
            if step_views and not isinstance(step_views, list):
                errs.append(f"vis_plan.routes[{i}].step_views must be an array")
            if isinstance(step_views, list):
                for j, sv in enumerate(step_views):
                    if not isinstance(sv, dict):
                        errs.append(f"vis_plan.routes[{i}].step_views[{j}] must be an object")
                        continue
                    if not isinstance(sv.get("step_index"), int):
                        errs.append(f"vis_plan.routes[{i}].step_views[{j}].step_index must be an integer")
                    if "caption" in sv and not isinstance(sv["caption"], str):
                        errs.append(f"vis_plan.routes[{i}].step_views[{j}].caption must be a string")

    return errs


def _safe_int_list(values: Any, atom_count: int) -> List[int]:
    if not isinstance(values, list):
        return []
    out: List[int] = []
    for v in values:
        if isinstance(v, int) and 0 <= v < atom_count:
            out.append(v)
    return out


def _pick_routes(all_routes: List[Dict[str, Any]], vis_plan: Dict[str, Any], top_k: int) -> List[Dict[str, Any]]:
    if not all_routes:
        return []

    selected_ids = vis_plan.get("selected_route_ids", []) if isinstance(vis_plan, dict) else []
    if isinstance(selected_ids, list) and selected_ids:
        selected = []
        id_set = set(selected_ids)
        for r in all_routes:
            if r.get("route_id") in id_set:
                selected.append(r)
        if selected:
            return selected[:3]

    ranked = sorted(
        all_routes,
        key=lambda x: (x.get("quality_score", x.get("score", 0)), x.get("score", 0)),
        reverse=True,
    )
    return ranked[: max(1, min(3, top_k))]


def _get_route_plan(vis_plan: Dict[str, Any], route_id: Any) -> Dict[str, Any]:
    for item in vis_plan.get("routes", []):
        if item.get("route_id") == route_id:
            return item
    return {}


def _step_caption_map(route_plan: Dict[str, Any]) -> Dict[int, str]:
    out: Dict[int, str] = {}
    for item in route_plan.get("step_views", []):
        idx = item.get("step_index")
        cap = item.get("caption")
        if isinstance(idx, int) and isinstance(cap, str):
            out[idx] = cap
    return out


def _canonicalize_smiles(smiles: str):
    from rdkit import Chem

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    for atom in mol.GetAtoms():
        if atom.GetAtomMapNum():
            atom.SetAtomMapNum(0)
    return Chem.MolToSmiles(mol, canonical=True)


# NOTE: legacy step-level audit helper functions removed from core_pipeline.

def _issue_code(message: str) -> str:
    text = (message or "").lower()
    if "acid_base_coexist" in text or "carboxylic acid and amine coexist" in text:
        return "ACID_BASE_COEXIST"
    if "competing_nucleophiles" in text or "multiple nucleophiles detected" in text:
        return "COMPETING_NUCLEOPHILES"
    if "atom conservation failed" in text:
        return "ATOM_CONSERVATION_FAIL"
    if "atom conservation satisfied after inferred byproducts" in text:
        return "ATOM_BALANCE_INFERRED_BYPRODUCT"
    if "reaction must contain '>>'" in text or "empty reactant or product side" in text:
        return "REACTION_FORMAT_INVALID"
    if "invalid molecule in reaction side" in text:
        return "REACTION_SMILES_INVALID"
    if "mapped reaction is invalid" in text:
        return "MAPPING_INVALID"
    if "mapping completeness is low" in text:
        return "MAPPING_INCOMPLETE"
    if "no bond changes detected" in text:
        return "NO_BOND_CHANGE_DETECTED"
    if "products do not connect" in text:
        return "ROUTE_COHERENCE_FAIL"
    if "final step product does not match target molecule" in text:
        return "TARGET_MISMATCH"
    if "target molecule analysis unavailable" in text:
        return "TARGET_ANALYSIS_MISSING"
    if "route precursors missing" in text:
        return "PRECURSOR_SET_MISSING"
    if "precursor analysis failed" in text:
        return "PRECURSOR_ANALYSIS_FAIL"
    if "functional-group analysis missing" in text:
        return "FG_ANALYSIS_MISSING"
    if "molecule analysis failed" in text:
        return "MOLECULE_ANALYSIS_FAIL"
    if "reaction class suggests" in text:
        return "REACTION_CLASS_INCONSISTENT"
    if "ring-count jump" in text:
        return "RING_TRANSFORMATION_INCONSISTENT"
    if "cannot explain target ring topology" in text:
        return "RING_TOPOLOGY_MISMATCH"
    if "heterocycle topology" in text:
        return "HETEROCYCLE_TOPOLOGY_MISMATCH"
    if "scaffold does not preserve target core skeleton" in text:
        return "SCAFFOLD_CONTINUITY_FAIL"
    if "do not align with strategy key_disconnections" in text:
        return "STRATEGY_MISMATCH"
    if "risk_budget" in text:
        return "RISK_BUDGET_EXCEEDED"
    if "quality below strategy risk budget threshold" in text:
        return "QUALITY_BELOW_THRESHOLD"
    if "route portfolio" in text or "style mix mismatch" in text:
        return "PORTFOLIO_STYLE_MISMATCH"
    return "AUDIT_ISSUE"


_ISSUE_TO_TEMPLATE_IDS: Dict[str, List[str]] = {
    "NO_BOND_CHANGE_DETECTED": ["suzuki_coupling", "buchwald_hartwig", "sn2_substitution"],
    "REACTION_CLASS_INCONSISTENT": ["suzuki_coupling", "buchwald_hartwig", "amidation"],
    "RING_TRANSFORMATION_INCONSISTENT": ["suzuki_coupling", "buchwald_hartwig"],
    "REACTION_SMILES_INVALID": ["suzuki_coupling", "buchwald_hartwig", "sn2_substitution"],
    "ATOM_CONSERVATION_FAIL": ["suzuki_coupling", "amidation", "esterification"],
    "ATOM_BALANCE_INFERRED_BYPRODUCT": ["amidation", "esterification", "suzuki_coupling"],
}


def _load_reaction_rules() -> List[Dict[str, Any]]:
    global _REACTION_RULES_CACHE
    if _REACTION_RULES_CACHE is not None:
        return _REACTION_RULES_CACHE
    skill_root = Path(__file__).resolve().parents[1]
    rules_path = skill_root / "internal" / "reaction_classifier" / "rules.json"
    try:
        payload = json.loads(rules_path.read_text(encoding="utf-8-sig"))
        if isinstance(payload, list):
            _REACTION_RULES_CACHE = [item for item in payload if isinstance(item, dict)]
        else:
            _REACTION_RULES_CACHE = []
    except Exception:
        _REACTION_RULES_CACHE = []
    return _REACTION_RULES_CACHE


def _reaction_rule_by_id(rule_id: str) -> Optional[Dict[str, Any]]:
    if not rule_id:
        return None
    for item in _load_reaction_rules():
        if str(item.get("id", "")).strip().lower() == rule_id.strip().lower():
            return item
    return None


def _repair_template_hints_for_step(step: Dict[str, Any]) -> List[str]:
    hints: List[str] = []
    candidate_ids: List[str] = []

    reaction_type = step.get("reaction_type", {})
    if isinstance(reaction_type, dict):
        rid = reaction_type.get("id")
        if isinstance(rid, str) and rid and rid.lower() != "unknown":
            candidate_ids.append(rid.lower())

    issue_codes: List[str] = []
    for msg in step.get("hard_issues", []) if isinstance(step.get("hard_issues", []), list) else []:
        if isinstance(msg, str):
            issue_codes.append(_issue_code(msg))
    for msg in step.get("soft_issues", []) if isinstance(step.get("soft_issues", []), list) else []:
        if isinstance(msg, str):
            issue_codes.append(_issue_code(msg))
    for code in issue_codes:
        for item in _ISSUE_TO_TEMPLATE_IDS.get(code, []):
            candidate_ids.append(item)

    reaction_class_text = str(step.get("reaction_class", "")).lower()
    if reaction_class_text:
        for rule in _load_reaction_rules():
            keywords = rule.get("keywords", [])
            if not isinstance(keywords, list):
                continue
            if any(isinstance(k, str) and k.lower() in reaction_class_text for k in keywords):
                rid = str(rule.get("id", "")).lower()
                if rid:
                    candidate_ids.append(rid)

    dedup_ids: List[str] = []
    for rid in candidate_ids:
        if rid and rid not in dedup_ids:
            dedup_ids.append(rid)

    for rid in dedup_ids[:3]:
        rule = _reaction_rule_by_id(rid)
        if not rule:
            continue
        byproducts = []
        for bp in rule.get("byproducts", []) if isinstance(rule.get("byproducts", []), list) else []:
            if isinstance(bp, dict) and isinstance(bp.get("smiles"), str):
                byproducts.append(bp["smiles"])
        byp_text = ", ".join(byproducts[:3]) if byproducts else "none"
        hints.append(
            f"  - template={rule.get('id')} ({rule.get('name')}), likely_byproducts={byp_text}"
        )
    return hints


def _analyze_molecule_bundle(smiles: str, analyze_mod, features_mod, fg_mod) -> Dict[str, Any]:
    # Try to use normalized SMILES with display format
    try:
        normalize_mod = _get_normalize_mod()
        normalized = normalize_mod.normalize_smiles(smiles)
    except Exception:
        normalized = {
            "original": smiles,
            "canonical": _canonicalize_smiles(smiles) or smiles,
            "display": smiles
        }
    
    bundle: Dict[str, Any] = {
        "smiles": smiles,
        "smiles_normalized": normalized,  # New: normalized SMILES with display format
        "canonical_smiles": normalized.get("canonical") or _canonicalize_smiles(smiles),
        "display_smiles": normalized.get("display", smiles),  # New: for visualization
        "structure": {},
        "features": {},
        "functional_groups": {},
        "ring_profile": {},
        "success": False,
        "error": None,
    }
    
    # Use canonical SMILES for analysis to avoid mapping issues
    analysis_smiles = bundle["canonical_smiles"] or smiles
    
    structure = analyze_mod.analyze_smiles(analysis_smiles, generate_3d=False)
    features = features_mod.calculate_features(analysis_smiles)
    fg = fg_mod.analyze_functional_groups(analysis_smiles)
    ring_profile = _ring_hetero_profile(analysis_smiles)
    bundle["structure"] = structure
    bundle["features"] = features
    bundle["functional_groups"] = fg
    bundle["ring_profile"] = ring_profile
    ok = bool(structure.get("success")) and bool(features.get("success")) and bool(fg.get("success"))
    bundle["success"] = ok
    if not ok:
        errors: List[str] = []
        if not structure.get("success"):
            errors.append(f"structure={structure.get('error', 'unknown')}")
        if not features.get("success"):
            errors.append(f"features={features.get('error', 'unknown')}")
        if not fg.get("success"):
            errors.append(f"functional_groups={fg.get('error', 'unknown')}")
        bundle["error"] = "; ".join(errors)
    return bundle


def _ring_hetero_profile(smiles: str) -> Dict[str, Any]:
    from rdkit import Chem

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return {
            "success": False,
            "error": "invalid_smiles",
            "ring_count": 0,
            "ring_sizes": [],
            "heterocycle_count": 0,
            "aromatic_ring_count": 0,
        }

    ring_info = mol.GetRingInfo()
    atom_rings = list(ring_info.AtomRings())
    ring_sizes = sorted([len(r) for r in atom_rings])
    heterocycle_count = 0
    aromatic_ring_count = 0
    for ring in atom_rings:
        is_hetero = any(mol.GetAtomWithIdx(idx).GetSymbol() != "C" for idx in ring)
        is_aromatic = all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring)
        if is_hetero:
            heterocycle_count += 1
        if is_aromatic:
            aromatic_ring_count += 1

    return {
        "success": True,
        "error": None,
        "ring_count": len(atom_rings),
        "ring_sizes": ring_sizes,
        "heterocycle_count": heterocycle_count,
        "aromatic_ring_count": aromatic_ring_count,
    }


def _build_target_chemistry_context(target_smiles: str, analyze_mod, features_mod, fg_mod) -> Dict[str, Any]:
    bundle = _analyze_molecule_bundle(target_smiles, analyze_mod, features_mod, fg_mod)
    structure_props = bundle.get("structure", {}).get("properties", {})
    features_payload = bundle.get("features", {})
    fg_payload = bundle.get("functional_groups", {})
    synthesis_critical = features_payload.get("synthesis_critical", {})
    supplementary = features_payload.get("supplementary", {})
    topological = supplementary.get("topological_descriptors", {})
    chiral_centers = structure_props.get("chiral_centers", [])
    ring_profile = bundle.get("ring_profile", {})
    return {
        "success": bool(bundle.get("success")),
        "error": bundle.get("error"),
        "canonical_smiles": bundle.get("canonical_smiles"),
        "summary": {
            "formula": structure_props.get("formula"),
            "molecular_weight": structure_props.get("molecular_weight"),
            "ring_count": topological.get("Number_of_Rings"),
            "chiral_center_count": len(chiral_centers) if isinstance(chiral_centers, list) else 0,
            "tpsa": synthesis_critical.get("TPSA"),
            "logp": supplementary.get("LogP"),
        },
        "ring_profile": ring_profile,
        "functional_groups": fg_payload.get("functional_groups", {}),
        "alerts": fg_payload.get("alerts", []),
    }


# NOTE: legacy in-file audit implementation removed.
# Audit logic now lives in planner_audit.py (audit_routes_with_gate).

def _write_schema_files(schemas_dir: Path) -> None:
    routes_schema = {
        "$schema": "https://json-schema.org/draft/2020-12/schema",
        "title": "AuditedRoutes",
        "type": "array",
        "items": {
            "type": "object",
            "required": ["route_id", "score", "audit_verdict", "critical_issues", "steps"],
            "properties": {
                "route_id": {"type": "integer"},
                "synthesis_style": {"type": "string", "enum": ["linear", "convergent"]},
                "score": {"type": "number"},
                "audit_verdict": {"type": "string"},
                "critical_issues": {"type": "array", "items": {"type": "string"}},
                "precursors": {"type": "array", "items": {"type": "string"}},
                "steps": {
                    "type": "array",
                    "items": {
                        "type": "object",
                        "properties": {
                            "reaction_class": {"type": "string"},
                            "reagents": {"type": "array", "items": {"type": "string"}},
                            "conditions": {"type": "string"},
                            "reaction_smiles": {"type": "string"},
                            "rxn_smiles": {"type": "string"},
                            "smirks": {"type": "string"}
                        },
                        "additionalProperties": True
                    }
                }
            },
            "additionalProperties": True
        }
    }
    strategy_schema = {
        "$schema": "https://json-schema.org/draft/2020-12/schema",
        "title": "StrategyAnalysis",
        "type": "object",
        "properties": {
            "analysis": {
                "type": "object",
                "properties": {
                    "core_skeleton": {"type": "string"},
                    "complexity_features": {"type": "array", "items": {"type": "string"}},
                    "strategy_type": {"type": "string"},
                    "key_disconnections": {"type": "array", "items": {"type": "string"}},
                    "global_constraints": {"type": "array", "items": {"type": "string"}},
                    "risk_budget": {
                        "type": "object",
                        "properties": {
                            "max_hard_issues": {"type": "integer", "minimum": 0},
                            "max_soft_issues": {"type": "integer", "minimum": 0},
                            "min_quality_score": {"type": "number", "minimum": 0, "maximum": 10}
                        },
                        "additionalProperties": True
                    }
                },
                "additionalProperties": True
            }
        },
        "additionalProperties": True
    }
    vis_plan_schema = {
        "$schema": "https://json-schema.org/draft/2020-12/schema",
        "title": "VisualizationPlan",
        "type": "object",
        "properties": {
            "selected_route_ids": {"type": "array", "items": {"type": "integer"}},
            "target_image": {
                "type": "object",
                "properties": {
                    "legend": {"type": "string"},
                    "highlight_atoms": {"type": "array", "items": {"type": "integer"}}
                },
                "additionalProperties": True
            },
            "routes": {
                "type": "array",
                "items": {
                    "type": "object",
                    "required": ["route_id"],
                    "properties": {
                        "route_id": {"type": "integer"},
                        "why_selected": {"type": "string"},
                        "overview_grid": {
                            "type": "object",
                            "properties": {
                                "enabled": {"type": "boolean"},
                                "precursor_legends": {"type": "array", "items": {"type": "string"}}
                            },
                            "additionalProperties": True
                        },
                        "tree_view": {
                            "type": "object",
                            "properties": {
                                "enabled": {"type": "boolean"},
                                "caption": {"type": "string"}
                            },
                            "additionalProperties": True
                        },
                        "step_views": {
                            "type": "array",
                            "items": {
                                "type": "object",
                                "required": ["step_index"],
                                "properties": {
                                    "step_index": {"type": "integer"},
                                    "caption": {"type": "string"}
                                },
                                "additionalProperties": True
                            }
                        },
                        "focus_molecules": {"type": "array", "items": {"type": "object"}}
                    },
                    "additionalProperties": True
                }
            },
            "global_notes": {"type": "array", "items": {"type": "string"}}
        },
        "additionalProperties": True
    }
    planner_request_schema = {
        "$schema": "https://json-schema.org/draft/2020-12/schema",
        "title": "PlannerRequest",
        "type": "object",
        "required": ["target_smiles", "planner_backend", "max_routes", "analysis_first_required", "required_preplanning_tools"],
        "properties": {
            "target_smiles": {"type": "string"},
            "planner_backend": {"type": "string", "enum": ["host"]},
            "max_routes": {"type": "integer", "minimum": MIN_REQUIRED_ROUTES, "maximum": MAX_AUTO_PROPOSE_ROUTES},
            "analysis_first_required": {"type": "boolean", "const": True},
            "required_preplanning_tools": {
                "type": "array",
                "items": {"type": "string"},
                "minItems": 1,
            },
            "strategy": {"type": "object"},
            "target_chemistry": {"type": "object"},
            "constraints": {"type": "array", "items": {"type": "string"}}
        },
        "additionalProperties": True
    }
    planner_routes_schema = {
        "$schema": "https://json-schema.org/draft/2020-12/schema",
        "title": "PlannerRoutes",
        "type": "array",
        "items": routes_schema["items"]
    }
    _write_json(schemas_dir / "routes.schema.json", routes_schema)
    _write_json(schemas_dir / "strategy.schema.json", strategy_schema)
    _write_json(schemas_dir / "vis_plan.schema.json", vis_plan_schema)
    _write_json(schemas_dir / "planner_request.schema.json", planner_request_schema)
    _write_json(schemas_dir / "planner_routes.schema.json", planner_routes_schema)


def _build_planner_request(target_smiles: str, max_routes: int, strategy: Any, target_chemistry: Optional[Dict[str, Any]] = None) -> Dict[str, Any]:
    return build_planner_request(
        target_smiles=target_smiles,
        max_routes=max_routes,
        strategy=strategy,
        target_chemistry=target_chemistry,
    )


def _build_planner_prompt(planner_request: Dict[str, Any]) -> str:
    return build_planner_prompt(planner_request)


def _build_repair_prompt(repair_request: Dict[str, Any]) -> str:
    """
    Build repair prompt by combining:
    1) repair-guidance issue suggestions
    2) reaction-template hints from classifier rules
    3) dedicated Repair Specialist persona prompt
    """
    rejected_routes = repair_request.get("rejected_routes", [])
    if not isinstance(rejected_routes, list):
        rejected_routes = []

    repair_guidance_lines = []
    local_repair_lines = []
    template_hint_lines = []

    try:
        analysis_mod = _get_analysis_mod()
        for route in rejected_routes:
            if not isinstance(route, dict):
                continue
            route_id = route.get("route_id", "?")
            audit_details = route.get("audit_details", {})
            issues = []
            
            # Collect issues with codes
            for msg in audit_details.get("hard_issues", []):
                issues.append({"code": _issue_code(msg), "message": msg})
            for msg in audit_details.get("soft_issues", []):
                issues.append({"code": _issue_code(msg), "message": msg})
            
            if issues:
                guidance = analysis_mod.generate_repair_guidance(issues)
                repair_guidance_lines.append(f"\n### Route {route_id} Repair Guidance")
                repair_guidance_lines.append(f"Strategy: {guidance.get('repair_strategy', 'Review and fix all issues')}")
                
                for g in guidance.get("guidance", [])[:5]:  # Limit to top 5 issues
                    repair_guidance_lines.append(f"\n**{g.get('issue_code', 'UNKNOWN')}**: {g.get('issue_message', '')}")
                    for suggestion in g.get("suggestions", [])[:2]:  # Top 2 suggestions
                        repair_guidance_lines.append(f"  - {suggestion}")

            # Step-level reaction_type guidance + template hints for LLM repair
            steps = audit_details.get("steps", [])
            if isinstance(steps, list) and steps:
                repair_guidance_lines.append(f"\n### Route {route_id} Step-level Reaction-Type Guidance")
                for step in steps:
                    if not isinstance(step, dict):
                        continue
                    step_idx = step.get("step_index", "?")
                    rtype = step.get("reaction_type", {})
                    rtype_id = rtype.get("id", "unknown") if isinstance(rtype, dict) else "unknown"
                    rtype_conf = 0.0
                    rtype_source = "unknown"
                    if isinstance(rtype, dict):
                        rtype_conf = float(rtype.get("confidence", 0.0) or 0.0)
                        rtype_source = str(rtype.get("source", "unknown"))
                    repair_guidance_lines.append(
                        f"- Step {step_idx}: reaction_type={rtype_id} (confidence={rtype_conf:.2f}, source={rtype_source})"
                    )

                    if isinstance(rtype, dict):
                        byp = rtype.get("byproduct_candidates", [])
                        if isinstance(byp, list) and byp:
                            byp_text = []
                            for c in byp[:3]:
                                if not isinstance(c, dict):
                                    continue
                                smi = c.get("smiles")
                                conf = float(c.get("confidence", 0.0) or 0.0)
                                if isinstance(smi, str) and smi:
                                    byp_text.append(f"{smi}({conf:.2f})")
                            if byp_text:
                                repair_guidance_lines.append(
                                    f"  - suggested explicit byproducts: {', '.join(byp_text)}"
                                )

                    mapping = step.get("mapping", {})
                    if isinstance(mapping, dict):
                        completeness = float(mapping.get("completeness", 0.0) or 0.0)
                        notes = mapping.get("notes", [])
                        if completeness < 0.99:
                            repair_guidance_lines.append(
                                f"  - mapping completeness={completeness:.2f}; output fully mapped rs>>ps for this step."
                            )
                        if isinstance(notes, list):
                            for note in notes[:2]:
                                if isinstance(note, str) and note.strip():
                                    repair_guidance_lines.append(f"  - mapping note: {note}")

                    ab = step.get("atom_balance", {})
                    if isinstance(ab, dict):
                        miss = ab.get("missing_product_delta", {})
                        excess = ab.get("excess_product_delta", {})
                        if miss:
                            repair_guidance_lines.append(f"  - missing_product_delta: {miss}")
                        if excess:
                            repair_guidance_lines.append(f"  - excess_product_delta: {excess}")

                    step_template_hints = _repair_template_hints_for_step(step)
                    if step_template_hints:
                        repair_guidance_lines.append("  - reaction template suggestions:")
                        repair_guidance_lines.extend(step_template_hints)

                    # Optional local one-step persona guidance for high-risk steps.
                    hard = step.get("hard_issues", [])
                    if isinstance(hard, list) and hard:
                        try:
                            prompts_mod = _get_prompts_personas_mod()
                            local_prompt = prompts_mod.get_local_chemistry_repair_prompt(
                                target=str(repair_request.get("target_smiles", "")),
                                problematic_step=step,
                                context={
                                    "route_id": route_id,
                                    "hard_issues": hard,
                                    "soft_issues": step.get("soft_issues", []),
                                },
                            )
                            local_repair_lines.append(f"\n### Local Repair Persona (Route {route_id} Step {step_idx})")
                            local_repair_lines.append(local_prompt)
                        except Exception:
                            pass
    except Exception:
        # Fall back to basic prompt if module unavailable
        pass

    try:
        for route in rejected_routes:
            if not isinstance(route, dict):
                continue
            route_id = route.get("route_id", "?")
            audit_details = route.get("audit_details", {})
            steps = audit_details.get("steps", [])
            if not isinstance(steps, list):
                continue
            for step in steps:
                if not isinstance(step, dict):
                    continue
                step_idx = step.get("step_index", "?")
                hints = _repair_template_hints_for_step(step)
                if hints:
                    template_hint_lines.append(f"- Route {route_id} Step {step_idx}:")
                    template_hint_lines.extend(hints)
    except Exception:
        pass

    guidance_section = "\n".join(repair_guidance_lines) if repair_guidance_lines else ""
    local_section = "\n".join(local_repair_lines[:4]) if local_repair_lines else ""
    template_section = (
        "Reaction-template matching hints (from reaction_classifier rules):\n"
        + "\n".join(template_hint_lines[:24])
    ) if template_hint_lines else ""

    persona_section = ""
    try:
        prompts_mod = _get_prompts_personas_mod()
        persona_section = prompts_mod.get_repair_specialist_prompt(
            target=str(repair_request.get("target_smiles", "")),
            failed_routes=rejected_routes,
            audit_findings={
                "audit_summary": repair_request.get("audit_summary", {}),
                "round": repair_request.get("round"),
                "instruction": repair_request.get("instruction"),
            },
            output_mode="planner_routes_array",
            output_template=PLANNER_ROUTES_TEMPLATE,
        )
    except Exception:
        persona_section = ""

    return (
        "You are the host LLM retrosynthesis repair planner in a strict audit-gated loop.\n"
        "Given failed candidate routes and audit findings, generate corrected routes directly.\n"
        "Priority rules (highest first):\n"
        "1) Final output MUST be a planner routes JSON array only.\n"
        "2) Keep route_id stable where possible.\n"
        "3) Repair local chemistry in rs>>ps: fix invalid bond disconnections, reaction class mismatch, and atom balance.\n"
        "4) Prefer reaction templates suggested by step-level hints when compatible.\n\n"
        "Requirements:\n"
        "- Return strict JSON array only (no markdown, no explanation, no planner option exploration).\n"
        "- Preserve route_id values where possible.\n"
        "- Ensure each step reaction is chemically valid and atom-conserving (including hydrogen counts).\n"
        "- Include explicit byproducts when required for atom balance (e.g., H2O/HCl/alcohol byproducts).\n"
        "- Ensure the final step product matches target_smiles.\n"
        "- Fill audit_verdict/critical_issues from your own reasoning (do not leave placeholder text).\n"
        f"{template_section}\n\n"
        f"{guidance_section}\n\n"
        f"{local_section}\n\n"
        "Repair persona guidance (content guidance only; output contract above takes precedence):\n"
        f"{persona_section}\n\n"
        f"Repair request:\n{json.dumps(repair_request, indent=2, ensure_ascii=False)}\n\n"
        f"Target output template:\n{json.dumps(PLANNER_ROUTES_TEMPLATE, indent=2, ensure_ascii=False)}\n"
    )


def _resolve_planner_routes_candidate(args: argparse.Namespace, output_dir: Path) -> Path:
    return resolve_planner_routes_candidate(args, output_dir)


def _render_host_command(command_template: str, values: Dict[str, Any]) -> str:
    return render_host_command(command_template, values)


def _invoke_host_planner(command_template: str, values: Dict[str, Any], log_path: Path) -> Dict[str, Any]:
    return invoke_host_planner(command_template, values, log_path, _write_json)


def _planner_request_payload(
    main_start: float,
    planner_request_path: Path,
    planner_routes_template_path: Path,
    planner_prompt_path: Path,
    planner_routes_candidate_path: Path,
    schema_dir: Path,
    mode: str,
    next_action: str,
    extra_outputs: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    return planner_request_payload(
        main_start=main_start,
        planner_request_path=planner_request_path,
        planner_routes_template_path=planner_routes_template_path,
        planner_prompt_path=planner_prompt_path,
        planner_routes_candidate_path=planner_routes_candidate_path,
        schema_dir=schema_dir,
        mode=mode,
        next_action=next_action,
        extra_outputs=extra_outputs,
    )


def _run_host_planner_repair_loop(
    args: argparse.Namespace,
    main_start: float,
    output_dir: Path,
    schema_dir: Path,
    planner_request: Dict[str, Any],
    planner_request_path: Path,
    planner_routes_template_path: Path,
    planner_prompt_path: Path,
    planner_routes_candidate_path: Path,
    strategy_for_planner: Dict[str, Any],
    require_all_routes_recommendable: bool = True,
    min_recommended_routes: int = 1,
) -> Tuple[Optional[Dict[str, Any]], Optional[Path], Dict[str, Any]]:
    return run_host_planner_repair_loop(
        args=args,
        main_start=main_start,
        output_dir=output_dir,
        schema_dir=schema_dir,
        planner_request=planner_request,
        planner_request_path=planner_request_path,
        planner_routes_template_path=planner_routes_template_path,
        planner_prompt_path=planner_prompt_path,
        planner_routes_candidate_path=planner_routes_candidate_path,
        strategy_for_planner=strategy_for_planner,
        write_json=_write_json,
        read_json=_read_json,
        validate_routes=_validate_routes,
        audit_routes_with_gate=audit_routes_with_gate,
        load_module=_load_module,
        build_repair_prompt_fn=_build_repair_prompt,
        require_all_routes_recommendable=require_all_routes_recommendable,
        min_recommended_routes=min_recommended_routes,
    )


def _ensure_required_args(args: argparse.Namespace) -> None:
    if not args.target_smiles:
        raise ValueError("--target-smiles is required unless --emit-vis-plan-template-only is used")
    if not args.routes_file and not args.auto_propose_routes:
        raise ValueError("--routes-file is required unless --emit-vis-plan-template-only is used")


def _enforce_min_route_count(routes: Any, min_routes: int = MIN_REQUIRED_ROUTES) -> None:
    if not isinstance(routes, list):
        raise ValueError("routes must be a JSON array")
    if len(routes) < min_routes:
        raise ValueError(
            f"At least {min_routes} routes are required, got {len(routes)}. "
            "Provide additional routes or use planner repair loop."
        )


def _allow_partial_routes_after_repair_limit(args: argparse.Namespace) -> bool:
    loop = getattr(args, "planner_loop_summary", None)
    if not isinstance(loop, dict):
        return False
    return str(loop.get("termination_reason", "")).strip().lower() == "max_iterations_reached"


def _generate_all_persona_prompts(
    prompts_mod,
    output_dir: Path,
    target_smiles: str,
    strategy: Dict[str, Any],
    audited_routes: List[Dict[str, Any]],
    functional_groups: Dict[str, Any],
    features: Dict[str, Any],
) -> Dict[str, str]:
    """
    Generate prompts for all 6 personas and save them to the output directory.
    
    Personas:
    1. Top-Level Designer - global strategy design
    2. Reaction Designer - reaction design details
    3. Auditor - audit and validation
    4. Integration Specialist - integrated recommendation report
    5. Visualization Specialist - visualization planning
    6. Organic Synthesis Repair Specialist - local chemistry repair
    
    Returns:
        Dict[str, str]: Mapping of persona name to prompt file path
    """
    prompts_dir = output_dir / "persona_prompts"
    prompts_dir.mkdir(parents=True, exist_ok=True)
    
    prompt_paths: Dict[str, str] = {}
    
    # Build context information
    broad_context = f"""
Target SMILES: {target_smiles}
Functional Groups: {functional_groups.get('functional_groups', {})}
Protecting Groups: {functional_groups.get('protecting_groups', {}).get('counts', {})}
Synthesis Critical Properties: {features.get('synthesis_critical', {})}
"""
    
    # 1. Top-Level Designer (Global Strategist)
    try:
        top_level_prompt = prompts_mod.get_top_level_designer_prompt(
            target=target_smiles,
            broad_context=broad_context
        )
        path = prompts_dir / "1_top_level_designer.txt"
        path.write_text(top_level_prompt, encoding="utf-8")
        prompt_paths["top_level_designer"] = str(path)
    except Exception as e:
        prompt_paths["top_level_designer_error"] = str(e)
    
    # 2. Reaction Designer (Tactician) - always generate (use fallback if no directions)
    try:
        analysis = strategy.get("analysis", {}) if isinstance(strategy, dict) else {}
        directions = strategy.get("directions", []) if isinstance(strategy, dict) else []
        # Use the first direction; if absent, build a fallback direction to avoid missing persona2.
        if directions and len(directions) > 0 and isinstance(directions[0], dict):
            direction = directions[0]
        else:
            key_disconnections = analysis.get("key_disconnections", []) if isinstance(analysis, dict) else []
            disconnection = key_disconnections[0] if isinstance(key_disconnections, list) and key_disconnections else "Strategic disconnection from audited route"
            direction = {
                "direction_id": 1,
                "disconnection_type": "Route-anchored Disconnection",
                "key_bond": "context-dependent",
                "reasoning": str(disconnection),
            }

        precursors = []
        if audited_routes and isinstance(audited_routes[0], dict):
            precursors = audited_routes[0].get("precursors", [])[:3]

        reaction_prompt = prompts_mod.get_reaction_designer_prompt(
            target=target_smiles,
            direction_info=direction,
            precursors=precursors
        )
        path = prompts_dir / "2_reaction_designer.txt"
        path.write_text(reaction_prompt, encoding="utf-8")
        prompt_paths["reaction_designer"] = str(path)
    except Exception as e:
        prompt_paths["reaction_designer_error"] = str(e)
    
    # 3. Auditor (The Critic) - audit routes
    try:
        if audited_routes:
            first_route = audited_routes[0] if audited_routes else {}
            reaction_design = first_route.get("steps", [{}])[0] if first_route.get("steps") else {}
            auditor_prompt = prompts_mod.get_auditor_prompt(
                target=target_smiles,
                reaction_design=reaction_design,
                history=""
            )
            path = prompts_dir / "3_auditor.txt"
            path.write_text(auditor_prompt, encoding="utf-8")
            prompt_paths["auditor"] = str(path)
    except Exception as e:
        prompt_paths["auditor_error"] = str(e)
    
    # 4. Integration Specialist (The Synthesizer)
    try:
        integration_prompt = prompts_mod.get_integration_prompt(audited_routes)
        path = prompts_dir / "4_integration_specialist.txt"
        path.write_text(integration_prompt, encoding="utf-8")
        prompt_paths["integration_specialist"] = str(path)
    except Exception as e:
        prompt_paths["integration_specialist_error"] = str(e)
    
    # 5. Visualization Specialist (Creative Director)
    try:
        vis_prompt = prompts_mod.get_visualization_specialist_prompt(
            target_smiles, strategy, audited_routes
        )
        path = prompts_dir / "5_visualization_specialist.txt"
        path.write_text(vis_prompt, encoding="utf-8")
        prompt_paths["visualization_specialist"] = str(path)
    except Exception as e:
        prompt_paths["visualization_specialist_error"] = str(e)

    # 6. Organic Synthesis Repair Specialist (The Chemist)
    try:
        failed_routes = []
        for route in audited_routes:
            if not isinstance(route, dict):
                continue
            verdict = str(route.get("audit_verdict", "")).upper()
            if verdict in {"FAIL", "CONDITIONAL"}:
                failed_routes.append(route)
        repair_prompt = prompts_mod.get_repair_specialist_prompt(
            target=target_smiles,
            failed_routes=failed_routes,
            audit_findings={"note": "Generated from audited routes for repair-loop guidance."},
            output_mode="planner_routes_array",
            output_template=PLANNER_ROUTES_TEMPLATE,
        )
        path = prompts_dir / "6_repair_specialist.txt"
        path.write_text(repair_prompt, encoding="utf-8")
        prompt_paths["repair_specialist"] = str(path)
    except Exception as e:
        prompt_paths["repair_specialist_error"] = str(e)
    
    return prompt_paths


def run_pipeline(args: argparse.Namespace) -> Dict[str, Any]:
    pipeline_start = time.perf_counter()
    stage_timers: Dict[str, float] = {}

    skill_root = Path(__file__).resolve().parents[1]
    output_dir = Path(args.output_dir).resolve()
    images_dir = output_dir / "images"
    images_dir.mkdir(parents=True, exist_ok=True)

    retroskill_dir = skill_root / "internal" / "retroskill"
    rdkit_utils_dir = skill_root / "internal" / "rdkit_utils"
    atom_mappers_dir = skill_root / "internal" / "atom_mappers"

    t = time.perf_counter()
    prompts_mod = _load_module("prompts_personas", retroskill_dir / "prompts_personas.py")
    report_mod = _load_module("report_generator", retroskill_dir / "report_generator.py")
    vis_mod = _load_module("visualize_routes", retroskill_dir / "scripts" / "visualize_routes.py")
    analyze_mod = _load_module("analyze_structure", rdkit_utils_dir / "analyze_structure.py")
    features_mod = _load_module("calculate_features", rdkit_utils_dir / "calculate_features.py")
    conformer_mod = _load_module("generate_conformer", rdkit_utils_dir / "generate_conformer.py")
    fg_mod = _load_module("functional_groups", rdkit_utils_dir / "functional_groups.py")
    rxn_mapper_mod = _load_module("rxn_mapper", atom_mappers_dir / "rxn_mapper.py")
    atom_utils_mod = _load_module("atom_utils", atom_mappers_dir / "utils.py")
    stage_timers["module_load_seconds"] = round(time.perf_counter() - t, 6)

    t = time.perf_counter()
    strategy = _read_json(args.strategy_file, {"analysis": {}})
    routes = _coerce_routes_payload(_read_json(args.routes_file, []))
    vis_plan = _read_json(args.vis_plan_file, {})
    stage_timers["input_read_seconds"] = round(time.perf_counter() - t, 6)

    t = time.perf_counter()
    route_ids = [r.get("route_id") for r in routes if isinstance(r, dict) and isinstance(r.get("route_id"), int)]
    validation_errors: List[str] = []
    validation_errors.extend(_validate_strategy(strategy))
    validation_errors.extend(_validate_routes(routes))
    validation_errors.extend(_validate_vis_plan(vis_plan, route_ids))
    if validation_errors:
        raise ValueError("Schema validation failed: " + " | ".join(validation_errors))
    stage_timers["validation_seconds"] = round(time.perf_counter() - t, 6)

    t = time.perf_counter()
    schema_dir = output_dir / "schemas"
    _write_schema_files(schema_dir)

    vis_template_path = Path(args.emit_vis_plan_template) if args.emit_vis_plan_template else (output_dir / "vis_plan.template.json")
    _write_json(vis_template_path, VIS_PLAN_TEMPLATE)
    stage_timers["schema_template_emit_seconds"] = round(time.perf_counter() - t, 6)

    t = time.perf_counter()
    structure = analyze_mod.analyze_smiles(args.target_smiles, generate_3d=False)
    features = features_mod.calculate_features(args.target_smiles)
    functional_groups = fg_mod.analyze_functional_groups(args.target_smiles)
    conformer = conformer_mod.generate_conformer(args.target_smiles, force_field=args.force_field, constrained_smarts=None)
    stage_timers["rdkit_analysis_seconds"] = round(time.perf_counter() - t, 6)

    t = time.perf_counter()
    audited_routes, recommendable_routes, rejected_routes, audit_summary = audit_routes_with_gate(
        routes=routes,
        target_smiles=args.target_smiles,
        strict_audit=args.strict_audit,
        rxn_mapper_mod=rxn_mapper_mod,
        atom_utils_mod=atom_utils_mod,
        strategy=strategy if isinstance(strategy, dict) else {},
        analyze_mod=analyze_mod,
        features_mod=features_mod,
        fg_mod=fg_mod,
        enforce_style_mix=bool(getattr(args, "enforce_style_mix", False)),
    )
    if (
        args.strict_audit
        and args.auto_propose_routes
        and rejected_routes
        and not isinstance(getattr(args, "planner_loop_summary", None), dict)
    ):
        raise ValueError(
            "Strict audit found rejected routes but no repair-loop trace was recorded. "
            "Rerun with --planner-auto-repair (default) and provide --host-planner-command if automatic repair is needed."
        )
    
    # Generate prompts for all 6 personas
    persona_prompt_paths = _generate_all_persona_prompts(
        prompts_mod=prompts_mod,
        output_dir=output_dir,
        target_smiles=args.target_smiles,
        strategy=strategy if isinstance(strategy, dict) else {},
        audited_routes=audited_routes,
        functional_groups=functional_groups,
        features=features,
    )
    
    # Backward compatibility: also save a standalone visualization_prompt.txt
    vis_prompt = prompts_mod.get_visualization_specialist_prompt(args.target_smiles, strategy, audited_routes)
    vis_prompt_path = output_dir / "visualization_prompt.txt"
    vis_prompt_path.write_text(vis_prompt, encoding="utf-8")
    stage_timers["prompt_write_seconds"] = round(time.perf_counter() - t, 6)

    selected_routes = _pick_routes(recommendable_routes, vis_plan, args.top_k)
    partial_after_limit = _allow_partial_routes_after_repair_limit(args)
    require_min_recommendable = bool(
        args.strict_audit
        and (
            bool(getattr(args, "auto_propose_routes", False))
            or bool(isinstance(getattr(args, "planner_loop_summary", None), dict))
        )
        and (not partial_after_limit)
    )
    if require_min_recommendable and len(recommendable_routes) < MIN_REQUIRED_ROUTES:
        raise ValueError(
            f"Strict audit requires at least {MIN_REQUIRED_ROUTES} recommendable routes for final report, "
            f"got {len(recommendable_routes)}."
        )
    render_result = render_and_report(
        args=args,
        output_dir=output_dir,
        images_dir=images_dir,
        vis_plan=vis_plan,
        vis_mod=vis_mod,
        report_mod=report_mod,
        audited_routes=audited_routes,
        selected_routes=selected_routes,
        rejected_routes=rejected_routes,
        strategy=strategy,
        audit_summary=audit_summary,
    )
    image_paths = render_result["image_paths"]
    report_path = render_result["report_path"]
    stage_timers["render_seconds"] = render_result["render_seconds"]
    stage_timers["report_seconds"] = render_result["report_seconds"]

    summary = {
        "success": True,
        "target_smiles": args.target_smiles,
        "selected_route_ids": [r.get("route_id") for r in selected_routes],
        "visualized_route_ids": [r.get("route_id") for r in (audited_routes if args.strict_audit else selected_routes)],
        "rejected_route_ids": [r.get("route_id") for r in rejected_routes],
        "audit": audit_summary,
        "outputs": {
            "report": str(report_path),
            "visualization_prompt": str(vis_prompt_path),
            "vis_plan_template": str(vis_template_path),
            "schemas": {
                "routes": str(schema_dir / "routes.schema.json"),
                "strategy": str(schema_dir / "strategy.schema.json"),
                "vis_plan": str(schema_dir / "vis_plan.schema.json"),
                "planner_request": str(schema_dir / "planner_request.schema.json"),
                "planner_routes": str(schema_dir / "planner_routes.schema.json")
            },
            "image_paths": image_paths
        },
        "rdkit": {
            "structure": structure,
            "features": features,
            "conformer": {
                "success": conformer.get("success", False),
                "force_field": conformer.get("conformer_data", {}).get("force_field"),
                "final_energy": conformer.get("conformer_data", {}).get("final_energy")
            },
            "functional_groups": functional_groups,
        },
        "timings": {
            "pipeline_seconds": round(time.perf_counter() - pipeline_start, 6),
            "stages": stage_timers,
        },
    }
    planner_loop = getattr(args, "planner_loop_summary", None)
    loop_enabled = bool(isinstance(planner_loop, dict) and planner_loop.get("rounds"))
    if isinstance(planner_loop, dict):
        if "final_round" not in planner_loop:
            rounds = planner_loop.get("rounds", [])
            if isinstance(rounds, list) and rounds:
                last = rounds[-1]
                if isinstance(last, dict) and isinstance(last.get("round"), int):
                    planner_loop["final_round"] = last.get("round")
        summary["planner_loop"] = planner_loop
    else:
        summary["planner_loop"] = {"enabled": loop_enabled, "final_round": 0 if loop_enabled else None, "rounds": []}
    if isinstance(getattr(args, "route_count_policy", None), dict):
        summary["route_count_policy"] = args.route_count_policy
    return summary


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run RDKit + retrosynthesis visualization/report closed loop.")
    parser.add_argument("--target-smiles", default="", help="Target molecule SMILES.")
    parser.add_argument("--routes-file", default="", help="Path to audited routes JSON.")
    parser.add_argument("--strategy-file", default="", help="Path to strategy JSON.")
    parser.add_argument("--vis-plan-file", default="", help="Path to visualization plan JSON.")
    parser.add_argument("--output-dir", default="closed_loop_output", help="Output directory.")
    parser.add_argument("--top-k", type=int, default=1, help="Fallback selected route count.")
    parser.add_argument("--target-legend", default="Target Molecule", help="Default legend when vis plan does not provide one.")
    parser.add_argument("--force-field", choices=["MMFF", "UFF"], default="MMFF", help="Force field for conformer generation.")
    parser.add_argument("--emit-vis-plan-template", default="", help="Write vis plan template JSON to this path (default: <output-dir>/vis_plan.template.json).")
    parser.add_argument("--emit-vis-plan-template-only", action="store_true", help="Only write vis plan template and schema files, then exit.")
    parser.add_argument("--validate-only", action="store_true", help="Validate input JSON and exit without rendering/report generation.")
    parser.add_argument("--auto-propose-routes", action="store_true", help="Enable auto route planning interface when routes file is absent.")
    parser.add_argument("--planner-backend", choices=["host"], default="host", help="Planner backend used by --auto-propose-routes.")
    parser.add_argument("--auto-propose-max-routes", type=int, default=5, help="Requested routes in auto planning (bounded to 3~5).")
    parser.add_argument("--planner-routes-file", default="", help="Optional planner returned routes file. Default: <output-dir>/planner_routes.json")
    parser.add_argument("--planner-request-only", action="store_true", help="Emit planner request/template outputs only, then exit.")
    parser.add_argument(
        "--planner-auto-repair",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Enable iterative planner->audit->repair loop for host backend (enabled by default; use --no-planner-auto-repair to disable).",
    )
    parser.add_argument("--planner-max-iterations", type=int, default=5, help="Max iterations for planner auto-repair loop (bounded to 5).")
    parser.add_argument("--host-planner-command", default="", help="Optional host command template to run LLM planner automatically. Placeholders: {prompt_file} {output_file} {request_file} {output_dir} {stage} {round}.")
    parser.add_argument("--strict-audit", action="store_true", help="Recompute chemistry audit verdict and gate FAIL routes from recommendation/report.")
    return parser.parse_args()


def main() -> int:
    main_start = time.perf_counter()
    args = parse_args()
    # Hard requirement: all operational runs use strict audit mode.
    args.strict_audit = True
    args.enforce_style_mix = bool(args.auto_propose_routes)
    out_path = Path(args.output_dir) / "run_summary.json"
    try:
        output_dir = Path(args.output_dir).resolve()
        schema_dir = output_dir / "schemas"
        _write_schema_files(schema_dir)

        vis_template_path = Path(args.emit_vis_plan_template) if args.emit_vis_plan_template else (output_dir / "vis_plan.template.json")
        _write_json(vis_template_path, VIS_PLAN_TEMPLATE)

        if args.emit_vis_plan_template_only:
            payload = {
                "success": True,
                "mode": "template_only",
                "timings": {"total_seconds": round(time.perf_counter() - main_start, 6)},
                "outputs": {
                    "vis_plan_template": str(vis_template_path),
                    "schemas": {
                        "routes": str(schema_dir / "routes.schema.json"),
                        "strategy": str(schema_dir / "strategy.schema.json"),
                        "vis_plan": str(schema_dir / "vis_plan.schema.json"),
                        "planner_request": str(schema_dir / "planner_request.schema.json"),
                        "planner_routes": str(schema_dir / "planner_routes.schema.json")
                    }
                }
            }
            _write_json(out_path, payload)
            print(json.dumps(payload, indent=2, ensure_ascii=False))
            return 0

        strategy_for_planner = _read_json(args.strategy_file, {"analysis": {}})
        target_chemistry_for_planner: Dict[str, Any] = {}
        if args.target_smiles:
            skill_root = Path(__file__).resolve().parents[1]
            rdkit_utils_dir = skill_root / "internal" / "rdkit_utils"
            analyze_plan_mod = _load_module("analyze_structure_planner", rdkit_utils_dir / "analyze_structure.py")
            features_plan_mod = _load_module("calculate_features_planner", rdkit_utils_dir / "calculate_features.py")
            fg_plan_mod = _load_module("functional_groups_planner", rdkit_utils_dir / "functional_groups.py")
            target_chemistry_for_planner = _build_target_chemistry_context(
                args.target_smiles,
                analyze_plan_mod,
                features_plan_mod,
                fg_plan_mod,
            )
            if not bool(target_chemistry_for_planner.get("success")):
                raise ValueError(
                    "Pre-planning chemistry analysis failed; cannot continue planner generation. "
                    f"target_chemistry error: {target_chemistry_for_planner.get('error')}"
                )
        requested_planner_routes = max(1, min(args.auto_propose_max_routes, 10))
        if args.auto_propose_routes:
            requested_planner_routes = max(
                MIN_REQUIRED_ROUTES,
                requested_planner_routes,
            )
        planner_request = build_planner_request(
            target_smiles=args.target_smiles,
            max_routes=requested_planner_routes,
            strategy=strategy_for_planner,
            target_chemistry=target_chemistry_for_planner,
        )
        planner_request_path = output_dir / "planner_request.json"
        planner_routes_template_path = output_dir / "planner_routes.template.json"
        planner_prompt_path = output_dir / "planner_prompt.txt"
        planner_routes_candidate_path = _resolve_planner_routes_candidate(args, output_dir)

        if args.planner_request_only:
            if not args.target_smiles:
                raise ValueError("--target-smiles is required for --planner-request-only")
            _write_json(planner_request_path, planner_request)
            _write_json(planner_routes_template_path, PLANNER_ROUTES_TEMPLATE)
            planner_prompt_path.write_text(build_planner_prompt(planner_request), encoding="utf-8")
            payload = {
                "success": True,
                "mode": "planner_request_only",
                "timings": {"total_seconds": round(time.perf_counter() - main_start, 6)},
                "outputs": {
                    "planner_request": str(planner_request_path),
                    "planner_routes_template": str(planner_routes_template_path),
                    "planner_prompt": str(planner_prompt_path),
                    "schemas": {
                        "planner_request": str(schema_dir / "planner_request.schema.json"),
                        "planner_routes": str(schema_dir / "planner_routes.schema.json")
                    },
                    "planner_routes_expected": str(planner_routes_candidate_path)
                }
            }
            _write_json(out_path, payload)
            print(json.dumps(payload, indent=2, ensure_ascii=False))
            return 0

        _ensure_required_args(args)
        if args.auto_propose_routes and args.planner_auto_repair:
            host_cmd_text = str(args.host_planner_command or "").strip()
            if host_cmd_text and _is_non_production_host_command(host_cmd_text):
                raise ValueError(
                    "planner auto-repair requires a real host planner command "
                    "(placeholder/fallback command is not allowed)."
                )

        if args.auto_propose_routes and not args.routes_file:
            if args.planner_auto_repair or args.host_planner_command:
                payload, resolved_routes_path, planner_loop_summary = _run_host_planner_repair_loop(
                    args=args,
                    main_start=main_start,
                    output_dir=output_dir,
                    schema_dir=schema_dir,
                    planner_request=planner_request,
                    planner_request_path=planner_request_path,
                    planner_routes_template_path=planner_routes_template_path,
                    planner_prompt_path=planner_prompt_path,
                    planner_routes_candidate_path=planner_routes_candidate_path,
                    strategy_for_planner=strategy_for_planner if isinstance(strategy_for_planner, dict) else {},
                    require_all_routes_recommendable=False,
                    min_recommended_routes=MIN_REQUIRED_ROUTES,
                )
                if payload is not None:
                    _write_json(out_path, payload)
                    print(json.dumps(payload, indent=2, ensure_ascii=False))
                    return 0
                if resolved_routes_path is None:
                    raise ValueError("planner auto-repair loop did not return a usable routes file")
                args.routes_file = str(resolved_routes_path)
                args.strict_audit = True
                args.planner_loop_summary = planner_loop_summary
            elif planner_routes_candidate_path.exists():
                args.routes_file = str(planner_routes_candidate_path)
            else:
                _write_json(planner_request_path, planner_request)
                _write_json(planner_routes_template_path, PLANNER_ROUTES_TEMPLATE)
                planner_prompt_path.write_text(build_planner_prompt(planner_request), encoding="utf-8")
                payload = {
                    "success": True,
                    "mode": "planner_request_only",
                    "timings": {"total_seconds": round(time.perf_counter() - main_start, 6)},
                    "outputs": {
                        "planner_request": str(planner_request_path),
                        "planner_routes_template": str(planner_routes_template_path),
                        "planner_prompt": str(planner_prompt_path),
                        "schemas": {
                            "planner_request": str(schema_dir / "planner_request.schema.json"),
                            "planner_routes": str(schema_dir / "planner_routes.schema.json")
                        },
                        "planner_routes_expected": str(planner_routes_candidate_path)
                    },
                    "next_action": "Have host LLM generate and execute the full host_planner_command, write planner routes to planner_routes_expected, then rerun the same command."
                }
                _write_json(out_path, payload)
                print(json.dumps(payload, indent=2, ensure_ascii=False))
                return 0

        if (
            args.strict_audit
            and args.planner_auto_repair
            and (not args.validate_only)
            and args.routes_file
            and not isinstance(getattr(args, "planner_loop_summary", None), dict)
            and (args.auto_propose_routes or bool(str(args.host_planner_command or "").strip()))
        ):
            payload, resolved_routes_path, planner_loop_summary = _run_host_planner_repair_loop(
                args=args,
                main_start=main_start,
                output_dir=output_dir,
                schema_dir=schema_dir,
                planner_request=planner_request,
                planner_request_path=planner_request_path,
                planner_routes_template_path=planner_routes_template_path,
                planner_prompt_path=planner_prompt_path,
                planner_routes_candidate_path=Path(args.routes_file).resolve(),
                strategy_for_planner=strategy_for_planner if isinstance(strategy_for_planner, dict) else {},
                require_all_routes_recommendable=False,
                min_recommended_routes=MIN_REQUIRED_ROUTES,
            )
            if payload is not None:
                _write_json(out_path, payload)
                print(json.dumps(payload, indent=2, ensure_ascii=False))
                return 0
            if resolved_routes_path is None:
                raise ValueError("strict repair loop did not return a usable routes file")
            args.routes_file = str(resolved_routes_path)
            args.strict_audit = True
            args.planner_loop_summary = planner_loop_summary

        strategy = _read_json(args.strategy_file, {"analysis": {}})
        routes = _coerce_routes_payload(_read_json(args.routes_file, []))
        vis_plan = _read_json(args.vis_plan_file, {})
        partial_routes_allowed = _allow_partial_routes_after_repair_limit(args)
        if partial_routes_allowed:
            if not isinstance(routes, list):
                raise ValueError("routes must be a JSON array")
            args.route_count_policy = {
                "minimum_required": MIN_REQUIRED_ROUTES,
                "provided": len(routes),
                "partial_allowed_due_to": "max_iterations_reached",
            }
        elif args.auto_propose_routes and not args.validate_only:
            _enforce_min_route_count(routes, min_routes=MIN_REQUIRED_ROUTES)
        route_ids = [r.get("route_id") for r in routes if isinstance(r, dict) and isinstance(r.get("route_id"), int)]

        validation_errors: List[str] = []
        validation_errors.extend(_validate_strategy(strategy))
        validation_errors.extend(_validate_routes(routes))
        validation_errors.extend(_validate_vis_plan(vis_plan, route_ids))
        if validation_errors:
            raise ValueError("Schema validation failed: " + " | ".join(validation_errors))

        if args.validate_only:
            validate_payload: Dict[str, Any] = {}
            if args.strict_audit:
                skill_root = Path(__file__).resolve().parents[1]
                rxn_mapper_mod = _load_module("rxn_mapper", skill_root / "internal" / "atom_mappers" / "rxn_mapper.py")
                atom_utils_mod = _load_module("atom_utils", skill_root / "internal" / "atom_mappers" / "utils.py")
                analyze_mod = _load_module("analyze_structure_validate", skill_root / "internal" / "rdkit_utils" / "analyze_structure.py")
                features_mod = _load_module("calculate_features_validate", skill_root / "internal" / "rdkit_utils" / "calculate_features.py")
                fg_mod = _load_module("functional_groups_validate", skill_root / "internal" / "rdkit_utils" / "functional_groups.py")
                _, _, _, validate_audit_summary = audit_routes_with_gate(
                    routes=routes,
                    target_smiles=args.target_smiles,
                    strict_audit=True,
                    rxn_mapper_mod=rxn_mapper_mod,
                    atom_utils_mod=atom_utils_mod,
                    strategy=strategy if isinstance(strategy, dict) else {},
                    analyze_mod=analyze_mod,
                    features_mod=features_mod,
                    fg_mod=fg_mod,
                    enforce_style_mix=bool(getattr(args, "enforce_style_mix", False)),
                )
                validate_payload["audit"] = validate_audit_summary
            payload = {
                "success": True,
                "mode": "validate_only",
                "timings": {"total_seconds": round(time.perf_counter() - main_start, 6)},
                "validated": {
                    "strategy_file": args.strategy_file,
                    "routes_file": args.routes_file,
                    "vis_plan_file": args.vis_plan_file
                },
                "outputs": {
                    "vis_plan_template": str(vis_template_path),
                    "schemas": {
                        "routes": str(schema_dir / "routes.schema.json"),
                        "strategy": str(schema_dir / "strategy.schema.json"),
                        "vis_plan": str(schema_dir / "vis_plan.schema.json"),
                        "planner_request": str(schema_dir / "planner_request.schema.json"),
                        "planner_routes": str(schema_dir / "planner_routes.schema.json")
                    }
                },
                **validate_payload,
            }
            _write_json(out_path, payload)
            print(json.dumps(payload, indent=2, ensure_ascii=False))
            return 0

        summary = run_pipeline(args)
        summary.setdefault("timings", {})
        summary["timings"]["total_seconds"] = round(time.perf_counter() - main_start, 6)
        _write_json(out_path, summary)
        print(json.dumps(summary, indent=2, ensure_ascii=False))
        return 0
    except Exception as exc:
        payload = {
            "success": False,
            "error": str(exc),
            "timings": {"total_seconds": round(time.perf_counter() - main_start, 6)},
        }
        _write_json(out_path, payload)
        print(json.dumps(payload, indent=2, ensure_ascii=False))
        return 1


if __name__ == "__main__":
    sys.exit(main())
