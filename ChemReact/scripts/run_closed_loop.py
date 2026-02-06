import argparse
import json
import sys
from pathlib import Path
import importlib.util
from typing import Any, Dict, List


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


def _load_module(module_name: str, file_path: Path):
    spec = importlib.util.spec_from_file_location(module_name, str(file_path))
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot load module from {file_path}")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


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
        if "critical_issues" in route and not isinstance(route["critical_issues"], list):
            errs.append(f"routes[{i}].critical_issues must be an array")
        if "precursors" in route and not isinstance(route["precursors"], list):
            errs.append(f"routes[{i}].precursors must be an array")
        if "steps" in route:
            errs.extend(_validate_steps(route["steps"], i))
    return errs


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

    ranked = sorted(all_routes, key=lambda x: x.get("score", 0), reverse=True)
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
                    "strategy_type": {"type": "string"}
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
    _write_json(schemas_dir / "routes.schema.json", routes_schema)
    _write_json(schemas_dir / "strategy.schema.json", strategy_schema)
    _write_json(schemas_dir / "vis_plan.schema.json", vis_plan_schema)


def _ensure_required_args(args: argparse.Namespace) -> None:
    if not args.target_smiles:
        raise ValueError("--target-smiles is required unless --emit-vis-plan-template-only is used")
    if not args.routes_file:
        raise ValueError("--routes-file is required unless --emit-vis-plan-template-only is used")


def run_pipeline(args: argparse.Namespace) -> Dict[str, Any]:
    skill_root = Path(__file__).resolve().parents[1]
    output_dir = Path(args.output_dir).resolve()
    images_dir = output_dir / "images"
    images_dir.mkdir(parents=True, exist_ok=True)

    retroskill_dir = skill_root / "internal" / "retroskill"
    rdkit_utils_dir = skill_root / "internal" / "rdkit_utils"

    prompts_mod = _load_module("prompts_personas", retroskill_dir / "prompts_personas.py")
    report_mod = _load_module("report_generator", retroskill_dir / "report_generator.py")
    vis_mod = _load_module("visualize_routes", retroskill_dir / "scripts" / "visualize_routes.py")
    analyze_mod = _load_module("analyze_structure", rdkit_utils_dir / "analyze_structure.py")
    features_mod = _load_module("calculate_features", rdkit_utils_dir / "calculate_features.py")
    conformer_mod = _load_module("generate_conformer", rdkit_utils_dir / "generate_conformer.py")

    strategy = _read_json(args.strategy_file, {"analysis": {}})
    routes = _read_json(args.routes_file, [])
    vis_plan = _read_json(args.vis_plan_file, {})

    route_ids = [r.get("route_id") for r in routes if isinstance(r, dict) and isinstance(r.get("route_id"), int)]
    validation_errors: List[str] = []
    validation_errors.extend(_validate_strategy(strategy))
    validation_errors.extend(_validate_routes(routes))
    validation_errors.extend(_validate_vis_plan(vis_plan, route_ids))
    if validation_errors:
        raise ValueError("Schema validation failed: " + " | ".join(validation_errors))

    schema_dir = output_dir / "schemas"
    _write_schema_files(schema_dir)

    vis_template_path = Path(args.emit_vis_plan_template) if args.emit_vis_plan_template else (output_dir / "vis_plan.template.json")
    _write_json(vis_template_path, VIS_PLAN_TEMPLATE)

    structure = analyze_mod.analyze_smiles(args.target_smiles, generate_3d=False)
    features = features_mod.calculate_features(args.target_smiles)
    conformer = conformer_mod.generate_conformer(args.target_smiles, force_field=args.force_field, constrained_smarts=None)

    vis_prompt = prompts_mod.get_visualization_specialist_prompt(args.target_smiles, strategy, routes)
    vis_prompt_path = output_dir / "visualization_prompt.txt"
    vis_prompt_path.write_text(vis_prompt, encoding="utf-8")

    selected_routes = _pick_routes(routes, vis_plan, args.top_k)
    image_paths: Dict[str, str] = {}

    from rdkit import Chem

    target_mol = Chem.MolFromSmiles(args.target_smiles)
    target_atom_count = target_mol.GetNumAtoms() if target_mol else 0
    target_plan = vis_plan.get("target_image", {}) if isinstance(vis_plan, dict) else {}
    target_highlights = _safe_int_list(target_plan.get("highlight_atoms", []), target_atom_count)
    target_legend = target_plan.get("legend", args.target_legend)
    target_path = images_dir / "target.png"
    vis_mod.generate_molecule_image(args.target_smiles, str(target_path), legend=target_legend, highlight_atoms=target_highlights)
    image_paths["target_image"] = "images/target.png"

    for route in selected_routes:
        route_id = route.get("route_id")
        if route_id is None:
            continue

        route_plan = _get_route_plan(vis_plan, route_id) if isinstance(vis_plan, dict) else {}
        captions = _step_caption_map(route_plan)

        precursors = route.get("precursors", [])
        if isinstance(precursors, list) and precursors and route_plan.get("overview_grid", {}).get("enabled", True):
            legends = route_plan.get("overview_grid", {}).get("precursor_legends", [])
            out_name = f"route_{route_id}_overview.png"
            vis_mod.generate_route_grid(precursors, legends, str(images_dir / out_name))
            image_paths[f"route_{route_id}_overview"] = f"images/{out_name}"

        if isinstance(precursors, list) and precursors and route_plan.get("tree_view", {}).get("enabled", True):
            out_name = f"route_{route_id}_tree.png"
            vis_mod.generate_reaction_tree_image(args.target_smiles, precursors, str(images_dir / out_name))
            image_paths[f"route_{route_id}_tree"] = f"images/{out_name}"

        steps = route.get("steps", [])
        for step_i, step in enumerate(steps, start=1):
            rxn = step.get("reaction_smiles") or step.get("rxn_smiles") or step.get("smirks")
            if not rxn:
                continue
            out_name = f"route_{route_id}_step_{step_i}.png"
            vis_mod.generate_reaction_image(rxn, str(images_dir / out_name))
            image_paths[f"route_{route_id}_step_{step_i}"] = f"images/{out_name}"
            if step_i in captions:
                step["visual_caption"] = captions[step_i]

    report_path = output_dir / "RETRO_REPORT.md"
    report_mod.generate_report(args.target_smiles, strategy, selected_routes, image_paths, str(report_path))

    summary = {
        "success": True,
        "target_smiles": args.target_smiles,
        "selected_route_ids": [r.get("route_id") for r in selected_routes],
        "outputs": {
            "report": str(report_path),
            "visualization_prompt": str(vis_prompt_path),
            "vis_plan_template": str(vis_template_path),
            "schemas": {
                "routes": str(schema_dir / "routes.schema.json"),
                "strategy": str(schema_dir / "strategy.schema.json"),
                "vis_plan": str(schema_dir / "vis_plan.schema.json")
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
            }
        }
    }
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
    return parser.parse_args()


def main() -> int:
    args = parse_args()
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
                "outputs": {
                    "vis_plan_template": str(vis_template_path),
                    "schemas": {
                        "routes": str(schema_dir / "routes.schema.json"),
                        "strategy": str(schema_dir / "strategy.schema.json"),
                        "vis_plan": str(schema_dir / "vis_plan.schema.json")
                    }
                }
            }
            _write_json(out_path, payload)
            print(json.dumps(payload, indent=2, ensure_ascii=False))
            return 0

        _ensure_required_args(args)

        strategy = _read_json(args.strategy_file, {"analysis": {}})
        routes = _read_json(args.routes_file, [])
        vis_plan = _read_json(args.vis_plan_file, {})
        route_ids = [r.get("route_id") for r in routes if isinstance(r, dict) and isinstance(r.get("route_id"), int)]

        validation_errors: List[str] = []
        validation_errors.extend(_validate_strategy(strategy))
        validation_errors.extend(_validate_routes(routes))
        validation_errors.extend(_validate_vis_plan(vis_plan, route_ids))
        if validation_errors:
            raise ValueError("Schema validation failed: " + " | ".join(validation_errors))

        if args.validate_only:
            payload = {
                "success": True,
                "mode": "validate_only",
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
                        "vis_plan": str(schema_dir / "vis_plan.schema.json")
                    }
                }
            }
            _write_json(out_path, payload)
            print(json.dumps(payload, indent=2, ensure_ascii=False))
            return 0

        summary = run_pipeline(args)
        _write_json(out_path, summary)
        print(json.dumps(summary, indent=2, ensure_ascii=False))
        return 0
    except Exception as exc:
        payload = {"success": False, "error": str(exc)}
        _write_json(out_path, payload)
        print(json.dumps(payload, indent=2, ensure_ascii=False))
        return 1


if __name__ == "__main__":
    sys.exit(main())
