"""Visualization + report stage for closed-loop pipeline."""

from pathlib import Path
from typing import Any, Dict, List
import time


def _safe_int_list(values: Any, atom_count: int) -> List[int]:
    if not isinstance(values, list):
        return []
    out: List[int] = []
    for v in values:
        if isinstance(v, int) and 0 <= v < atom_count:
            out.append(v)
    return out


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


def _route_rs_precursors(route_obj: Dict[str, Any]) -> List[str]:
    values: List[str] = []
    raw = route_obj.get("precursors", [])
    if isinstance(raw, list):
        for item in raw:
            if isinstance(item, str) and item.strip():
                values.append(item.strip())
    if values:
        return values
    steps_local = route_obj.get("steps", [])
    if isinstance(steps_local, list) and steps_local:
        first = steps_local[0] if isinstance(steps_local[0], dict) else {}
        rxn = first.get("reaction_smiles") or first.get("rxn_smiles") or first.get("smirks")
        if isinstance(rxn, str) and ">>" in rxn:
            left = rxn.split(">>", 1)[0]
            for frag in left.split("."):
                frag = frag.strip()
                if frag:
                    values.append(frag)
    return values


def render_and_report(
    *,
    args,
    output_dir: Path,
    images_dir: Path,
    vis_plan: Dict[str, Any],
    vis_mod,
    report_mod,
    audited_routes: List[Dict[str, Any]],
    selected_routes: List[Dict[str, Any]],
    rejected_routes: List[Dict[str, Any]],
    strategy: Dict[str, Any],
    audit_summary: Dict[str, Any],
) -> Dict[str, Any]:
    image_paths: Dict[str, str] = {}

    from rdkit import Chem

    t = time.perf_counter()
    routes_for_visualization = audited_routes if args.strict_audit else selected_routes
    target_mol = Chem.MolFromSmiles(args.target_smiles)
    target_atom_count = target_mol.GetNumAtoms() if target_mol else 0
    target_plan = vis_plan.get("target_image", {}) if isinstance(vis_plan, dict) else {}
    target_highlights = _safe_int_list(target_plan.get("highlight_atoms", []), target_atom_count)
    target_legend = target_plan.get("legend", args.target_legend)
    target_path = images_dir / "target.png"
    vis_mod.generate_molecule_image(args.target_smiles, str(target_path), legend=target_legend, highlight_atoms=target_highlights)
    image_paths["target_image"] = "images/target.png"

    for route in routes_for_visualization:
        route_id = route.get("route_id")
        if route_id is None:
            continue

        route_plan = _get_route_plan(vis_plan, route_id) if isinstance(vis_plan, dict) else {}
        captions = _step_caption_map(route_plan)

        precursors = _route_rs_precursors(route)
        precursor_image_files: List[str] = []
        if isinstance(precursors, list):
            for idx, smi in enumerate(precursors, start=1):
                out_name = f"route_{route_id}_precursor_{idx}.png"
                out_path = images_dir / out_name
                size = vis_mod.suggest_molecule_image_size(smi) if hasattr(vis_mod, "suggest_molecule_image_size") else (520, 420)
                if vis_mod.generate_molecule_image(smi, str(out_path), legend=f"RS {idx}", size=size):
                    image_paths[f"route_{route_id}_precursor_{idx}"] = f"images/{out_name}"
                    precursor_image_files.append(str(out_path))

        if isinstance(precursors, list) and precursors and route_plan.get("overview_grid", {}).get("enabled", True):
            out_name = f"route_{route_id}_overview.png"
            ok = False
            if hasattr(vis_mod, "generate_precursor_target_orthogonal_tree_image"):
                ok = vis_mod.generate_precursor_target_orthogonal_tree_image(
                    args.target_smiles,
                    precursors,
                    str(images_dir / out_name),
                    precursor_image_paths=precursor_image_files,
                )
            if not ok:
                ok = vis_mod.generate_rsps_route_pathway_image(args.target_smiles, route, str(images_dir / out_name))
            if not ok:
                legends = route_plan.get("overview_grid", {}).get("precursor_legends", [])
                vis_mod.generate_route_grid(precursors, legends, str(images_dir / out_name))
            image_paths[f"route_{route_id}_overview"] = f"images/{out_name}"

        if isinstance(precursors, list) and precursors and route_plan.get("tree_view", {}).get("enabled", True):
            out_name = f"route_{route_id}_tree.png"
            ok = False
            if hasattr(vis_mod, "generate_precursor_target_orthogonal_tree_image"):
                ok = vis_mod.generate_precursor_target_orthogonal_tree_image(
                    args.target_smiles,
                    precursors,
                    str(images_dir / out_name),
                    precursor_image_paths=precursor_image_files,
                )
            if not ok:
                ok = vis_mod.generate_rsps_route_pathway_image(args.target_smiles, route, str(images_dir / out_name))
            if not ok:
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
    render_seconds = round(time.perf_counter() - t, 6)

    t = time.perf_counter()
    report_path = output_dir / "RETRO_REPORT.md"
    report_mod.generate_report(
        target_smiles=args.target_smiles,
        strategy_analysis=strategy,
        top_routes=audited_routes if args.strict_audit else selected_routes,
        image_paths=image_paths,
        output_file=str(report_path),
        rejected_routes=rejected_routes,
        audit_summary=audit_summary,
        all_routes_visualized=bool(args.strict_audit),
        planner_loop_summary=getattr(args, "planner_loop_summary", None),
        planner_loop_enabled=bool(
            isinstance(getattr(args, "planner_loop_summary", None), dict)
            and getattr(args, "planner_loop_summary", {}).get("rounds")
        ),
        route_count_policy=getattr(args, "route_count_policy", None),
    )
    report_seconds = round(time.perf_counter() - t, 6)

    return {
        "image_paths": image_paths,
        "report_path": report_path,
        "render_seconds": render_seconds,
        "report_seconds": report_seconds,
    }
