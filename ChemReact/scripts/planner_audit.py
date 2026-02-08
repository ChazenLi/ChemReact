"""Planner + audit utilities for closed-loop pipeline."""

import argparse
import json
import re
import importlib.util
import subprocess
import time
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Tuple
from collections import Counter


MIN_REQUIRED_ROUTES = 3
MAX_AUTO_PROPOSE_ROUTES = 5
MAX_REPAIR_ATTEMPTS = 5
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
                "reaction_smiles": "<REACTION_SMILES>",
            }
        ],
    }
]


def build_planner_request(
    target_smiles: str,
    max_routes: int,
    strategy: Any,
    target_chemistry: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    bounded_routes = max(MIN_REQUIRED_ROUTES, min(int(max_routes), MAX_AUTO_PROPOSE_ROUTES))
    portfolio_note = (
        f"Generate {MIN_REQUIRED_ROUTES} to {MAX_AUTO_PROPOSE_ROUTES} routes with a mix of convergent and "
        f"linear strategies. Include both synthesis styles when target complexity allows."
    )
    payload: Dict[str, Any] = {
        "target_smiles": target_smiles,
        "planner_backend": "host",
        "max_routes": bounded_routes,
        "analysis_first_required": True,
        "required_preplanning_tools": [
            "analyze_structure",
            "calculate_features",
            "functional_groups",
        ],
        "constraints": [
            "First perform detailed chemistry analysis of the target molecule before planning routes.",
            "Return strict JSON only.",
            "Each route must include route_id, score, audit_verdict, critical_issues, steps.",
            "Use reaction_smiles/rxn_smiles/smirks for step-level renderability.",
            "All routes must be chemically plausible and strategy-consistent.",
            "Route-level topological disconnections must be derived from target_chemistry and functional_groups.",
            "Target ring/heterocycle topology must be preserved or explicitly justified by ring-forming/heterocycle-forming steps.",
            "Final step product must equal target_smiles.",
            "Every molecule (target, precursors, step reactants, step products) must remain analyzable and chemically valid.",
            "Every reaction step must pass chemistry feasibility checks: atom conservation (including H), bond-change plausibility, class-mechanism consistency.",
            "If atom balance requires leaving-group products, include explicit byproducts in reaction_smiles.",
            "Each route must include synthesis_style=linear|convergent.",
            portfolio_note,
        ],
    }
    if isinstance(strategy, dict) and strategy:
        payload["strategy"] = strategy
    if isinstance(target_chemistry, dict) and target_chemistry:
        payload["target_chemistry"] = target_chemistry
    return payload


def build_planner_prompt(planner_request: Dict[str, Any]) -> str:
    return (
        "You are the host LLM retrosynthesis planner.\n"
        "You must first analyze the input molecule with available chemistry tools and extract detailed chemistry context.\n"
        "Then generate candidate routes that satisfy the JSON contract and reflect that analysis.\n"
        "Execution contract: analysis_first_required=true and required_preplanning_tools are mandatory.\n"
        "CRITICAL: Read and use 'target_chemistry' in the request as hard planning context. "
        "Respect functional groups, ring topology, key disconnections, and risk constraints.\n"
        "Return a JSON array only (no markdown, no explanation, no planner-option exploration).\n\n"
        f"Planner request:\n{json.dumps(planner_request, indent=2, ensure_ascii=False)}\n\n"
        f"Target output template:\n{json.dumps(PLANNER_ROUTES_TEMPLATE, indent=2, ensure_ascii=False)}\n"
    )


def resolve_planner_routes_candidate(args: argparse.Namespace, output_dir: Path) -> Path:
    if args.planner_routes_file:
        return Path(args.planner_routes_file).resolve()
    return (output_dir / "planner_routes.json").resolve()


def render_host_command(command_template: str, values: Dict[str, Any]) -> str:
    try:
        return command_template.format(**values)
    except KeyError as exc:
        raise ValueError(f"--host-planner-command missing placeholder value: {exc}") from exc


def invoke_host_planner(
    command_template: str,
    values: Dict[str, Any],
    log_path: Path,
    write_json: Callable[[Path, Any], None],
) -> Dict[str, Any]:
    command_text = render_host_command(command_template, values)
    proc = subprocess.run(command_text, shell=True, capture_output=True, text=True, check=False)
    payload = {
        "command": command_text,
        "returncode": proc.returncode,
        "stdout": proc.stdout,
        "stderr": proc.stderr,
    }
    write_json(log_path, payload)
    return payload


def planner_request_payload(
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
    outputs: Dict[str, Any] = {
        "planner_request": str(planner_request_path),
        "planner_routes_template": str(planner_routes_template_path),
        "planner_prompt": str(planner_prompt_path),
        "schemas": {
            "planner_request": str(schema_dir / "planner_request.schema.json"),
            "planner_routes": str(schema_dir / "planner_routes.schema.json"),
        },
        "planner_routes_expected": str(planner_routes_candidate_path),
    }
    if extra_outputs:
        outputs.update(extra_outputs)
    return {
        "success": True,
        "mode": mode,
        "timings": {"total_seconds": round(time.perf_counter() - main_start, 6)},
        "outputs": outputs,
        "next_action": next_action,
    }


def _route_precursors(route: Dict[str, Any]) -> List[str]:
    values: List[str] = []
    for item in route.get("precursors", []):
        if isinstance(item, str) and item.strip():
            values.append(item.strip())
    for step in route.get("steps", []):
        if not isinstance(step, dict):
            continue
        rxn = step.get("reaction_smiles") or step.get("rxn_smiles") or step.get("smirks")
        if isinstance(rxn, str) and ">>" in rxn:
            left = rxn.split(">>", 1)[0]
            for frag in left.split("."):
                frag = frag.strip()
                if frag:
                    values.append(frag)
    dedup: List[str] = []
    seen = set()
    for value in values:
        if value not in seen:
            seen.add(value)
            dedup.append(value)
    return dedup


def _extract_invalid_smiles_hints(rejected_routes: List[Dict[str, Any]]) -> List[str]:
    hints: List[str] = []
    patterns = [
        r"Invalid molecule in reaction side:\s*([^\s,;]+)",
        r"mapped reaction is invalid:\s*([^\s,;]+)",
        r"final step product does not match target molecule:\s*([^\s,;]+)",
    ]
    for route in rejected_routes:
        if not isinstance(route, dict):
            continue
        audit = route.get("audit_details", {})
        if not isinstance(audit, dict):
            continue
        issue_messages: List[str] = []
        for key in ("hard_issues", "soft_issues"):
            values = audit.get(key, [])
            if isinstance(values, list):
                issue_messages.extend([v for v in values if isinstance(v, str)])
        for step in audit.get("steps", []):
            if not isinstance(step, dict):
                continue
            for key in ("hard_issues", "soft_issues"):
                values = step.get(key, [])
                if isinstance(values, list):
                    issue_messages.extend([v for v in values if isinstance(v, str)])
        for msg in issue_messages:
            for pattern in patterns:
                match = re.search(pattern, msg)
                if match:
                    smi = match.group(1).strip()
                    if smi and smi not in hints:
                        hints.append(smi)
    return hints


def _compact_route(route: Dict[str, Any]) -> Dict[str, Any]:
    steps = route.get("steps", [])
    compact_steps: List[Dict[str, Any]] = []
    if isinstance(steps, list):
        for idx, step in enumerate(steps, start=1):
            if not isinstance(step, dict):
                continue
            compact_steps.append(
                {
                    "step_index": idx,
                    "reaction_class": step.get("reaction_class"),
                    "reagents": step.get("reagents", []),
                    "conditions": step.get("conditions"),
                    "reaction_smiles": step.get("reaction_smiles") or step.get("rxn_smiles") or step.get("smirks"),
                }
            )

    audit = route.get("audit_details", {})
    compact_audit = {}
    if isinstance(audit, dict):
        compact_audit = {
            "issue_codes": audit.get("issue_codes", []),
            "hard_issues": audit.get("hard_issues", []),
            "soft_issues": audit.get("soft_issues", []),
        }

    return {
        "route_id": route.get("route_id"),
        "synthesis_style": route.get("synthesis_style"),
        "score": route.get("score"),
        "quality_score": route.get("quality_score"),
        "audit_verdict": route.get("audit_verdict"),
        "critical_issues": route.get("critical_issues", []),
        "precursors": _route_precursors(route),
        "steps": compact_steps,
        "audit_details": compact_audit,
    }


def _coerce_routes_payload(payload: Any) -> Any:
    if isinstance(payload, list):
        return payload
    if isinstance(payload, dict) and isinstance(payload.get("routes"), list):
        return payload["routes"]
    return payload


def run_host_planner_repair_loop(
    *,
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
    write_json: Callable[[Path, Any], None],
    read_json: Callable[[str, Any], Any],
    validate_routes: Callable[[Any], List[str]],
    audit_routes_with_gate: Callable[..., Tuple[List[Dict[str, Any]], List[Dict[str, Any]], List[Dict[str, Any]], Dict[str, Any]]],
    load_module: Callable[[str, Path], Any],
    build_repair_prompt_fn: Callable[[Dict[str, Any]], str],
    require_all_routes_recommendable: bool = True,
    min_recommended_routes: int = 1,
) -> Tuple[Optional[Dict[str, Any]], Optional[Path], Dict[str, Any]]:
    planner_prompt_path.write_text(build_planner_prompt(planner_request), encoding="utf-8")
    write_json(planner_request_path, planner_request)
    write_json(planner_routes_template_path, PLANNER_ROUTES_TEMPLATE)

    rounds = max(1, min(int(args.planner_max_iterations), MAX_REPAIR_ATTEMPTS))
    host_command = (args.host_planner_command or "").strip()
    records: List[Dict[str, Any]] = []
    current_routes_path = planner_routes_candidate_path

    if not current_routes_path.exists():
        if not host_command:
            payload = planner_request_payload(
                main_start=main_start,
                planner_request_path=planner_request_path,
                planner_routes_template_path=planner_routes_template_path,
                planner_prompt_path=planner_prompt_path,
                planner_routes_candidate_path=planner_routes_candidate_path,
                schema_dir=schema_dir,
                mode="planner_request_only",
                next_action="Have the host LLM generate and execute the full host_planner_command to produce planner_routes_expected.",
            )
            return payload, None, {"rounds": records}
        log_path = output_dir / "planner_host_invocation_round_1.json"
        invoke_result = invoke_host_planner(
            host_command,
            {
                "prompt_file": str(planner_prompt_path),
                "output_file": str(current_routes_path),
                "request_file": str(planner_request_path),
                "output_dir": str(output_dir),
                "stage": "plan",
                "round": 1,
            },
            log_path=log_path,
            write_json=write_json,
        )
        if invoke_result["returncode"] != 0:
            raise RuntimeError(f"Host planner command failed at round 1: {invoke_result['stderr']}")
        if not current_routes_path.exists():
            payload = planner_request_payload(
                main_start=main_start,
                planner_request_path=planner_request_path,
                planner_routes_template_path=planner_routes_template_path,
                planner_prompt_path=planner_prompt_path,
                planner_routes_candidate_path=planner_routes_candidate_path,
                schema_dir=schema_dir,
                mode="planner_generation_pending",
                next_action="Host planner command executed but planner_routes_expected is missing. Regenerate and execute a full host_planner_command from host LLM.",
                extra_outputs={"planner_host_log_round_1": str(log_path)},
            )
            return payload, None, {"rounds": records}

    skill_root = Path(__file__).resolve().parents[1]
    atom_mappers_dir = skill_root / "internal" / "atom_mappers"
    rdkit_utils_dir = skill_root / "internal" / "rdkit_utils"
    rxn_mapper_mod = load_module("rxn_mapper_loop", atom_mappers_dir / "rxn_mapper.py")
    atom_utils_mod = load_module("atom_utils_loop", atom_mappers_dir / "utils.py")
    analyze_mod = load_module("analyze_structure_loop", rdkit_utils_dir / "analyze_structure.py")
    features_mod = load_module("calculate_features_loop", rdkit_utils_dir / "calculate_features.py")
    fg_mod = load_module("functional_groups_loop", rdkit_utils_dir / "functional_groups.py")

    for round_idx in range(1, rounds + 1):
        routes = _coerce_routes_payload(read_json(str(current_routes_path), []))
        schema_errors = validate_routes(routes)
        if schema_errors:
            audited_routes = routes if isinstance(routes, list) else []
            recommendable_routes: List[Dict[str, Any]] = []
            rejected_routes = audited_routes
            audit_summary = {
                "strict_mode": True,
                "total_routes": len(audited_routes),
                "recommended_routes": 0,
                "rejected_routes": len(audited_routes),
                "verdict_counts": {"FAIL": len(audited_routes)},
                "schema_errors": schema_errors,
            }
        else:
            audited_routes, recommendable_routes, rejected_routes, audit_summary = audit_routes_with_gate(
                routes=routes,
                target_smiles=args.target_smiles,
                strict_audit=True,
                rxn_mapper_mod=rxn_mapper_mod,
                atom_utils_mod=atom_utils_mod,
                strategy=strategy_for_planner if isinstance(strategy_for_planner, dict) else {},
                analyze_mod=analyze_mod,
                features_mod=features_mod,
                fg_mod=fg_mod,
                enforce_style_mix=True,
            )

        audited_path = output_dir / f"planner_routes.audited.round_{round_idx}.json"
        write_json(audited_path, audited_routes)
        round_record = {
            "round": round_idx,
            "input_routes": str(current_routes_path),
            "audited_routes": str(audited_path),
            "audit": audit_summary,
            "recommended_route_ids": [r.get("route_id") for r in recommendable_routes],
            "rejected_route_ids": [r.get("route_id") for r in rejected_routes],
        }
        current_quality = float(audit_summary.get("best_quality_score", 0.0) or 0.0)
        round_record["best_quality_score"] = current_quality
        round_record["requires_repair"] = bool(rejected_routes) if require_all_routes_recommendable else (not bool(recommendable_routes))
        records.append(round_record)

        enough_recommended = len(recommendable_routes) >= max(1, int(min_recommended_routes))
        if enough_recommended and (not require_all_routes_recommendable or not rejected_routes):
            return None, audited_path, {
                "rounds": records,
                "final_round": round_idx,
                "max_iterations": rounds,
                "termination_reason": "recommendable_routes_found",
                "reached_repair_limit": False,
            }
        if round_idx >= rounds:
            return None, audited_path, {
                "rounds": records,
                "final_round": round_idx,
                "max_iterations": rounds,
                "termination_reason": "max_iterations_reached",
                "reached_repair_limit": True,
            }

        repair_request = {
            "target_smiles": args.target_smiles,
            "round": round_idx,
            "max_repair_attempts": rounds,
            "audit_summary": audit_summary,
            "strategy": strategy_for_planner,
            "planner_constraints": planner_request.get("constraints", []),
            "target_chemistry": planner_request.get("target_chemistry", {}),
            "candidate_routes": [_compact_route(r) for r in audited_routes if isinstance(r, dict)],
            "candidate_precursors": [
                {
                    "route_id": r.get("route_id"),
                    "precursors": _route_precursors(r),
                }
                for r in audited_routes
                if isinstance(r, dict)
            ],
            "invalid_smiles_hints": _extract_invalid_smiles_hints(rejected_routes),
            "previous_rounds": records,
            "rejected_routes": rejected_routes,
            "instruction": "Repair failed routes and return a corrected planner_routes JSON array.",
        }
        repair_request_path = output_dir / f"planner_repair_request.round_{round_idx}.json"
        repair_prompt_path = output_dir / f"planner_repair_prompt.round_{round_idx}.txt"
        write_json(repair_request_path, repair_request)
        repair_prompt_path.write_text(build_repair_prompt_fn(repair_request), encoding="utf-8")

        next_routes_path = output_dir / f"planner_routes.repaired.round_{round_idx + 1}.json"
        if not host_command:
            payload = planner_request_payload(
                main_start=main_start,
                planner_request_path=planner_request_path,
                planner_routes_template_path=planner_routes_template_path,
                planner_prompt_path=planner_prompt_path,
                planner_routes_candidate_path=next_routes_path,
                schema_dir=schema_dir,
                mode="planner_repair_requested",
                next_action="Have the host LLM generate and execute a full repair-stage host_planner_command using planner_repair_prompt, then rerun.",
                extra_outputs={
                    "planner_repair_request": str(repair_request_path),
                    "planner_repair_prompt": str(repair_prompt_path),
                    "planner_loop": {"rounds": records, "max_iterations": rounds, "reached_repair_limit": False},
                },
            )
            return payload, None, {"rounds": records}

        log_path = output_dir / f"planner_host_invocation_round_{round_idx + 1}.json"
        invoke_result = invoke_host_planner(
            host_command,
            {
                "prompt_file": str(repair_prompt_path),
                "output_file": str(next_routes_path),
                "request_file": str(repair_request_path),
                "output_dir": str(output_dir),
                "stage": "repair",
                "round": round_idx + 1,
            },
            log_path=log_path,
            write_json=write_json,
        )
        if invoke_result["returncode"] != 0:
            raise RuntimeError(f"Host planner command failed at round {round_idx + 1}: {invoke_result['stderr']}")
        if not next_routes_path.exists():
            payload = planner_request_payload(
                main_start=main_start,
                planner_request_path=planner_request_path,
                planner_routes_template_path=planner_routes_template_path,
                planner_prompt_path=planner_prompt_path,
                planner_routes_candidate_path=next_routes_path,
                schema_dir=schema_dir,
                mode="planner_repair_pending",
                next_action="Host planner command returned but repaired routes file is missing. Regenerate and execute a full repair-stage host_planner_command from host LLM.",
                extra_outputs={
                    "planner_repair_request": str(repair_request_path),
                    "planner_repair_prompt": str(repair_prompt_path),
                    "planner_host_log": str(log_path),
                    "planner_loop": {"rounds": records, "max_iterations": rounds, "reached_repair_limit": False},
                },
            )
            return payload, None, {"rounds": records}
        current_routes_path = next_routes_path

    return None, current_routes_path, {
        "rounds": records,
        "final_round": len(records),
        "max_iterations": rounds,
        "termination_reason": "loop_exhausted",
        "reached_repair_limit": len(records) >= rounds,
    }

# --- migrated from core_pipeline.py: audit helpers + audit gate ---
def _load_module(module_name: str, file_path: Path):
    spec = importlib.util.spec_from_file_location(module_name, str(file_path))
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot load module from {file_path}")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


_REACTION_CLASSIFIER_MOD = None
_NORMALIZE_MOD = None
_ANALYSIS_MOD = None


def _get_reaction_classifier_mod():
    global _REACTION_CLASSIFIER_MOD
    if _REACTION_CLASSIFIER_MOD is None:
        skill_root = Path(__file__).resolve().parents[1]
        classifier_path = skill_root / "internal" / "reaction_classifier" / "classifier.py"
        _REACTION_CLASSIFIER_MOD = _load_module("reaction_classifier_loop", classifier_path)
    return _REACTION_CLASSIFIER_MOD


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
        # Load repair_guidance directly as it contains core functions used
        repair_guidance_path = skill_root / "internal" / "analysis" / "repair_guidance.py"
        _ANALYSIS_MOD = _load_module("repair_guidance_loop", repair_guidance_path)
    return _ANALYSIS_MOD


_PRECURSOR_ANALYZER_MOD = None
_REACTION_VALIDATOR_MOD = None


def _get_precursor_analyzer_mod():
    """Lazy load the precursor analyzer module."""
    global _PRECURSOR_ANALYZER_MOD
    if _PRECURSOR_ANALYZER_MOD is None:
        skill_root = Path(__file__).resolve().parents[1]
        path = skill_root / "internal" / "analysis" / "precursor_analyzer.py"
        _PRECURSOR_ANALYZER_MOD = _load_module("precursor_analyzer_loop", path)
    return _PRECURSOR_ANALYZER_MOD


def _get_reaction_validator_mod():
    """Lazy load the reaction validator module."""
    global _REACTION_VALIDATOR_MOD
    if _REACTION_VALIDATOR_MOD is None:
        skill_root = Path(__file__).resolve().parents[1]
        path = skill_root / "internal" / "analysis" / "reaction_validator.py"
        _REACTION_VALIDATOR_MOD = _load_module("reaction_validator_loop", path)

def _canonicalize_smiles(smiles: str):
    from rdkit import Chem

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    for atom in mol.GetAtoms():
        if atom.GetAtomMapNum():
            atom.SetAtomMapNum(0)
    return Chem.MolToSmiles(mol, canonical=True)


def _parse_reaction_sides(rxn_text: str) -> Tuple[Optional[str], Optional[str], Optional[str]]:
    if not isinstance(rxn_text, str) or ">>" not in rxn_text:
        return None, None, "Reaction must contain '>>' separator"
    left, right = rxn_text.split(">>", 1)
    left = left.strip()
    right = right.strip()
    if not left or not right:
        return None, None, "Reaction has empty reactant or product side"
    return left, right, None


def _element_counter_from_side(side: str) -> Tuple[Optional[Counter], Optional[str]]:
    from rdkit import Chem

    counts: Counter = Counter()
    for frag in side.split("."):
        frag = frag.strip()
        if not frag:
            continue
        mol = Chem.MolFromSmiles(frag)
        if mol is None:
            return None, f"Invalid molecule in reaction side: {frag}"
        # Include implicit hydrogens in atom conservation by expanding them explicitly.
        mol = Chem.AddHs(mol)
        for atom in mol.GetAtoms():
            counts[atom.GetSymbol()] += 1
    return counts, None


def _classify_reaction_type(left_side: str, right_side: str, reaction_class: str) -> Dict[str, Any]:
    classifier_mod = _get_reaction_classifier_mod()
    return classifier_mod.classify_reaction(
        left_side=left_side,
        right_side=right_side,
        reaction_class=reaction_class,
    )


def _infer_missing_byproducts_for_delta(
    reaction_class: str,
    delta: Counter,
    reaction_type_info: Optional[Dict[str, Any]] = None,
) -> List[Dict[str, Any]]:
    classifier_mod = _get_reaction_classifier_mod()
    inferred = classifier_mod.infer_missing_byproducts_for_delta(
        reaction_class=reaction_class,
        delta=delta,
        reaction_type_info=reaction_type_info,
        max_byproducts=2,
    )
    if inferred:
        return inferred
    return _infer_common_small_byproducts(delta, max_byproducts=2)


def _infer_common_small_byproducts(delta: Counter, max_byproducts: int = 2) -> List[Dict[str, Any]]:
    """Fallback byproduct inference for common omitted small molecules."""
    if not delta or any(v < 0 for v in delta.values()):
        return []
    library = [
        {"smiles": "O", "formula": Counter({"H": 2, "O": 1}), "confidence": 0.35},
        {"smiles": "Cl", "formula": Counter({"H": 1, "Cl": 1}), "confidence": 0.30},
        {"smiles": "Br", "formula": Counter({"H": 1, "Br": 1}), "confidence": 0.30},
        {"smiles": "I", "formula": Counter({"H": 1, "I": 1}), "confidence": 0.30},
        {"smiles": "N", "formula": Counter({"N": 1, "H": 3}), "confidence": 0.28},
        {"smiles": "O=C=O", "formula": Counter({"C": 1, "O": 2}), "confidence": 0.25},
        {"smiles": "[H][H]", "formula": Counter({"H": 2}), "confidence": 0.22},
        {"smiles": "CO", "formula": Counter({"C": 1, "H": 4, "O": 1}), "confidence": 0.20},
        {"smiles": "CCO", "formula": Counter({"C": 2, "H": 6, "O": 1}), "confidence": 0.20},
    ]

    best_combo: List[Dict[str, Any]] = []
    best_score = -1.0
    for item in library:
        if item["formula"] == delta and item["confidence"] > best_score:
            best_combo = [item]
            best_score = float(item["confidence"])

    if max_byproducts >= 2:
        for item1 in library:
            for item2 in library:
                if item1["formula"] + item2["formula"] != delta:
                    continue
                score = float(item1["confidence"]) + float(item2["confidence"])
                if score > best_score:
                    best_combo = [item1, item2]
                    best_score = score

    return [
        {
            "smiles": item["smiles"],
            "confidence": round(float(item["confidence"]), 3),
            "source": "common_small_molecule_fallback",
        }
        for item in best_combo
    ]


def _is_likely_small_molecule_imbalance(
    missing_product_delta: Counter,
    excess_product_delta: Counter,
) -> bool:
    if excess_product_delta:
        return False
    if not missing_product_delta:
        return False
    allowed = {"H", "O", "N", "Cl", "Br", "I", "S", "P", "F", "Na", "K", "Li", "C"}
    if not set(missing_product_delta.keys()).issubset(allowed):
        return False
    total = sum(abs(v) for v in missing_product_delta.values())
    carbon = abs(missing_product_delta.get("C", 0))
    return total <= 12 and carbon <= 2


def _normalized_side_smiles_set(side: str) -> Tuple[Optional[List[str]], Optional[str]]:
    values: List[str] = []
    for frag in side.split("."):
        frag = frag.strip()
        if not frag:
            continue
        norm = _canonicalize_smiles(frag)
        if norm is None:
            return None, f"Invalid molecule in reaction side: {frag}"
        values.append(norm)
    return values, None


def _looks_like_mapped_reaction(rxn_text: str) -> bool:
    return bool(re.search(r":\d+\]", rxn_text))


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


def _route_ring_topology_issues(
    target_profile: Dict[str, Any],
    precursor_profiles: List[Dict[str, Any]],
    step_audits: List[Dict[str, Any]],
    route: Dict[str, Any],
) -> Tuple[List[str], Dict[str, Any]]:
    issues: List[str] = []
    target_sizes = set(target_profile.get("ring_sizes", [])) if isinstance(target_profile, dict) else set()
    target_ring_count = int(target_profile.get("ring_count", 0) or 0) if isinstance(target_profile, dict) else 0
    target_heterocycle_count = int(target_profile.get("heterocycle_count", 0) or 0) if isinstance(target_profile, dict) else 0

    precursor_sizes: set = set()
    precursor_ring_count_total = 0
    precursor_heterocycle_total = 0
    for p in precursor_profiles:
        if not isinstance(p, dict):
            continue
        precursor_sizes.update(p.get("ring_sizes", []))
        precursor_ring_count_total += int(p.get("ring_count", 0) or 0)
        precursor_heterocycle_total += int(p.get("heterocycle_count", 0) or 0)

    class_text = " ".join(
        (s.get("reaction_class", "") or "").lower() for s in route.get("steps", []) if isinstance(s, dict)
    )
    ring_forming_keywords = ("cycl", "annulation", "ring")
    hetero_keywords = ("aza", "oxa", "thia", "heterocycle", "hetero")
    has_ring_forming_step = any(
        int(step.get("side_metrics", {}).get("ring_delta", 0) or 0) > 0 for step in step_audits
    ) or any(k in class_text for k in ring_forming_keywords)
    has_hetero_introduction_hint = any(k in class_text for k in hetero_keywords)

    missing_sizes = sorted([x for x in target_sizes if x not in precursor_sizes])
    if target_ring_count > 0 and missing_sizes and not has_ring_forming_step:
        issues.append(
            f"Precursor ring-size set cannot explain target ring topology (missing ring sizes: {missing_sizes})"
        )

    if target_ring_count > 0 and precursor_ring_count_total == 0 and not has_ring_forming_step:
        issues.append("Precursor set has no rings but route has no ring-forming step for target topology")

    if target_heterocycle_count > 0 and precursor_heterocycle_total == 0 and not has_hetero_introduction_hint:
        issues.append("Precursor heterocycle topology cannot explain target heterocycle requirements")

    diag = {
        "target": target_profile,
        "precursors": {
            "ring_sizes_union": sorted(precursor_sizes),
            "ring_count_total": precursor_ring_count_total,
            "heterocycle_count_total": precursor_heterocycle_total,
        },
        "missing_target_ring_sizes": missing_sizes,
        "has_ring_forming_step": has_ring_forming_step,
        "has_hetero_introduction_hint": has_hetero_introduction_hint,
    }
    return issues, diag


def _parse_side_mols(side: str):
    from rdkit import Chem

    mols = []
    for frag in side.split("."):
        frag = frag.strip()
        if not frag:
            continue
        mol = Chem.MolFromSmiles(frag)
        if mol is None:
            return [], f"Invalid molecule in reaction side: {frag}"
        mols.append(mol)
    return mols, None


def _side_metrics(side: str) -> Dict[str, int]:
    from rdkit import Chem

    mols, _ = _parse_side_mols(side)
    return {
        "molecule_count": len(mols),
        "total_formal_charge": int(sum(Chem.GetFormalCharge(m) for m in mols)) if mols else 0,
        "total_rings": int(sum(m.GetRingInfo().NumRings() for m in mols)) if mols else 0,
        "heavy_atoms": int(sum(m.GetNumHeavyAtoms() for m in mols)) if mols else 0,
    }


def _reaction_class_consistency_issues(reaction_class: str, bond_changes: Dict[str, List[Dict[str, Any]]]) -> List[str]:
    issues: List[str] = []
    name = (reaction_class or "").lower()
    formed = len(bond_changes.get("formed", []))
    broken = len(bond_changes.get("broken", []))

    if any(k in name for k in ["coupling", "alkylation", "arylation", "acylation"]) and formed < 1:
        issues.append("Reaction class suggests bond formation but no formed bond detected")
    if any(k in name for k in ["substitution", "hydrolysis", "deprotection", "cleavage"]) and (formed < 1 or broken < 1):
        issues.append("Reaction class suggests exchange/cleavage but formed/broken bond pattern is inconsistent")
    if "protection" in name and broken > 0:
        issues.append("Protection step unexpectedly shows bond cleavage")
    return issues


def _murcko_scaffold(smi: str) -> Optional[str]:
    try:
        from rdkit import Chem
        from rdkit.Chem.Scaffolds import MurckoScaffold

        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            return None
        scf = MurckoScaffold.GetScaffoldForMol(mol)
        if scf is None:
            return None
        return Chem.MolToSmiles(scf, canonical=True)
    except Exception:
        return None


def _strategy_alignment_issues(
    strategy: Dict[str, Any],
    route: Dict[str, Any],
    step_audits: List[Dict[str, Any]],
) -> List[str]:
    issues: List[str] = []
    analysis = strategy.get("analysis", {}) if isinstance(strategy, dict) else {}
    if not isinstance(analysis, dict):
        return issues

    key_disconnections = analysis.get("key_disconnections", [])
    if isinstance(key_disconnections, list) and key_disconnections:
        step_classes = " ".join(
            (s.get("reaction_class", "") or "").lower() for s in route.get("steps", []) if isinstance(s, dict)
        )
        matched = 0
        for item in key_disconnections:
            if isinstance(item, str) and item.strip() and item.lower() in step_classes:
                matched += 1
        if matched == 0:
            issues.append("Route reaction classes do not align with strategy key_disconnections")

    strategy_type = str(analysis.get("strategy_type", "")).lower()
    if strategy_type == "convergent":
        precursor_count = len(route.get("precursors", [])) if isinstance(route.get("precursors", []), list) else 0
        if precursor_count < 2:
            issues.append("Convergent strategy expects multi-precursor assembly")
    if strategy_type == "linear":
        if len(step_audits) > 3:
            issues.append("Linear strategy route has excessive step count")

    return issues


def _route_quality(
    route: Dict[str, Any],
    route_hard: List[str],
    route_soft: List[str],
    strategy_issues: List[str],
) -> Dict[str, Any]:
    base = float(route.get("score", 0.0) or 0.0)
    hard_penalty = 3.0 * len(route_hard)
    soft_penalty = 1.0 * len(route_soft)
    strategy_penalty = 1.5 * len(strategy_issues)
    quality = max(0.0, min(10.0, base - hard_penalty - soft_penalty - strategy_penalty))
    return {
        "base_score": base,
        "hard_penalty": hard_penalty,
        "soft_penalty": soft_penalty,
        "strategy_penalty": strategy_penalty,
        "quality_score": round(quality, 3),
    }


def _route_style(route: Dict[str, Any]) -> str:
    style = route.get("synthesis_style")
    if isinstance(style, str):
        norm = style.strip().lower()
        if norm in {"linear", "convergent"}:
            return norm
    precursors = route.get("precursors", [])
    if isinstance(precursors, list) and len(precursors) >= 2:
        return "convergent"
    return "linear"


def _portfolio_style_counts(routes: List[Dict[str, Any]]) -> Dict[str, int]:
    counts = {"linear": 0, "convergent": 0}
    for r in routes:
        style = _route_style(r)
        if style in counts:
            counts[style] += 1
    return counts


def _audit_single_step(
    step: Dict[str, Any],
    step_index: int,
    route_precursors_norm: List[str],
    mapper,
    atom_utils_mod,
) -> Dict[str, Any]:
    rxn = step.get("reaction_smiles") or step.get("rxn_smiles") or step.get("smirks")
    reaction_class = step.get("reaction_class", "")
    hard_issues: List[str] = []
    soft_issues: List[str] = []
    reactants_norm: List[str] = []
    products_norm: List[str] = []
    side_metrics: Dict[str, Any] = {}
    reaction_type_info: Dict[str, Any] = {
        "id": "unknown",
        "name": "Unknown",
        "confidence": 0.0,
        "source": "none",
        "matched_features": [],
        "candidate_byproducts": [],
        "byproduct_candidates": [],
    }
    atom_balance: Dict[str, Any] = {
        "balanced": None,
        "reactants": {},
        "products": {},
        "missing_product_delta": {},
        "excess_product_delta": {},
        "inferred_byproducts": [],
        "inferred_byproduct_details": [],
        "inferred_byproduct_confidence": 0.0,
    }
    mapping_summary: Dict[str, Any] = {
        "used": False,
        "success": False,
        "confidence": 0.0,
        "reaction_center": [],
        "completeness": 0.0,
        "bond_changes": {"formed": [], "broken": []},
        "notes": [],
    }

    if not isinstance(rxn, str) or not rxn.strip():
        hard_issues.append("Missing reaction_smiles/rxn_smiles/smirks")
        return {
            "step_index": step_index,
            "reaction": rxn or "",
            "hard_issues": hard_issues,
            "soft_issues": soft_issues,
            "reactants_norm": reactants_norm,
            "products_norm": products_norm,
            "reaction_type": reaction_type_info,
            "atom_balance": atom_balance,
            "mapping": mapping_summary,
            "side_metrics": side_metrics,
        }

    left, right, parse_err = _parse_reaction_sides(rxn.strip())
    if parse_err:
        hard_issues.append(parse_err)
        return {
            "step_index": step_index,
            "reaction": rxn,
            "hard_issues": hard_issues,
            "soft_issues": soft_issues,
            "reactants_norm": reactants_norm,
            "products_norm": products_norm,
            "reaction_type": reaction_type_info,
            "atom_balance": atom_balance,
            "mapping": mapping_summary,
            "side_metrics": side_metrics,
        }

    reaction_type_info = _classify_reaction_type(left, right, str(reaction_class))

    reactant_counts, count_err = _element_counter_from_side(left)
    if count_err:
        hard_issues.append(count_err)
    product_counts, count_err = _element_counter_from_side(right)
    if count_err:
        hard_issues.append(count_err)
    if reactant_counts is not None and product_counts is not None and reactant_counts != product_counts:
        missing_product_delta = reactant_counts - product_counts
        excess_product_delta = product_counts - reactant_counts
        inferred_byproducts: List[Dict[str, Any]] = []
        if not excess_product_delta:
            inferred_byproducts = _infer_missing_byproducts_for_delta(
                str(reaction_class),
                missing_product_delta,
                reaction_type_info=reaction_type_info,
            )
        if inferred_byproducts:
            inferred_smiles = [str(item.get("smiles")) for item in inferred_byproducts if isinstance(item, dict)]
            inferred_conf = round(
                max(
                    [float(item.get("confidence", 0.0) or 0.0) for item in inferred_byproducts if isinstance(item, dict)]
                    or [0.0]
                ),
                3,
            )
            soft_issues.append(
                "Atom conservation satisfied after inferred byproducts: " + ".".join(inferred_smiles)
            )
            mapping_summary["notes"].append("Inferred byproducts for atom balance: " + ".".join(inferred_smiles))
            atom_balance = {
                "balanced": True,
                "reactants": dict(reactant_counts),
                "products": dict(product_counts),
                "missing_product_delta": dict(missing_product_delta),
                "excess_product_delta": dict(excess_product_delta),
                "inferred_byproducts": inferred_smiles,
                "inferred_byproduct_details": inferred_byproducts,
                "inferred_byproduct_confidence": inferred_conf,
            }
        elif _is_likely_small_molecule_imbalance(missing_product_delta, excess_product_delta):
            soft_issues.append(
                "Atom conservation imbalance likely due to omitted small-molecule byproducts; provide explicit byproducts in mapped rs>>ps for repair"
            )
            atom_balance = {
                "balanced": False,
                "reactants": dict(reactant_counts),
                "products": dict(product_counts),
                "missing_product_delta": dict(missing_product_delta),
                "excess_product_delta": dict(excess_product_delta),
                "inferred_byproducts": [],
                "inferred_byproduct_details": [],
                "inferred_byproduct_confidence": 0.0,
            }
        else:
            hard_issues.append("Atom conservation failed between reactants and products")
            atom_balance = {
                "balanced": False,
                "reactants": dict(reactant_counts),
                "products": dict(product_counts),
                "missing_product_delta": dict(missing_product_delta),
                "excess_product_delta": dict(excess_product_delta),
                "inferred_byproducts": [],
                "inferred_byproduct_details": [],
                "inferred_byproduct_confidence": 0.0,
            }
    elif reactant_counts is not None and product_counts is not None:
        atom_balance = {
            "balanced": True,
            "reactants": dict(reactant_counts),
            "products": dict(product_counts),
            "missing_product_delta": {},
            "excess_product_delta": {},
            "inferred_byproducts": [],
            "inferred_byproduct_details": [],
            "inferred_byproduct_confidence": 0.0,
        }

    left_metrics = _side_metrics(left)
    right_metrics = _side_metrics(right)
    side_metrics = {
        "reactants": left_metrics,
        "products": right_metrics,
        "ring_delta": right_metrics["total_rings"] - left_metrics["total_rings"],
        "charge_delta": right_metrics["total_formal_charge"] - left_metrics["total_formal_charge"],
    }
    if abs(side_metrics["charge_delta"]) > 1:
        soft_issues.append("Large net formal charge change detected across reaction sides")
    if (
        abs(side_metrics["ring_delta"]) > 1
        and not any(k in str(reaction_class).lower() for k in ["ring", "cycl", "metathesis"])
    ):
        hard_issues.append("Ring-count jump is inconsistent with declared reaction class")

    reactants_norm, norm_err = _normalized_side_smiles_set(left)
    if norm_err:
        hard_issues.append(norm_err)
        reactants_norm = []
    products_norm, norm_err = _normalized_side_smiles_set(right)
    if norm_err:
        hard_issues.append(norm_err)
        products_norm = []

    if step_index == 1 and route_precursors_norm and reactants_norm:
        if not set(reactants_norm).issubset(set(route_precursors_norm)):
            soft_issues.append("Step-1 reactants are not fully covered by route precursors")

    mapped_rxn = None
    if mapper is not None and getattr(mapper, "is_available", lambda: False)():
        mapping_summary["used"] = True
        map_result = mapper.map_reaction(rxn)
        mapping_summary["success"] = bool(map_result.get("success"))
        mapping_summary["confidence"] = float(map_result.get("confidence", 0.0) or 0.0)
        mapping_summary["reaction_center"] = map_result.get("reaction_center", [])
        if mapping_summary["success"] and map_result.get("mapped_rxn"):
            mapped_rxn = map_result.get("mapped_rxn")
        else:
            mapping_summary["notes"].append("RXNMapper failed to map this reaction")
    elif _looks_like_mapped_reaction(rxn):
        mapped_rxn = rxn
        mapping_summary["notes"].append("Using provided atom-mapped reaction text")
    else:
        mapping_summary["notes"].append("Atom mapping unavailable for this step")

    if mapped_rxn:
        mapping_validation = atom_utils_mod.validate_mapped_reaction(mapped_rxn)
        mapping_summary["completeness"] = float(mapping_validation.get("completeness", 0.0) or 0.0)
        if not mapping_validation.get("valid", False):
            hard_issues.append("Mapped reaction is invalid")
        if mapping_summary["completeness"] < 0.99:
            mapping_summary["notes"].append(
                f"Atom mapping completeness is low ({mapping_summary['completeness']:.2f}); include fully mapped rs>>ps for LLM repair/verification"
            )
        mapping_summary["bond_changes"] = atom_utils_mod.get_bond_changes(mapped_rxn)
        if (
            not mapping_summary["bond_changes"].get("formed")
            and not mapping_summary["bond_changes"].get("broken")
        ):
            hard_issues.append("No bond changes detected in mapped reaction")
        soft_issues.extend(_reaction_class_consistency_issues(str(reaction_class), mapping_summary["bond_changes"]))

    # Reject identity/no-op reactions in strict mode: reactant and product sets are identical.
    if reactants_norm and products_norm and set(reactants_norm) == set(products_norm):
        hard_issues.append("Reaction is identity/no-op (RS equals PS), no retrosynthetic transformation provided")

    # Generate selectivity notes based on reaction type
    selectivity_notes: List[str] = []
    reaction_type_name = str(reaction_type_info.get("id", reaction_type_info.get("name", ""))).lower()
    
    if "esterification" in reaction_type_name:
        selectivity_notes.append("Carboxylic acid is more reactive than alcohol")
    if "amide" in reaction_type_name or "coupling" in reaction_type_name:
        selectivity_notes.append("Primary amines more reactive than secondary")
    if "substitution" in reaction_type_name:
        selectivity_notes.append("Steric hindrance affects regioselectivity")
    if "oxidation" in reaction_type_name:
        selectivity_notes.append("Primary alcohols may over-oxidize to carboxylic acids")
    if "reduction" in reaction_type_name:
        selectivity_notes.append("Consider chemoselective reducing agents")
    
    # Normalize reaction SMILES for visualization
    reaction_normalized: Dict[str, str] = {"original": rxn}
    try:
        normalize_mod = _get_normalize_mod()
        reaction_normalized = normalize_mod.normalize_reaction_smiles(rxn)
    except Exception:
        reaction_normalized["display"] = rxn
        reaction_normalized["canonical"] = rxn

    return {
        "step_index": step_index,
        "reaction": rxn,
        "reaction_normalized": reaction_normalized,  # New: normalized reaction SMILES
        "hard_issues": hard_issues,
        "soft_issues": soft_issues,
        "reactants_norm": reactants_norm,
        "products_norm": products_norm,
        "reaction_type": reaction_type_info,
        "atom_balance": atom_balance,
        "mapping": mapping_summary,
        "side_metrics": side_metrics,
        "selectivity_notes": selectivity_notes,  # New: selectivity analysis
    }


def _audit_routes_with_gate(
    routes: List[Dict[str, Any]],
    target_smiles: str,
    strict_audit: bool,
    rxn_mapper_mod,
    atom_utils_mod,
    strategy: Optional[Dict[str, Any]] = None,
    analyze_mod=None,
    features_mod=None,
    fg_mod=None,
    enforce_style_mix: bool = False,
) -> Tuple[List[Dict[str, Any]], List[Dict[str, Any]], List[Dict[str, Any]], Dict[str, Any]]:
    if not strict_audit:
        summary = {
            "strict_mode": False,
            "total_routes": len(routes),
            "recommended_routes": len(routes),
            "rejected_routes": 0,
            "verdict_counts": dict(Counter(r.get("audit_verdict", "UNKNOWN") for r in routes)),
        }
        return routes, routes, [], summary

    target_norm = _canonicalize_smiles(target_smiles)
    target_scaffold = _murcko_scaffold(target_smiles)
    target_ring_profile = _ring_hetero_profile(target_smiles)
    mapper = rxn_mapper_mod.RXNMapperWrapper() if rxn_mapper_mod is not None else None
    analysis_enabled = analyze_mod is not None and features_mod is not None and fg_mod is not None
    target_analysis = _analyze_molecule_bundle(target_smiles, analyze_mod, features_mod, fg_mod) if analysis_enabled else None
    precursor_analysis_cache: Dict[str, Dict[str, Any]] = {}
    molecule_analysis_cache: Dict[str, Dict[str, Any]] = {}
    portfolio_counts = _portfolio_style_counts(routes)
    portfolio_issue = None
    if enforce_style_mix:
        # Flexible constraint: require at least MIN_REQUIRED_ROUTES and both styles.
        if len(routes) < MIN_REQUIRED_ROUTES:
            portfolio_issue = (
                f"Route portfolio must contain at least {MIN_REQUIRED_ROUTES} routes, "
                f"got {len(routes)}"
            )
        elif portfolio_counts["convergent"] == 0 or portfolio_counts["linear"] == 0:
            # Require both styles to be present; exact counts are not fixed.
            portfolio_issue = (
                f"Route portfolio must include both convergent and linear styles. "
                f"Got convergent={portfolio_counts['convergent']}, linear={portfolio_counts['linear']}"
            )

    audited: List[Dict[str, Any]] = []
    recommended: List[Dict[str, Any]] = []
    rejected: List[Dict[str, Any]] = []
    verdict_counts: Counter = Counter()

    for route in routes:
        route_copy = dict(route)
        route_copy["input_audit_verdict"] = route.get("audit_verdict", "UNKNOWN")
        route_copy["synthesis_style"] = _route_style(route_copy)
        route_hard: List[str] = []
        route_soft: List[str] = []
        step_audits: List[Dict[str, Any]] = []
        precursor_analysis: Dict[str, Dict[str, Any]] = {}
        step_molecule_analysis: Dict[str, Dict[str, Any]] = {}
        precursor_ring_profiles: List[Dict[str, Any]] = []
        precursor_compatibility_issues: List[Dict[str, Any]] = []
        ring_topology_diag: Dict[str, Any] = {}
        if portfolio_issue:
            route_hard.append(portfolio_issue)

        precursors = route.get("precursors", [])
        route_precursors_norm: List[str] = []
        if not isinstance(precursors, list) or not precursors:
            route_hard.append("Route precursors missing; RS analysis unavailable")
        if isinstance(precursors, list):
            for smi in precursors:
                if isinstance(smi, str):
                    # Disallow unresolved synthons/dummy atoms in final route precursors.
                    if "*" in smi:
                        route_hard.append(
                            f"Precursor contains unresolved synthon marker '*': {smi}"
                        )
                    norm = _canonicalize_smiles(smi)
                    if norm is not None:
                        route_precursors_norm.append(norm)
                    if analysis_enabled:
                        if smi not in precursor_analysis_cache:
                            precursor_analysis_cache[smi] = _analyze_molecule_bundle(smi, analyze_mod, features_mod, fg_mod)
                        precursor_analysis[smi] = precursor_analysis_cache[smi]
                        if isinstance(precursor_analysis[smi].get("ring_profile"), dict):
                            precursor_ring_profiles.append(precursor_analysis[smi]["ring_profile"])
                        if not precursor_analysis[smi].get("success"):
                            route_hard.append(f"Precursor analysis failed: {smi}")
                        elif not precursor_analysis[smi].get("functional_groups", {}).get("success"):
                            route_hard.append(f"Precursor functional-group analysis missing: {smi}")

            # Use precursor_analyzer for precursor-set compatibility analysis.
            try:
                precursor_analyzer_mod = _get_precursor_analyzer_mod()
                precursor_set_result = precursor_analyzer_mod.analyze_precursor_set(
                    [s for s in precursors if isinstance(s, str)],
                    include_compatibility_check=True
                )
                if precursor_set_result.get("compatibility_issues"):
                    for compat_issue in precursor_set_result["compatibility_issues"]:
                        if isinstance(compat_issue, dict):
                            precursor_compatibility_issues.append(compat_issue)
                        severity = compat_issue.get("severity", "info")
                        code = compat_issue.get("code", "COMPATIBILITY")
                        msg = compat_issue.get("message", "Unknown compatibility issue")
                        guidance = ""
                        if code == "ACID_BASE_COEXIST":
                            guidance = (
                                "Explain salt-formation risk, specify protection/deprotection order, "
                                "or provide staged addition/acid-base control."
                            )
                        elif code == "COMPETING_NUCLEOPHILES":
                            guidance = (
                                "Explain chemoselectivity rationale and why chosen nucleophile dominates."
                            )
                        detail = f"Precursor compatibility [{code}]: {msg}"
                        if guidance:
                            detail += f" Guidance: {guidance}"
                        if severity == "warning":
                            route_soft.append(detail)
                        # info-level findings are recorded but not added to issues
            except Exception:
                pass  # compatibility analysis failure does not block the main flow

        if not analysis_enabled:
            route_hard.append("Target molecule analysis unavailable")
        elif not target_analysis or not target_analysis.get("success"):
            route_hard.append("Target molecule analysis unavailable")

        steps = route.get("steps", [])
        if not isinstance(steps, list) or not steps:
            route_hard.append("Route has no reaction steps")
        else:
            for idx, step in enumerate(steps, start=1):
                if not isinstance(step, dict):
                    route_hard.append(f"Step {idx} is not an object")
                    continue
                step_audit = _audit_single_step(
                    step=step,
                    step_index=idx,
                    route_precursors_norm=route_precursors_norm,
                    mapper=mapper,
                    atom_utils_mod=atom_utils_mod,
                )
                step_audits.append(step_audit)
                route_hard.extend([f"Step {idx}: {msg}" for msg in step_audit["hard_issues"]])
                route_soft.extend([f"Step {idx}: {msg}" for msg in step_audit["soft_issues"]])
                if analysis_enabled:
                    step_smiles = list(dict.fromkeys(step_audit.get("reactants_norm", []) + step_audit.get("products_norm", [])))
                    for smi in step_smiles:
                        if smi not in molecule_analysis_cache:
                            molecule_analysis_cache[smi] = _analyze_molecule_bundle(smi, analyze_mod, features_mod, fg_mod)
                        analyzed = molecule_analysis_cache[smi]
                        step_molecule_analysis[smi] = analyzed
                        if not analyzed.get("success"):
                            route_hard.append(f"Step {idx} molecule analysis failed: {smi}")
                        elif not analyzed.get("functional_groups", {}).get("success"):
                            route_hard.append(f"Step {idx} functional-group analysis missing: {smi}")

            if len(step_audits) >= 2:
                for i in range(1, len(step_audits)):
                    prev_products = set(step_audits[i - 1].get("products_norm", []))
                    curr_reactants = set(step_audits[i].get("reactants_norm", []))
                    if prev_products and curr_reactants and prev_products.isdisjoint(curr_reactants):
                        route_hard.append(f"Step {i} products do not connect to Step {i + 1} reactants")

            if target_norm and step_audits:
                last_products = step_audits[-1].get("products_norm", [])
                if target_norm not in last_products:
                    route_hard.append("Final step product does not match target molecule")
                if target_scaffold and last_products:
                    matched = False
                    for prod in last_products:
                        if _murcko_scaffold(prod) == target_scaffold:
                            matched = True
                            break
                    if not matched:
                        route_hard.append("Final-step scaffold does not preserve target core skeleton")

            # Cross-step ring continuity sanity
            if step_audits:
                ring_deltas = [int(step.get("side_metrics", {}).get("ring_delta", 0) or 0) for step in step_audits]
                if sum(abs(x) for x in ring_deltas) > max(2, len(step_audits)):
                    route_soft.append("Route shows unstable ring-transform profile across steps")

            ring_topology_issues, ring_topology_diag = _route_ring_topology_issues(
                target_profile=target_ring_profile,
                precursor_profiles=precursor_ring_profiles,
                step_audits=step_audits,
                route=route,
            )
            route_hard.extend(ring_topology_issues)

        strategy_issues = _strategy_alignment_issues(strategy or {}, route, step_audits)
        # NOTE: do not append strategy_issues to route_soft to avoid double-penalizing in _route_quality.
        # strategy_issues are handled separately via strategy_penalty in _route_quality.

        risk_budget = {}
        if isinstance(strategy, dict):
            analysis = strategy.get("analysis", {})
            if isinstance(analysis, dict) and isinstance(analysis.get("risk_budget"), dict):
                risk_budget = analysis.get("risk_budget", {})
        if risk_budget:
            max_hard = int(risk_budget.get("max_hard_issues", 0) or 0)
            max_soft = int(risk_budget.get("max_soft_issues", 999999) or 999999)
            if max_hard >= 0 and len(route_hard) > max_hard:
                route_hard.append("Route exceeds strategy risk_budget.max_hard_issues")
            if max_soft >= 0 and len(route_soft) > max_soft:
                route_soft.append("Route exceeds strategy risk_budget.max_soft_issues")

        quality = _route_quality(route, route_hard, route_soft, strategy_issues)
        min_quality = None
        if risk_budget and isinstance(risk_budget.get("min_quality_score"), (int, float)):
            min_quality = float(risk_budget["min_quality_score"])
            if quality["quality_score"] < min_quality:
                route_hard.append("Route quality below strategy risk budget threshold")

        if route_hard:
            verdict = "FAIL"
        elif route_soft:
            verdict = "CONDITIONAL"
        else:
            verdict = "PASS"

        merged_issues: List[str] = []
        existing_issues = route.get("critical_issues", [])
        if isinstance(existing_issues, list):
            for issue in existing_issues:
                if isinstance(issue, str):
                    merged_issues.append(issue)
        merged_issues.extend(route_hard)
        merged_issues.extend(route_soft)

        dedup_issues: List[str] = []
        seen = set()
        for issue in merged_issues:
            if issue not in seen:
                seen.add(issue)
                dedup_issues.append(issue)
        issue_codes = sorted({_issue_code(msg) for msg in dedup_issues})

        route_copy["audit_verdict"] = verdict
        route_copy["critical_issues"] = dedup_issues
        route_copy["quality_score"] = quality["quality_score"]
        route_copy["audit_details"] = {
            "hard_issues": route_hard,
            "soft_issues": route_soft,
            "issue_codes": issue_codes,
            "steps": step_audits,
            "strategy_issues": strategy_issues,
            "risk_budget": risk_budget,
            "quality": quality,
            "target_scaffold": target_scaffold,
            "target_analysis": {
                "success": bool(target_analysis and target_analysis.get("success")),
                "error": None if not target_analysis else target_analysis.get("error"),
                "functional_groups": {} if not target_analysis else target_analysis.get("functional_groups", {}).get("functional_groups", {}),
                "ring_profile": target_ring_profile,
            },
            "precursor_analysis": {
                smi: {
                    "success": data.get("success", False),
                    "error": data.get("error"),
                    "functional_groups": data.get("functional_groups", {}).get("functional_groups", {}),
                    "ring_profile": data.get("ring_profile", {}),
                }
                for smi, data in precursor_analysis.items()
            },
            "precursor_compatibility_issues": precursor_compatibility_issues,
            "step_molecule_analysis": {
                smi: {
                    "success": data.get("success", False),
                    "error": data.get("error"),
                    "functional_groups": data.get("functional_groups", {}).get("functional_groups", {}),
                    "ring_profile": data.get("ring_profile", {}),
                }
                for smi, data in step_molecule_analysis.items()
            },
            "ring_topology": ring_topology_diag,
        }

        verdict_counts[verdict] += 1
        audited.append(route_copy)
        if verdict == "FAIL":
            rejected.append(route_copy)
        else:
            recommended.append(route_copy)

    summary = {
        "strict_mode": True,
        "total_routes": len(audited),
        "recommended_routes": len(recommended),
        "rejected_routes": len(rejected),
        "verdict_counts": dict(verdict_counts),
        "analysis_required_for_pass": True,
        "target_analysis_success": bool(target_analysis and target_analysis.get("success")),
        "best_quality_score": max([float(r.get("quality_score", 0.0) or 0.0) for r in audited], default=0.0),
        "style_counts": portfolio_counts,
        "portfolio_style_requirement": (
            f"at least {MIN_REQUIRED_ROUTES} routes with mixed convergent/linear styles"
            if enforce_style_mix
            else "not_enforced"
        ),
        "portfolio_style_issue": portfolio_issue,
    }
    return audited, recommended, rejected, summary

def audit_routes_with_gate(
    *,
    routes: List[Dict[str, Any]],
    target_smiles: str,
    strict_audit: bool,
    rxn_mapper_mod: Any,
    atom_utils_mod: Any,
    strategy: Dict[str, Any],
    analyze_mod: Any,
    features_mod: Any,
    fg_mod: Any,
    enforce_style_mix: bool,
) -> Tuple[List[Dict[str, Any]], List[Dict[str, Any]], List[Dict[str, Any]], Dict[str, Any]]:
    return _audit_routes_with_gate(
        routes=routes,
        target_smiles=target_smiles,
        strict_audit=strict_audit,
        rxn_mapper_mod=rxn_mapper_mod,
        atom_utils_mod=atom_utils_mod,
        strategy=strategy,
        analyze_mod=analyze_mod,
        features_mod=features_mod,
        fg_mod=fg_mod,
        enforce_style_mix=enforce_style_mix,
    )

