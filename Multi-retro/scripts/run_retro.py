#!/usr/bin/env python
"""
MultiRetro 统一 CLI 入口
========================

Replaces all legacy scripts (run_skill, run_taskflow, run_mesh,
run_closed_loop, retro_batch, etc.) with a single entry point.

Subcommands
-----------
  skill     — Call a single skill by name with JSON args
  plan      — Create a retrosynthesis plan for a target SMILES
  run       — Execute tasks in a session
  decide    — Submit a decision for a paused session
  explore   — Exploration tool calls at decision points
  compare   — Multi-route comparison
  finalize  — Finalize and generate reports
  status    — View session status
  batch     — Batch execution from a JSON task file

Usage::

    python scripts/run_retro.py skill analyze_molecule --args '{"smiles":"CCO"}'
    python scripts/run_retro.py plan --target_smiles "c1ccccc1"
    python scripts/run_retro.py run --session_dir outputs/session1 --route route_main
    python scripts/run_retro.py batch tasks.json
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any, Dict

# ---------------------------------------------------------------------------
# Path setup — ensure project root is importable
# ---------------------------------------------------------------------------
_HERE = Path(__file__).resolve().parent
_ROOT = _HERE.parent
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))


# ---------------------------------------------------------------------------
# Output helper
# ---------------------------------------------------------------------------

def _json_output(data: Dict[str, Any]) -> None:
    """Write *data* as UTF-8 JSON to stdout (Windows-safe)."""
    output = json.dumps(data, indent=2, ensure_ascii=False, default=str)
    sys.stdout.buffer.write(output.encode("utf-8", errors="replace"))
    sys.stdout.buffer.write(b"\n")
    sys.stdout.buffer.flush()


# ---------------------------------------------------------------------------
# Skill registry (lazy — only built when 'skill' subcommand is used)
# ---------------------------------------------------------------------------

def _get_skill_registry() -> Dict[str, Dict[str, Any]]:
    """Return the skill name → {skill_cls, description} mapping."""
    from tools.skills import (
        AnalyzeMoleculeSkill,
        AnalyzeScaffoldSkill,
        AnalyzeSelectivitySkill,
        BreakBondSkill,
        CheckAvailabilitySkill,
        DecideStrategySkill,
        GetGlobalStrategySkill,
        MapAtomsSkill,
        PlanRingStrategySkill,
        ProposeDisconnectionSkill,
        RenderReportSkill,
        RepairReactionSkill,
        ValidateReactionSkill,
    )
    return {
        "analyze_molecule": {
            "skill_cls": AnalyzeMoleculeSkill,
            "description": "Analyze target molecule structure and properties",
        },
        "analyze_scaffold": {
            "skill_cls": AnalyzeScaffoldSkill,
            "description": "Analyze molecular scaffold for strategic disconnections",
        },
        "analyze_selectivity": {
            "skill_cls": AnalyzeSelectivitySkill,
            "description": "Analyze chemoselectivity and protection requirements",
        },
        "break_bond": {
            "skill_cls": BreakBondSkill,
            "description": "Break bond by atom indices to generate precursors",
        },
        "check_availability": {
            "skill_cls": CheckAvailabilitySkill,
            "description": "Check precursor availability via SA score",
        },
        "decide_strategy": {
            "skill_cls": DecideStrategySkill,
            "description": "Comprehensive strategy decision (linear/convergent)",
        },
        "get_global_strategy": {
            "skill_cls": GetGlobalStrategySkill,
            "description": "Get global strategy context for LLM prompts",
        },
        "map_atoms": {
            "skill_cls": MapAtomsSkill,
            "description": "Add atom-to-atom mapping using RXNMapper",
        },
        "plan_ring_strategy": {
            "skill_cls": PlanRingStrategySkill,
            "description": "Plan ring formation strategy",
        },
        "propose_disconnection": {
            "skill_cls": ProposeDisconnectionSkill,
            "description": "Propose retrosynthetic disconnection points",
        },
        "render_report": {
            "skill_cls": RenderReportSkill,
            "description": "Generate synthesis report from RetroGraph",
        },
        "repair_reaction": {
            "skill_cls": RepairReactionSkill,
            "description": "Repair invalid reaction SMILES",
        },
        "validate_reaction": {
            "skill_cls": ValidateReactionSkill,
            "description": "Validate reaction for atom balance and plausibility",
        },
    }


# ---------------------------------------------------------------------------
# Subcommand handlers
# ---------------------------------------------------------------------------

def _cmd_skill(args: argparse.Namespace) -> int:
    """Call a single skill by name with JSON args."""
    registry = _get_skill_registry()
    skill_name = args.skill_name

    # List available skills
    if skill_name == "list":
        skills = {k: v["description"] for k, v in registry.items()}
        _json_output({"success": True, "skills": skills})
        return 0

    config = registry.get(skill_name)
    if config is None:
        _json_output({"success": False, "error": f"Unknown skill: {skill_name}"})
        return 1

    # Parse JSON args
    try:
        skill_args = json.loads(args.args) if args.args else {}
    except json.JSONDecodeError as e:
        _json_output({"success": False, "error": f"Invalid JSON args: {e}"})
        return 1

    skill_cls = config["skill_cls"]
    try:
        result = skill_cls().execute(skill_args)
        _json_output(result)
        return 0 if result.get("success", False) else 1
    except Exception as exc:
        _json_output({"success": False, "error": str(exc)})
        return 1


def _cmd_plan(args: argparse.Namespace) -> int:
    """Create a retrosynthesis plan for a target SMILES."""
    from orchestrator.orchestrator import RetroOrchestrator

    session_dir = Path(args.output_dir)
    orch = RetroOrchestrator(session_dir=session_dir)
    result = orch.plan(target_smiles=args.target_smiles)
    _json_output(result)
    return 0 if result.get("success", False) else 1


def _cmd_run(args: argparse.Namespace) -> int:
    """Execute tasks in a session."""
    from orchestrator.orchestrator import RetroOrchestrator

    session_dir = Path(args.session_dir)
    if not session_dir.exists():
        _json_output({"success": False, "error": f"Session directory not found: {session_dir}"})
        return 1

    orch = RetroOrchestrator(session_dir=session_dir)
    orch.load_session()
    result = orch.run_all(route_id=args.route)
    _json_output(result)
    return 0 if result.get("success", False) else 1


def _cmd_decide(args: argparse.Namespace) -> int:
    """Submit a decision for a paused session."""
    from orchestrator.orchestrator import RetroOrchestrator
    from tools.models.workflow_models import DecisionInstruction

    session_dir = Path(args.session_dir)
    if not session_dir.exists():
        _json_output({"success": False, "error": f"Session directory not found: {session_dir}"})
        return 1

    # Parse decision JSON
    decision_str = args.decision
    if decision_str.startswith("@"):
        fpath = Path(decision_str[1:])
        if not fpath.exists():
            _json_output({"success": False, "error": f"Decision file not found: {fpath}"})
            return 1
        decision_str = fpath.read_text(encoding="utf-8")

    try:
        decision_data = json.loads(decision_str)
    except json.JSONDecodeError as e:
        _json_output({"success": False, "error": f"Invalid decision JSON: {e}"})
        return 1

    try:
        instruction = DecisionInstruction.from_dict(decision_data)
    except (ValueError, KeyError) as e:
        _json_output({"success": False, "error": f"Invalid decision format: {e}"})
        return 1

    orch = RetroOrchestrator(session_dir=session_dir)
    orch.load_session()
    result = orch.decide(route_id=args.route, instruction=instruction)
    _json_output(result)
    return 0 if result.get("success", False) else 1


def _cmd_explore(args: argparse.Namespace) -> int:
    """Exploration tool calls at decision points."""
    from orchestrator.orchestrator import RetroOrchestrator

    session_dir = Path(args.session_dir)
    if not session_dir.exists():
        _json_output({"success": False, "error": f"Session directory not found: {session_dir}"})
        return 1

    # Parse params
    try:
        params = json.loads(args.params) if args.params else {}
    except json.JSONDecodeError as e:
        _json_output({"success": False, "error": f"Invalid params JSON: {e}"})
        return 1

    orch = RetroOrchestrator(session_dir=session_dir)
    orch.load_session()
    result = orch.explore(route_id=args.route, tool_id=args.tool, params=params)
    _json_output(result)
    return 0 if result.get("success", False) else 1


def _cmd_compare(args: argparse.Namespace) -> int:
    """Multi-route comparison."""
    from orchestrator.orchestrator import RetroOrchestrator

    session_dir = Path(args.session_dir)
    if not session_dir.exists():
        _json_output({"success": False, "error": f"Session directory not found: {session_dir}"})
        return 1

    orch = RetroOrchestrator(session_dir=session_dir)
    orch.load_session()
    result = orch.compare()
    _json_output(result)
    return 0 if result.get("success", False) else 1


def _cmd_finalize(args: argparse.Namespace) -> int:
    """Finalize and generate reports."""
    from orchestrator.orchestrator import RetroOrchestrator

    session_dir = Path(args.session_dir)
    if not session_dir.exists():
        _json_output({"success": False, "error": f"Session directory not found: {session_dir}"})
        return 1

    orch = RetroOrchestrator(session_dir=session_dir)
    orch.load_session()
    result = orch.finalize(route_id=args.route)
    _json_output(result)
    return 0 if result.get("success", False) else 1


def _cmd_status(args: argparse.Namespace) -> int:
    """View session status."""
    from orchestrator.orchestrator import RetroOrchestrator

    session_dir = Path(args.session_dir)
    if not session_dir.exists():
        _json_output({"success": False, "error": f"Session directory not found: {session_dir}"})
        return 1

    orch = RetroOrchestrator(session_dir=session_dir)
    orch.load_session()
    result = orch.get_status()
    _json_output(result)
    return 0 if result.get("success", False) else 1


def _cmd_batch(args: argparse.Namespace) -> int:
    """Batch execution from a JSON task file."""
    task_file = Path(args.task_file)
    if not task_file.exists():
        _json_output({"success": False, "error": f"Task file not found: {task_file}"})
        return 1

    try:
        data = json.loads(task_file.read_text(encoding="utf-8"))
    except json.JSONDecodeError as e:
        _json_output({"success": False, "error": f"Invalid JSON in task file: {e}"})
        return 1

    # Accept both {"tasks": [...]} and bare [...]
    tasks = data.get("tasks", data) if isinstance(data, dict) else data
    if not isinstance(tasks, list):
        _json_output({"success": False, "error": "Task file must contain {\"tasks\": [...]} or a JSON array"})
        return 1

    registry = _get_skill_registry()
    results = []
    has_failure = False

    for i, task_spec in enumerate(tasks, 1):
        skill_name = task_spec.get("skill", "")
        task_args = task_spec.get("args", {})

        config = registry.get(skill_name)
        if config is None:
            result = {"success": False, "error": f"Unknown skill: {skill_name}"}
            has_failure = True
        else:
            try:
                result = config["skill_cls"]().execute(task_args)
                if not result.get("success", False):
                    has_failure = True
            except Exception as exc:
                result = {"success": False, "error": str(exc)}
                has_failure = True

        results.append({"task_index": i, "skill": skill_name, "result": result})

    _json_output({
        "success": not has_failure,
        "total": len(tasks),
        "results": results,
    })
    return 1 if has_failure else 0


# ---------------------------------------------------------------------------
# CLI setup
# ---------------------------------------------------------------------------

def build_parser() -> argparse.ArgumentParser:
    """Build and return the argument parser (exposed for testing)."""
    parser = argparse.ArgumentParser(
        prog="run_retro",
        description="MultiRetro 逆合成分析 — 统一 CLI 入口",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    subparsers = parser.add_subparsers(dest="command", help="Available commands")

    # -- skill --
    skill_p = subparsers.add_parser(
        "skill",
        help="Call a single skill by name (use 'list' to see available skills)",
    )
    skill_p.add_argument("skill_name", help="Skill name or 'list'")
    skill_p.add_argument("--args", default=None, help="JSON string of skill arguments")

    # -- plan --
    plan_p = subparsers.add_parser("plan", help="Create a retrosynthesis plan")
    plan_p.add_argument("--target_smiles", required=True, help="Target molecule SMILES")
    plan_p.add_argument("--output_dir", default="outputs/retro", help="Output directory")

    # -- run --
    run_p = subparsers.add_parser("run", help="Execute tasks in a session")
    run_p.add_argument("--session_dir", required=True, help="Session directory path")
    run_p.add_argument("--route", required=True, help="Route ID")

    # -- decide --
    decide_p = subparsers.add_parser("decide", help="Submit a decision for a paused session")
    decide_p.add_argument("--session_dir", required=True, help="Session directory path")
    decide_p.add_argument("--route", required=True, help="Route ID")
    decide_p.add_argument(
        "--decision", required=True,
        help="Decision JSON string, or @filepath to load from file",
    )

    # -- explore --
    explore_p = subparsers.add_parser("explore", help="Exploration tool calls at decision points")
    explore_p.add_argument("--session_dir", required=True, help="Session directory path")
    explore_p.add_argument("--route", required=True, help="Route ID")
    explore_p.add_argument("--tool", required=True, help="Exploration tool ID")
    explore_p.add_argument("--params", default=None, help="Tool parameters as JSON string")

    # -- compare --
    compare_p = subparsers.add_parser("compare", help="Multi-route comparison")
    compare_p.add_argument("--session_dir", required=True, help="Session directory path")

    # -- finalize --
    finalize_p = subparsers.add_parser("finalize", help="Finalize and generate reports")
    finalize_p.add_argument("--session_dir", required=True, help="Session directory path")
    finalize_p.add_argument("--route", required=True, help="Route ID")

    # -- status --
    status_p = subparsers.add_parser("status", help="View session status")
    status_p.add_argument("--session_dir", required=True, help="Session directory path")

    # -- batch --
    batch_p = subparsers.add_parser("batch", help="Batch execution from a JSON task file")
    batch_p.add_argument("task_file", help="Path to JSON task file")

    return parser


def main(argv: list[str] | None = None) -> int:
    """CLI entry point.

    Parameters
    ----------
    argv:
        Argument list to parse.  Defaults to ``sys.argv[1:]``.
    """
    parser = build_parser()
    args = parser.parse_args(argv)

    if not args.command:
        parser.print_help()
        return 1

    handlers = {
        "skill": _cmd_skill,
        "plan": _cmd_plan,
        "run": _cmd_run,
        "decide": _cmd_decide,
        "explore": _cmd_explore,
        "compare": _cmd_compare,
        "finalize": _cmd_finalize,
        "status": _cmd_status,
        "batch": _cmd_batch,
    }

    handler = handlers.get(args.command)
    if handler is None:
        _json_output({"success": False, "error": f"Unknown command: {args.command}"})
        return 1

    try:
        return handler(args)
    except Exception as exc:
        _json_output({"success": False, "error": str(exc)})
        return 1


if __name__ == "__main__":
    sys.exit(main())
