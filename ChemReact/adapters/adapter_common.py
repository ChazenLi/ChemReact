import argparse
import json
import os
import subprocess
import sys
from pathlib import Path
from typing import Dict, Any


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


def _load_request(path: str) -> Dict[str, Any]:
    with open(path, "r", encoding="utf-8") as f:
        payload = json.load(f)
    if not isinstance(payload, dict):
        raise ValueError("request_file must contain a JSON object")
    return payload


def _validate_request(payload: Dict[str, Any]) -> None:
    if "target_smiles" not in payload:
        raise ValueError("request missing required field: target_smiles")
    if not isinstance(payload["target_smiles"], str):
        raise ValueError("request field target_smiles must be str")
    if payload["target_smiles"].strip() in {"", "__REPLACE_WITH_TARGET_SMILES__"}:
        raise ValueError("request field target_smiles must be a real SMILES, not placeholder")

    has_routes_file = isinstance(payload.get("routes_file"), str) and bool(payload.get("routes_file"))
    auto_propose = bool(payload.get("auto_propose_routes"))
    if not has_routes_file and not auto_propose:
        raise ValueError("request must provide routes_file or set auto_propose_routes=true")
    optional_types = {
        "strategy_file": str,
        "vis_plan_file": str,
        "output_dir": str,
        "top_k": int,
        "target_legend": str,
        "force_field": str,
        "emit_vis_plan_template": str,
        "planner_backend": str,
        "planner_routes_file": str,
        "auto_propose_routes": bool,
        "auto_propose_max_routes": int,
        "planner_request_only": bool,
        "strict_audit": bool,
        "planner_auto_repair": bool,
        "planner_max_iterations": int,
        "host_planner_command": str,
    }
    for key, typ in optional_types.items():
        if key in payload and payload[key] is not None and not isinstance(payload[key], typ):
            raise ValueError(f"request field {key} must be {typ.__name__}")
    planner_backend = payload.get("planner_backend")
    if planner_backend is not None and planner_backend != "host":
        raise ValueError("request field planner_backend must be 'host'")
    auto_max = payload.get("auto_propose_max_routes")
    if auto_max is not None and not (3 <= auto_max <= 5):
        raise ValueError("request field auto_propose_max_routes must be in [3, 5]")
    planner_iters = payload.get("planner_max_iterations")
    if planner_iters is not None and not (1 <= planner_iters <= 5):
        raise ValueError("request field planner_max_iterations must be in [1, 5]")
    cmd = payload.get("host_planner_command")
    if cmd is not None:
        required = ("{prompt_file}", "{output_file}", "{request_file}", "{stage}", "{round}")
        missing = [key for key in required if key not in cmd]
        if missing:
            raise ValueError(
                "request field host_planner_command must include placeholders: "
                + ", ".join(required)
            )
        if _is_non_production_host_command(cmd):
            raise ValueError(
                "request field host_planner_command is placeholder/fallback. "
                "Provide a real host planner command from opencode/claudecode integration."
            )


def _resolve_path(request_file: Path, value: str) -> str:
    p = Path(value)
    if p.is_absolute():
        return str(p)
    return str((request_file.parent / p).resolve())


def _build_command(run_script: Path, request_file: Path, payload: Dict[str, Any], default_output_dir: str):
    output_dir = payload.get("output_dir", default_output_dir)
    output_dir = _resolve_path(request_file, output_dir)
    cmd = [
        sys.executable,
        str(run_script),
        "--target-smiles",
        payload["target_smiles"],
        "--output-dir",
        output_dir,
    ]
    if payload.get("routes_file"):
        cmd.extend(["--routes-file", _resolve_path(request_file, payload["routes_file"])])
    if payload.get("strategy_file"):
        cmd.extend(["--strategy-file", _resolve_path(request_file, payload["strategy_file"])])
    if payload.get("vis_plan_file"):
        cmd.extend(["--vis-plan-file", _resolve_path(request_file, payload["vis_plan_file"])])
    if payload.get("top_k") is not None:
        cmd.extend(["--top-k", str(payload["top_k"])])
    if payload.get("target_legend"):
        cmd.extend(["--target-legend", payload["target_legend"]])
    if payload.get("force_field"):
        cmd.extend(["--force-field", payload["force_field"]])
    if payload.get("emit_vis_plan_template"):
        cmd.extend(["--emit-vis-plan-template", _resolve_path(request_file, payload["emit_vis_plan_template"])])
    if payload.get("auto_propose_routes"):
        cmd.append("--auto-propose-routes")
    if payload.get("planner_backend"):
        cmd.extend(["--planner-backend", payload["planner_backend"]])
    if payload.get("auto_propose_max_routes") is not None:
        cmd.extend(["--auto-propose-max-routes", str(payload["auto_propose_max_routes"])])
    if payload.get("planner_routes_file"):
        cmd.extend(["--planner-routes-file", _resolve_path(request_file, payload["planner_routes_file"])])
    if payload.get("planner_request_only"):
        cmd.append("--planner-request-only")
    if payload.get("strict_audit", True):
        cmd.append("--strict-audit")
    planner_auto_repair = payload.get("planner_auto_repair")
    if planner_auto_repair is None or planner_auto_repair is True:
        cmd.append("--planner-auto-repair")
    else:
        cmd.append("--no-planner-auto-repair")
    if payload.get("planner_max_iterations") is not None:
        cmd.extend(["--planner-max-iterations", str(payload["planner_max_iterations"])])
    if payload.get("host_planner_command"):
        cmd.extend(["--host-planner-command", payload["host_planner_command"]])
    return cmd


def run_adapter(host_name: str, default_output_dir: str):
    parser = argparse.ArgumentParser(
        description=f"{host_name} adapter for ChemReact pipeline."
    )
    parser.add_argument("--request-file", required=True, help="JSON request payload path.")
    parser.add_argument("--print-command", action="store_true", help="Only print resolved command.")
    args = parser.parse_args()

    request_file = Path(args.request_file).resolve()
    payload = _load_request(str(request_file))

    # Default to fully automated strict closed loop.
    payload.setdefault("strict_audit", True)
    payload.setdefault("planner_auto_repair", True)

    env_cmd = os.environ.get("CHEMREACT_HOST_PLANNER_COMMAND", "").strip()
    if (not payload.get("host_planner_command")) and env_cmd:
        payload["host_planner_command"] = env_cmd

    if payload.get("planner_auto_repair"):
        cmd = str(payload.get("host_planner_command", "") or "").strip()
        if not cmd:
            raise ValueError(
                "planner_auto_repair is enabled but host_planner_command is missing. "
                "Provide request.host_planner_command or set CHEMREACT_HOST_PLANNER_COMMAND."
            )

    _validate_request(payload)

    run_script = Path(__file__).resolve().parents[1] / "scripts" / "run_closed_loop.py"
    cmd = _build_command(run_script, request_file, payload, default_output_dir)

    if args.print_command:
        print(" ".join(cmd))
        return

    proc = subprocess.run(cmd, capture_output=True, text=True)
    print(proc.stdout)
    if proc.returncode != 0:
        raise RuntimeError(proc.stderr.strip() or f"{host_name} adapter execution failed")
