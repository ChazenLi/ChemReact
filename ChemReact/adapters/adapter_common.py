import argparse
import json
import subprocess
import sys
from pathlib import Path
from typing import Dict, Any


def _load_request(path: str) -> Dict[str, Any]:
    with open(path, "r", encoding="utf-8") as f:
        payload = json.load(f)
    if not isinstance(payload, dict):
        raise ValueError("request_file must contain a JSON object")
    return payload


def _validate_request(payload: Dict[str, Any]) -> None:
    required = {
        "target_smiles": str,
        "routes_file": str,
    }
    for key, typ in required.items():
        if key not in payload:
            raise ValueError(f"request missing required field: {key}")
        if not isinstance(payload[key], typ):
            raise ValueError(f"request field {key} must be {typ.__name__}")
    optional_types = {
        "strategy_file": str,
        "vis_plan_file": str,
        "output_dir": str,
        "top_k": int,
        "target_legend": str,
        "force_field": str,
        "emit_vis_plan_template": str,
    }
    for key, typ in optional_types.items():
        if key in payload and payload[key] is not None and not isinstance(payload[key], typ):
            raise ValueError(f"request field {key} must be {typ.__name__}")


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
        "--routes-file",
        _resolve_path(request_file, payload["routes_file"]),
        "--output-dir",
        output_dir,
    ]
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
