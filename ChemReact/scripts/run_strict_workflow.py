"""Single-entry strict retrosynthesis workflow runner.

Usage:
  python skills/ChemReact/scripts/run_strict_workflow.py --target-smiles "<SMILES>" [--output-dir ...]
"""

import argparse
import json
import os
import subprocess
import sys
from pathlib import Path


INTERMEDIATE_MODES = {
    "planner_request_only",
    "planner_generation_pending",
    "planner_repair_requested",
    "planner_repair_pending",
}


def _default_host_command() -> str:
    return os.environ.get("CHEMREACT_HOST_PLANNER_COMMAND", "").strip()


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


def _load_json(path: Path):
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)


def _has_route_images(images_dir: Path) -> bool:
    return any(images_dir.glob("route_*_overview.png")) and any(images_dir.glob("route_*_tree.png"))


def main() -> int:
    parser = argparse.ArgumentParser(description="Run strict ChemReact workflow from target SMILES.")
    parser.add_argument("--target-smiles", required=True, help="Target molecule SMILES.")
    parser.add_argument("--output-dir", default="outputs/retro_strict_run", help="Output directory.")
    parser.add_argument(
        "--host-planner-command",
        default="",
        help="Optional host planner command. If omitted, use CHEMREACT_HOST_PLANNER_COMMAND.",
    )
    args = parser.parse_args()

    if not args.target_smiles.strip():
        raise ValueError("--target-smiles must be non-empty")

    output_dir = Path(args.output_dir).resolve()
    host_cmd = (args.host_planner_command or "").strip() or _default_host_command()
    if not host_cmd:
        print(
            "Missing host planner command. Provide a concrete configured host_planner_command via CHEMREACT_HOST_PLANNER_COMMAND or --host-planner-command.",
            file=sys.stderr,
        )
        return 6
    if _is_non_production_host_command(host_cmd):
        print(
            "host planner command is placeholder/fallback. "
            "Provide a real configured host_planner_command from your host software integration.",
            file=sys.stderr,
        )
        return 7

    run_script = Path(__file__).resolve().parent / "run_closed_loop.py"
    cmd = [
        sys.executable,
        str(run_script),
        "--target-smiles",
        args.target_smiles.strip(),
        "--output-dir",
        str(output_dir),
        "--auto-propose-routes",
        "--planner-backend",
        "host",
        "--auto-propose-max-routes",
        "5",
        "--strict-audit",
        "--planner-auto-repair",
        "--planner-max-iterations",
        "5",
        "--host-planner-command",
        host_cmd,
    ]

    proc = subprocess.run(cmd, capture_output=True, text=True, check=False)
    if proc.stdout:
        print(proc.stdout)
    if proc.returncode != 0:
        if proc.stderr:
            print(proc.stderr, file=sys.stderr)
        return proc.returncode

    summary_path = output_dir / "run_summary.json"
    if not summary_path.exists():
        print(f"Missing run summary: {summary_path}", file=sys.stderr)
        return 2
    summary = _load_json(summary_path)

    mode = str(summary.get("mode", "")).strip()
    if mode in INTERMEDIATE_MODES:
        print(f"Workflow stopped in intermediate mode: {mode}", file=sys.stderr)
        return 3

    audit = summary.get("audit", {})
    recommended_routes = 0
    if isinstance(audit, dict):
        recommended_routes = int(audit.get("recommended_routes", 0) or 0)
    loop = summary.get("planner_loop", {})
    reached_repair_limit = bool(isinstance(loop, dict) and loop.get("reached_repair_limit"))
    if recommended_routes <= 0 and not reached_repair_limit:
        print("Strict workflow finished with zero recommended routes.", file=sys.stderr)
        return 8

    report_path = output_dir / "RETRO_REPORT.md"
    images_dir = output_dir / "images"
    if not report_path.exists():
        print(f"Missing report: {report_path}", file=sys.stderr)
        return 4
    if not images_dir.exists() or not _has_route_images(images_dir):
        print(f"Missing route overview/tree images in: {images_dir}", file=sys.stderr)
        return 5

    print(
        json.dumps(
            {
                "success": True,
                "output_dir": str(output_dir),
                "report": str(report_path),
                "run_summary": str(summary_path),
                "recommended_routes": recommended_routes,
                "reached_repair_limit": reached_repair_limit,
            },
            ensure_ascii=False,
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
