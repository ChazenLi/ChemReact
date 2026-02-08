"""Host-bridge planner for ChemReact strict closed-loop workflow.

This script is provider-agnostic. It does not call any model vendor API directly.
It delegates generation to a host-supplied command via:
- --bridge-command
- or CHEMREACT_LLM_BRIDGE_COMMAND

Bridge command placeholders:
- {prompt_file}
- {request_file}
- {stage}
- {round}

Bridge command must print JSON to stdout. Accepted payload:
- JSON array of routes
- or object with key `routes` as array
"""

from __future__ import annotations

import argparse
import json
import os
import re
import subprocess
import sys
from pathlib import Path
from typing import Any, Dict, List


def _read_json(path: str) -> Dict[str, Any]:
    with open(path, "r", encoding="utf-8") as f:
        data = json.load(f)
    if not isinstance(data, dict):
        raise ValueError("request-file must contain a JSON object")
    return data


def _write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2, ensure_ascii=False)


def _extract_json(text: str) -> Any:
    cleaned = (text or "").strip()
    if not cleaned:
        raise ValueError("Bridge command returned empty output")
    if cleaned.startswith("```"):
        cleaned = re.sub(r"^```(?:json)?\\s*", "", cleaned)
        cleaned = re.sub(r"\\s*```$", "", cleaned)
    try:
        return json.loads(cleaned)
    except Exception:
        pass

    decoder = json.JSONDecoder()
    for marker in ("[", "{"):
        idx = cleaned.find(marker)
        while idx >= 0:
            try:
                obj, end = decoder.raw_decode(cleaned[idx:])
                if end > 0:
                    return obj
            except Exception:
                pass
            idx = cleaned.find(marker, idx + 1)
    raise ValueError("Could not parse JSON from bridge command output")


def _normalize_routes(raw: Any) -> List[Dict[str, Any]]:
    payload = raw
    if isinstance(raw, dict) and isinstance(raw.get("routes"), list):
        payload = raw["routes"]
    if not isinstance(payload, list):
        raise ValueError("Planner output must be a JSON array (or object with routes array)")

    routes: List[Dict[str, Any]] = []
    for idx, item in enumerate(payload, start=1):
        if not isinstance(item, dict):
            continue
        steps_raw = item.get("steps", [])
        if not isinstance(steps_raw, list) or not steps_raw:
            continue

        steps: List[Dict[str, Any]] = []
        for step in steps_raw:
            if not isinstance(step, dict):
                continue
            rxn = step.get("reaction_smiles") or step.get("rxn_smiles") or step.get("smirks")
            if not isinstance(rxn, str) or ">>" not in rxn:
                continue
            reagents = step.get("reagents", [])
            if not isinstance(reagents, list):
                reagents = []
            steps.append(
                {
                    "reaction_class": str(step.get("reaction_class", "Unspecified Step")),
                    "reagents": [str(x) for x in reagents if isinstance(x, str)],
                    "conditions": str(step.get("conditions", "TBD")),
                    "reaction_smiles": rxn,
                }
            )
        if not steps:
            continue

        route_id = item.get("route_id") if isinstance(item.get("route_id"), int) else idx
        style = str(item.get("synthesis_style", "linear")).lower()
        if style not in {"linear", "convergent"}:
            style = "linear" if idx % 2 == 1 else "convergent"
        score = item.get("score", 6.0)
        if not isinstance(score, (int, float)):
            score = 6.0
        verdict = str(item.get("audit_verdict", "CONDITIONAL")).upper()
        if verdict not in {"PASS", "FAIL", "CONDITIONAL"}:
            verdict = "CONDITIONAL"
        issues = item.get("critical_issues", [])
        if not isinstance(issues, list):
            issues = []
        precursors = item.get("precursors", [])
        if not isinstance(precursors, list):
            precursors = []

        routes.append(
            {
                "route_id": route_id,
                "synthesis_style": style,
                "score": float(score),
                "audit_verdict": verdict,
                "critical_issues": [str(x) for x in issues if isinstance(x, str)],
                "precursors": [str(x) for x in precursors if isinstance(x, str)],
                "steps": steps,
            }
        )

    if not routes:
        raise ValueError("Planner produced no valid routes after normalization")
    return routes


def _run_bridge_command(command_template: str, prompt_file: str, request_file: str, stage: str, round_idx: int) -> str:
    required = ("{prompt_file}", "{request_file}", "{stage}", "{round}")
    missing = [k for k in required if k not in command_template]
    if missing:
        raise ValueError(
            "bridge command missing placeholders: " + ", ".join(required)
        )
    cmd = command_template.format(
        prompt_file=prompt_file,
        request_file=request_file,
        stage=stage,
        round=round_idx,
    )
    proc = subprocess.run(cmd, shell=True, capture_output=True, text=True, check=False)
    if proc.returncode != 0:
        raise RuntimeError(
            f"bridge command failed (code={proc.returncode}): {proc.stderr.strip()}"
        )
    return proc.stdout


def main() -> int:
    parser = argparse.ArgumentParser(description="ChemReact host-bridge planner runner.")
    parser.add_argument("--prompt-file", required=True)
    parser.add_argument("--output-file", required=True)
    parser.add_argument("--request-file", required=True)
    parser.add_argument("--stage", required=True, choices=["plan", "repair"])
    parser.add_argument("--round", type=int, required=True)
    parser.add_argument("--bridge-command", default="")
    args = parser.parse_args()

    try:
        _ = _read_json(args.request_file)
        bridge_cmd = (args.bridge_command or "").strip() or os.environ.get("CHEMREACT_LLM_BRIDGE_COMMAND", "").strip()
        if not bridge_cmd:
            raise RuntimeError(
                "Missing bridge command. Set --bridge-command or CHEMREACT_LLM_BRIDGE_COMMAND."
            )

        raw_text = _run_bridge_command(
            bridge_cmd,
            prompt_file=str(Path(args.prompt_file).resolve()),
            request_file=str(Path(args.request_file).resolve()),
            stage=args.stage,
            round_idx=args.round,
        )
        routes = _normalize_routes(_extract_json(raw_text))

        out_path = Path(args.output_file).resolve()
        _write_json(out_path, routes)
        return 0
    except Exception as exc:
        print(f"llm_planner failed: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
