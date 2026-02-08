"""Built-in fallback host planner runner for ChemReact closed-loop execution.

This runner exists to keep the adapter contract executable even when an external
LLM host runner is not configured yet.
"""

import argparse
import json
from pathlib import Path
from typing import Any, Dict, List


def _read_json(path: str) -> Dict[str, Any]:
    with open(path, "r", encoding="utf-8") as f:
        data = json.load(f)
    if not isinstance(data, dict):
        raise ValueError("request-file must contain a JSON object")
    return data


def _write_routes(path: str, routes: List[Dict[str, Any]]) -> None:
    out = Path(path)
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w", encoding="utf-8") as f:
        json.dump(routes, f, indent=2, ensure_ascii=False)


def _identity_routes(target: str) -> List[Dict[str, Any]]:
    """Generate a minimal 3-route portfolio so the pipeline can continue."""
    target = (target or "").strip() or "CCO"
    return [
        {
            "route_id": 1,
            "synthesis_style": "linear",
            "score": 4.0,
            "audit_verdict": "CONDITIONAL",
            "critical_issues": ["Fallback planner route: replace with real host LLM runner for production."],
            "precursors": [target],
            "steps": [
                {
                    "reaction_class": "Fallback Identity Step",
                    "reagents": [],
                    "conditions": "N/A",
                    "reaction_smiles": f"{target}>>{target}",
                }
            ],
        },
        {
            "route_id": 2,
            "synthesis_style": "convergent",
            "score": 3.8,
            "audit_verdict": "CONDITIONAL",
            "critical_issues": ["Fallback planner route: replace with real host LLM runner for production."],
            "precursors": [target],
            "steps": [
                {
                    "reaction_class": "Fallback Identity Step",
                    "reagents": [],
                    "conditions": "N/A",
                    "reaction_smiles": f"{target}>>{target}",
                }
            ],
        },
        {
            "route_id": 3,
            "synthesis_style": "convergent",
            "score": 3.6,
            "audit_verdict": "CONDITIONAL",
            "critical_issues": ["Fallback planner route: replace with real host LLM runner for production."],
            "precursors": [target],
            "steps": [
                {
                    "reaction_class": "Fallback Identity Step",
                    "reagents": [],
                    "conditions": "N/A",
                    "reaction_smiles": f"{target}>>{target}",
                }
            ],
        },
    ]


def main() -> int:
    parser = argparse.ArgumentParser(description="ChemReact fallback host planner runner.")
    parser.add_argument("--prompt-file", required=True)
    parser.add_argument("--output-file", required=True)
    parser.add_argument("--request-file", required=True)
    parser.add_argument("--stage", required=True, choices=["plan", "repair"])
    parser.add_argument("--round", type=int, required=True)
    args = parser.parse_args()

    payload = _read_json(args.request_file)

    # Prefer route candidates from repair request when available.
    if args.stage == "repair":
        candidate_routes = payload.get("candidate_routes")
        if isinstance(candidate_routes, list) and candidate_routes:
            _write_routes(args.output_file, candidate_routes)
            return 0
        rejected_routes = payload.get("rejected_routes")
        if isinstance(rejected_routes, list) and rejected_routes:
            _write_routes(args.output_file, rejected_routes)
            return 0

    target = str(payload.get("target_smiles", "")).strip()
    _write_routes(args.output_file, _identity_routes(target))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
