"""
RouteManager — multi-route management and comparison for retrosynthesis.

Manages creation, tracking, abandonment, and comparison of multiple
synthesis routes for the same target molecule.  Persists state to
``{session_dir}/routes_state.json``.

Migrated from ``core/taskflow/route_manager.py``.
Key changes:

* Imports ``RetroRoute``, ``RetroTaskList`` from ``tools.models.workflow_models``
  (the new canonical data models).
* Imports ``TaskStatus``, ``TaskType`` from ``tools.common.status``.
* Route-level fields that were previously top-level attributes on the old
  ``RetroRoute`` (``route_name``, ``strategy``, ``status``, ``retro_graph``,
  ``score``, ``images``) are now stored inside ``RetroRoute.metadata``.
"""

from __future__ import annotations

import json
import logging
import uuid
from pathlib import Path
from typing import Any, Dict, List, Optional

from tools.common.status import RouteStatus, TaskStatus, TaskType
from tools.models.workflow_models import RetroRoute, RetroTaskList

__all__ = ["RouteManager"]

logger = logging.getLogger(__name__)


class RouteManager:
    """Manage multiple synthesis routes for a single target molecule."""

    def __init__(self, target_smiles: str, session_dir: Path) -> None:
        self.target_smiles = target_smiles
        self.session_dir = Path(session_dir)
        self.routes: Dict[str, RetroRoute] = {}

    # ------------------------------------------------------------------
    # Route lifecycle
    # ------------------------------------------------------------------

    def create_route(self, strategy: str, route_name: str) -> RetroRoute:
        """Create a new route with a unique id and its own subdirectory.

        The subdirectory ``{session_dir}/{route_name}/`` and an ``images/``
        child directory are created on disk.
        """
        route_id = uuid.uuid4().hex[:8]

        route = RetroRoute(
            route_id=route_id,
            target_smiles=self.target_smiles,
            metadata={
                "route_name": route_name,
                "strategy": strategy,
                "status": RouteStatus.PLANNING.value,
                "retro_graph": None,
                "score": None,
                "images": {},
            },
        )

        # Create route directory and images/ subdirectory
        route_dir = self.session_dir / route_name
        route_dir.mkdir(parents=True, exist_ok=True)
        (route_dir / "images").mkdir(parents=True, exist_ok=True)

        self.routes[route_id] = route
        return route

    # ------------------------------------------------------------------
    # Queries
    # ------------------------------------------------------------------

    def get_route(self, route_id: str) -> Optional[RetroRoute]:
        """Return the route with the given *route_id*, or ``None``."""
        return self.routes.get(route_id)

    def list_routes(self) -> List[RetroRoute]:
        """Return all routes (active and abandoned)."""
        return list(self.routes.values())

    def get_active_routes(self) -> List[RetroRoute]:
        """Return all routes whose status is **not** ``abandoned``."""
        return [
            r for r in self.routes.values()
            if r.metadata.get("status") != RouteStatus.ABANDONED.value
        ]

    # ------------------------------------------------------------------
    # Abandon
    # ------------------------------------------------------------------

    def abandon_route(self, route_id: str, reason: str) -> None:
        """Mark a route as abandoned and persist the reason.

        The reason is stored in the route's ``metadata`` under the key
        ``abandon_reason``.
        """
        route = self.routes.get(route_id)
        if route is None:
            return
        route.metadata["status"] = RouteStatus.ABANDONED.value
        route.metadata["abandon_reason"] = reason

    # ------------------------------------------------------------------
    # Comparison
    # ------------------------------------------------------------------

    def compare_routes(self) -> Dict[str, Any]:
        """Compare all completed routes on key metrics.

        Returns a dict keyed by *route_id* with sub-dicts containing:
        ``total_steps``, ``terminal_precursors``, ``avg_sa_score``,
        ``max_depth``, ``validation_pass_rate``.

        If no completed routes exist, returns an empty dict.
        """
        completed = [
            r for r in self.routes.values()
            if r.metadata.get("status") == RouteStatus.COMPLETED.value
        ]
        if not completed:
            return {}

        comparison: Dict[str, Any] = {}
        for route in completed:
            tasks = route.tasks.tasks
            total_steps = len(tasks)

            # Terminal precursors: tasks marked SKIPPED (SA < threshold)
            terminal_precursors = sum(
                1 for t in tasks if t.status == TaskStatus.SKIPPED
            )

            # Average SA_Score from task results
            sa_scores: List[float] = []
            for t in tasks:
                if t.result and isinstance(t.result, dict):
                    sa = t.result.get("sa_score")
                    if sa is not None:
                        try:
                            sa_scores.append(float(sa))
                        except (TypeError, ValueError):
                            pass
            avg_sa = sum(sa_scores) / len(sa_scores) if sa_scores else 0.0

            # Max depth
            max_depth = max((t.depth for t in tasks), default=0)

            # Validation pass rate
            validate_tasks = [
                t for t in tasks if t.task_type == TaskType.VALIDATE
            ]
            if validate_tasks:
                passed = sum(
                    1
                    for t in validate_tasks
                    if t.result.get("verdict") == "PASS"
                )
                pass_rate = passed / len(validate_tasks)
            else:
                pass_rate = 0.0

            route_name = route.metadata.get("route_name", route.route_id)
            strategy = route.metadata.get("strategy", "")

            comparison[route.route_id] = {
                "route_name": route_name,
                "strategy": strategy,
                "total_steps": total_steps,
                "terminal_precursors": terminal_precursors,
                "avg_sa_score": round(avg_sa, 2),
                "max_depth": max_depth,
                "validation_pass_rate": round(pass_rate * 100, 1),
                "status": route.metadata.get("status", ""),
            }
        return comparison

    # ------------------------------------------------------------------
    # Comparison report
    # ------------------------------------------------------------------

    def generate_comparison_report(self) -> str:
        """Generate a ``comparison.md`` Markdown report.

        The report includes a summary table and per-route synthesis tree
        image references.  It is also written to
        ``{session_dir}/comparison.md``.
        """
        comparison = self.compare_routes()

        lines: List[str] = [
            "# 路线对比报告",
            "",
            f"**Target**: `{self.target_smiles}`",
            "",
            "## 路线概览",
            "",
            "| 路线 | 策略 | 总步数 | 终端前体 | 平均 SA_Score | 最大深度 | 验证通过率 | 状态 |",
            "|------|------|--------|----------|---------------|----------|------------|------|",
        ]

        for info in comparison.values():
            lines.append(
                f"| {info['route_name']} "
                f"| {info['strategy']} "
                f"| {info['total_steps']} "
                f"| {info['terminal_precursors']} "
                f"| {info['avg_sa_score']} "
                f"| {info['max_depth']} "
                f"| {info['validation_pass_rate']}% "
                f"| {info['status']} |"
            )

        # If no completed routes, add a note
        if not comparison:
            lines.append("| (无已完成路线) | - | - | - | - | - | - | - |")

        lines.append("")
        lines.append("## 合成树对比")
        lines.append("")

        for route in self.routes.values():
            if route.metadata.get("status") == RouteStatus.COMPLETED.value:
                route_name = route.metadata.get("route_name", route.route_id)
                strategy = route.metadata.get("strategy", "")
                lines.append(f"### {route_name}: {strategy}_route")
                lines.append(
                    f"![{route_name} 合成树]"
                    f"({route_name}/images/synthesis_tree.png)"
                )
                lines.append("")

        report = "\n".join(lines)

        # Write to disk
        report_path = self.session_dir / "comparison.md"
        report_path.parent.mkdir(parents=True, exist_ok=True)
        report_path.write_text(report, encoding="utf-8")

        return report

    # ------------------------------------------------------------------
    # Persistence
    # ------------------------------------------------------------------

    def save_state(self) -> None:
        """Persist all route states to ``routes_state.json``."""
        state = {
            "target_smiles": self.target_smiles,
            "routes": {
                rid: route.to_dict() for rid, route in self.routes.items()
            },
        }
        state_path = self.session_dir / "routes_state.json"
        state_path.parent.mkdir(parents=True, exist_ok=True)
        state_path.write_text(
            json.dumps(state, ensure_ascii=False, indent=2),
            encoding="utf-8",
        )

    def load_state(self) -> None:
        """Load route states from ``routes_state.json``.

        If the state file does not exist, this is a no-op (routes dict
        remains unchanged).
        """
        state_path = self.session_dir / "routes_state.json"
        if not state_path.exists():
            return
        data = json.loads(state_path.read_text(encoding="utf-8"))
        self.target_smiles = data.get("target_smiles", self.target_smiles)
        self.routes = {
            rid: RetroRoute.from_dict(rd)
            for rid, rd in data.get("routes", {}).items()
        }
