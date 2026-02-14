"""
RetroOrchestrator — flat pipeline orchestrator for retrosynthesis.

Two-phase execution: Phase 1 (decomposition) → Phase 2 (finalization).
Data flows through typed dataclasses; skills dispatched via skill_dispatch.
"""
from __future__ import annotations

import logging
from pathlib import Path
from typing import Any, Dict, List, Optional, Set

from tools.common.status import TaskStatus, TaskType
from tools.models.workflow_models import (
    DecisionContext, DecisionInstruction, RetroRoute, RetroTask, RetroTaskList,
)
from tools.output.graph_builder import build_retro_graph, collect_terminal_precursors
from tools.output.journal import ProcessJournal
from tools.output.report_generator import generate_synthesis_report
from tools.skills import (
    AnalyzeMoleculeSkill, AnalyzeScaffoldSkill, AnalyzeSelectivitySkill,
    BreakBondSkill, CheckAvailabilitySkill, DecideStrategySkill,
    GetGlobalStrategySkill, MapAtomsSkill, PlanRingStrategySkill,
    ProposeDisconnectionSkill, RenderReportSkill, RepairReactionSkill,
    ValidateReactionSkill,
)
from tools.workflow.decision_manager import DecisionManager
from tools.workflow.exploration import execute_exploration, ExplorationSession
from tools.workflow.route_manager import RouteManager
from tools.workflow.skill_dispatch import dispatch_skill, extract_result
from tools.workflow.snapshot import create_checkpoint, save_snapshot
from tools.workflow.subtask_manager import (
    append_subtasks, create_repair_subtask,
    create_validate_subtask, expand_precursor_subtasks,
)
from tools.workflow.task_state import (
    complete_task, create_task, fail_task, start_task, transition_status,
)

__all__ = ["RetroOrchestrator"]
logger = logging.getLogger(__name__)

_PHASE1_TYPES: Set[TaskType] = {
    TaskType.STRATEGY, TaskType.ANALYZE, TaskType.DISCONNECT,
    TaskType.VALIDATE, TaskType.REPAIR,
}
_PHASE2_TYPES: Set[TaskType] = {TaskType.AVAILABILITY, TaskType.REPORT}


class RetroOrchestrator:
    """Flat pipeline orchestrator — ties all tool modules together."""

    def __init__(self, session_dir: Path, config: Optional[Dict[str, Any]] = None) -> None:
        self.session_dir = Path(session_dir)
        self.config = config or {}
        self.route_manager: Optional[RouteManager] = None
        self.decision_manager: Optional[DecisionManager] = None
        self.journal: Optional[ProcessJournal] = None
        self.skill_registry: Dict[str, Any] = self._build_skill_registry()

    # -- Skill registry ----------------------------------------------------

    @staticmethod
    def _build_skill_registry() -> Dict[str, Any]:
        """Instantiate all skills keyed by TaskType value."""
        return {
            TaskType.STRATEGY.value: DecideStrategySkill(),
            TaskType.ANALYZE.value: AnalyzeMoleculeSkill(),
            TaskType.VALIDATE.value: ValidateReactionSkill(),
            TaskType.REPAIR.value: RepairReactionSkill(),
            TaskType.AVAILABILITY.value: CheckAvailabilitySkill(),
            TaskType.REPORT.value: RenderReportSkill(),
            "disconnect_break_bond": BreakBondSkill(),
            "disconnect_propose": ProposeDisconnectionSkill(),
            TaskType.DISCONNECT.value: ProposeDisconnectionSkill(),
            "analyze_scaffold": AnalyzeScaffoldSkill(),
            "analyze_selectivity": AnalyzeSelectivitySkill(),
            "get_global_strategy": GetGlobalStrategySkill(),
            "map_atoms": MapAtomsSkill(),
            "plan_ring_strategy": PlanRingStrategySkill(),
        }

    # -- load_session ------------------------------------------------------

    def load_session(self) -> None:
        """Load existing session state from disk.

        Restores the route_manager from ``routes_state.json`` so that
        CLI subcommands (run, decide, explore, compare, finalize, status)
        can operate on a previously created session.
        """
        state_path = self.session_dir / "routes_state.json"
        if not state_path.exists():
            raise RuntimeError(
                f"No routes_state.json found in {self.session_dir} — "
                "call plan() first to create a session"
            )
        import json
        data = json.loads(state_path.read_text(encoding="utf-8"))
        target_smiles = data.get("target_smiles", "")
        self.route_manager = RouteManager(target_smiles, self.session_dir)
        self.route_manager.load_state()

    # -- plan --------------------------------------------------------------

    def plan(self, target_smiles: str) -> Dict[str, Any]:
        """Analyse target molecule, create route with seed tasks."""
        self.route_manager = RouteManager(target_smiles, self.session_dir)
        route = self.route_manager.create_route(strategy="auto", route_name="route_1")

        for tid, ttype in [("1", TaskType.STRATEGY), ("2", TaskType.ANALYZE), ("3", TaskType.DISCONNECT)]:
            route.tasks.add_task(create_task(tid, ttype, target_smiles, depth=0))

        route_dir = self.session_dir / "route_1"
        self.decision_manager = DecisionManager(route, route_dir)
        self.journal = ProcessJournal(route_dir / "process.md")
        self.journal.write_header(route)
        save_snapshot(route, route_dir, task_id="initial")
        self.route_manager.save_state()
        logger.info("Plan created for %s — route %s", target_smiles, route.route_id)
        return {
            "success": True,
            "route_id": route.route_id,
            "target_smiles": target_smiles,
            "session_dir": str(self.session_dir),
            "tasks": len(route.tasks.tasks),
        }

    # -- run_all -----------------------------------------------------------

    def run_all(self, route_id: str) -> Dict[str, Any]:
        """Two-phase execution. Pauses at decision points."""
        route = self._get_route(route_id)
        route_dir = self._route_dir(route)
        self._ensure_managers(route, route_dir)
        if self._run_phase(route, route_dir, _PHASE1_TYPES):
            return {
                "success": True,
                "status": "paused",
                "route_id": route_id,
                "message": "Execution paused at decision point",
            }
        self._run_phase(route, route_dir, _PHASE2_TYPES)
        return {
            "success": True,
            "status": "completed",
            "route_id": route_id,
            "message": "All tasks completed",
        }

    # -- decide ------------------------------------------------------------

    def decide(self, route_id: str, instruction: DecisionInstruction) -> Dict[str, Any]:
        """Apply a decision and resume execution."""
        route = self._get_route(route_id)
        route_dir = self._route_dir(route)
        self._ensure_managers(route, route_dir)
        assert self.decision_manager is not None

        task = route.tasks.get_task(instruction.task_id)
        if task is None:
            raise ValueError(f"Task {instruction.task_id!r} not found")

        pending = self.decision_manager.load_pending_decision()
        ctx = pending.get("decision_context") if pending else None
        self.decision_manager.apply_decision(instruction, task, ctx)
        transition_status(task, TaskStatus.IN_PROGRESS, force=True)
        save_snapshot(route, route_dir)
        return self.run_all(route_id)

    # -- finalize ----------------------------------------------------------

    def finalize(self, route_id: str) -> Dict[str, Any]:
        """Build retro graph and synthesis report."""
        route = self._get_route(route_id)
        route_dir = self._route_dir(route)
        tasks = route.tasks.tasks
        graph = build_retro_graph(tasks, route.target_smiles)
        terminal = collect_terminal_precursors(tasks)
        route.metadata["retro_graph"] = graph.to_dict()
        route.metadata["status"] = TaskStatus.COMPLETED.value
        report_path = generate_synthesis_report(graph, route_dir)
        create_checkpoint(route, route_dir)
        if self.route_manager:
            self.route_manager.save_state()
        logger.info("Route %s finalized — %d terminal precursors", route_id, len(terminal))
        return {
            "success": True,
            "retro_graph": graph.to_dict(),
            "terminal_precursors": terminal,
            "report_path": str(report_path),
        }

    # -- explore -----------------------------------------------------------

    def explore(
        self,
        route_id: str,
        tool_id: str,
        params: Dict[str, Any],
    ) -> Dict[str, Any]:
        """Execute an exploration tool at a decision point (read-only).

        Args:
            route_id: The route currently paused at a decision point.
            tool_id: One of the exploration tool IDs.
            params: Parameters for the exploration tool.

        Returns:
            Exploration result dict with ``success`` key.
        """
        route = self._get_route(route_id)
        route_dir = self._route_dir(route)

        session = ExplorationSession(route_dir)
        if session.remaining <= 0:
            return {
                "success": False,
                "error": "Exploration budget exhausted",
                "remaining": 0,
            }

        result = execute_exploration(tool_id, params)
        session.record_call(tool_id, params, result)

        return {
            **result,
            "remaining_budget": session.remaining,
        }

    # -- compare -----------------------------------------------------------

    def compare(self) -> Dict[str, Any]:
        """Multi-route comparison.

        Returns:
            Comparison dict with per-route metrics.
        """
        if self.route_manager is None:
            raise RuntimeError("No route_manager — call plan() or load_session() first")

        comparison = self.route_manager.compare_routes()
        report = self.route_manager.generate_comparison_report()
        return {
            "success": True,
            "comparison": comparison,
            "report_path": str(self.session_dir / "comparison.md"),
        }

    # -- get_status --------------------------------------------------------

    def get_status(self) -> Dict[str, Any]:
        """View session status for all routes.

        Returns:
            Status dict with route summaries and rendered Markdown.
        """
        if self.route_manager is None:
            raise RuntimeError("No route_manager — call plan() or load_session() first")

        from tools.output.status_renderer import render_route_status

        routes_status: Dict[str, Any] = {}
        for route in self.route_manager.list_routes():
            routes_status[route.route_id] = {
                "target_smiles": route.target_smiles,
                "status": route.metadata.get("status", "planning"),
                "total_tasks": len(route.tasks.tasks),
                "completed": sum(
                    1 for t in route.tasks.tasks
                    if t.status in (TaskStatus.COMPLETED, TaskStatus.VALIDATED)
                ),
                "pending": sum(
                    1 for t in route.tasks.tasks
                    if t.status == TaskStatus.PENDING
                ),
                "failed": sum(
                    1 for t in route.tasks.tasks
                    if t.status == TaskStatus.FAILED
                ),
                "awaiting_decision": sum(
                    1 for t in route.tasks.tasks
                    if t.status == TaskStatus.AWAITING_DECISION
                ),
                "markdown": render_route_status(route),
            }

        return {
            "success": True,
            "session_dir": str(self.session_dir),
            "routes": routes_status,
        }

    # -- Internal: phase runner --------------------------------------------

    def _run_phase(self, route: RetroRoute, route_dir: Path, allowed: Set[TaskType]) -> bool:
        """Run pending tasks of *allowed* types. Returns True if paused."""
        while True:
            task = self._next_task(route, allowed)
            if task is None:
                return False
            result = self._execute_task(task, route, route_dir)
            if result is None:
                continue
            if self._check_decision(task, result, route, route_dir):
                return True
            complete_task(task)
            self._process_result(task, result, route)
            save_snapshot(route, route_dir)

    # -- Internal: task execution ------------------------------------------

    def _execute_task(self, task: RetroTask, route: RetroRoute, route_dir: Path) -> Optional[Dict[str, Any]]:
        """Dispatch task to skill. Returns result (task stays IN_PROGRESS) or None on failure."""
        start_task(task)
        try:
            result = dispatch_skill(task, self.skill_registry, fallback_target=route.target_smiles)
            extract_result(task, result)
            task.result = {**task.result, **result} if task.result else dict(result)
            self._log_step(task, result)
            return result
        except Exception as exc:
            logger.error("Task %s failed: %s", task.task_id, exc)
            fail_task(task)
            self._log_step(task, {"error": str(exc)})
            save_snapshot(route, route_dir)
            return None

    # -- Internal: decision check ------------------------------------------

    def _check_decision(self, task: RetroTask, result: Dict[str, Any],
                        route: RetroRoute, route_dir: Path) -> bool:
        """Returns True if execution should pause for a decision.

        Called while *task* is still IN_PROGRESS.
        """
        assert self.decision_manager is not None
        ctx = self.decision_manager.check_decision_point(task, result)
        if ctx is None:
            return False
        if self.config.get("auto_decide", False):
            fb = self.decision_manager.get_fallback_decision(ctx.decision_type, task, result)
            self.decision_manager.apply_decision(fb, task, ctx)
            self.decision_manager.record_decision(fb, ctx.decision_type, was_fallback=True)
            return False
        transition_status(task, TaskStatus.AWAITING_DECISION)
        self.decision_manager.save_pending_decision(ctx, task)
        save_snapshot(route, route_dir)
        logger.info("Paused at decision: %s (task %s)", ctx.decision_type, task.task_id)
        return True

    # -- Internal: result processing / subtask creation --------------------

    def _process_result(self, task: RetroTask, result: Dict[str, Any], route: RetroRoute) -> None:
        """Create subtasks based on task result."""
        new_tasks: List[RetroTask] = []

        if task.task_type == TaskType.DISCONNECT and result.get("success"):
            rxn = task.reaction_smiles or ""
            if rxn:
                new_tasks.append(create_validate_subtask(task, rxn))

        elif task.task_type == TaskType.VALIDATE:
            vs = (task.result or {}).get("validation_status", "")
            if vs == "PASS" and task.precursors:
                new_tasks.extend(expand_precursor_subtasks(task, task.precursors, result.get("sa_scores")))
            elif vs == "FAIL" and task.retry_count < task.max_retries:
                issues = result.get("issues") or []
                strs = [i.get("description", str(i)) if isinstance(i, dict) else str(i) for i in issues]
                new_tasks.append(create_repair_subtask(task, strs))

        elif task.task_type == TaskType.REPAIR and result.get("success"):
            repaired = result.get("repaired_smiles", task.reaction_smiles)
            if repaired:
                new_tasks.append(create_validate_subtask(task, repaired))

        if new_tasks:
            append_subtasks(route.tasks, new_tasks)

    # -- Internal: helpers -------------------------------------------------

    @staticmethod
    def _next_task(route: RetroRoute, allowed: Set[TaskType]) -> Optional[RetroTask]:
        for task in route.tasks.tasks:
            if task.status == TaskStatus.PENDING and task.task_type in allowed:
                return task
        return None

    def _get_route(self, route_id: str) -> RetroRoute:
        if self.route_manager is None:
            raise RuntimeError("No route_manager — call plan() first")
        route = self.route_manager.get_route(route_id)
        if route is None:
            raise ValueError(f"Route {route_id!r} not found")
        return route

    def _route_dir(self, route: RetroRoute) -> Path:
        return self.session_dir / route.metadata.get("route_name", route.route_id)

    def _ensure_managers(self, route: RetroRoute, route_dir: Path) -> None:
        if self.decision_manager is None:
            self.decision_manager = DecisionManager(route, route_dir)
        if self.journal is None:
            self.journal = ProcessJournal(route_dir / "process.md")

    def _log_step(self, task: RetroTask, result: Dict[str, Any]) -> None:
        if self.journal is None:
            return
        ok = result.get("success", False)
        err = result.get("error", "")
        self.journal.append_step(
            task,
            tool_call_summary=f"dispatch_skill({task.task_type.value})",
            tool_result_summary=f"{task.task_type.value} → {'OK' if ok else 'FAIL'}",
            llm_thinking="",
            conclusion=err if err else f"Completed {task.task_type.value}",
        )
