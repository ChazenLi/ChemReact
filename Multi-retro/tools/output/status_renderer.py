"""
Status renderer — Markdown status tables for retrosynthesis tasks and routes.

Migrated from ``core/taskflow/renderer.py``.  Key changes vs. the original:

* Works with the new ``tools.models.workflow_models`` data classes (no
  ``route_name``, ``strategy``, ``images``, ``title``, ``subtasks`` fields).
* Provides two focused public functions instead of the monolithic
  ``TaskRenderer`` class:
  - ``render_status_table(task_list)`` — Markdown table for a task list.
  - ``render_route_status(route)`` — full route-level status report.
* Uses ``TaskStatus`` / ``TaskType`` from ``tools.common.status``.

Usage::

    from tools.output.status_renderer import render_status_table, render_route_status

    md = render_status_table(task_list)
    md = render_route_status(route)
"""

from __future__ import annotations

from typing import Dict, List, Tuple

from tools.common.status import TaskStatus, TaskType
from tools.models.workflow_models import RetroRoute, RetroTask, RetroTaskList

__all__ = ["render_status_table", "render_route_status"]

# Phase grouping: (phase number, display name, set of TaskType members)
_PHASE_GROUPS: List[Tuple[int, str, set]] = [
    (0, "侦察", {TaskType.STRATEGY, TaskType.ANALYZE}),
    (1, "断裂", {TaskType.DISCONNECT, TaskType.VALIDATE, TaskType.REPAIR}),
    (2, "可用性", {TaskType.AVAILABILITY}),
    (3, "收尾", {TaskType.REPORT}),
]

_CHECKED_STATUSES = {TaskStatus.COMPLETED, TaskStatus.VALIDATED}

# Human-readable status labels
_STATUS_LABELS: Dict[TaskStatus, str] = {
    TaskStatus.PENDING: "⬜ 待执行",
    TaskStatus.IN_PROGRESS: "🔄 进行中",
    TaskStatus.AWAITING_DECISION: "⏸️ 等待决策",
    TaskStatus.VALIDATED: "✅ 已验证",
    TaskStatus.COMPLETED: "✅ 已完成",
    TaskStatus.FAILED: "❌ 失败",
    TaskStatus.SKIPPED: "⏭️ 跳过",
    TaskStatus.BLOCKED: "🚫 阻塞",
}


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _task_id_sort_key(task_id: str) -> List[int]:
    """Convert ``'1.2.3'`` to ``[1, 2, 3]`` for natural sorting."""
    try:
        return [int(x) for x in task_id.split(".")]
    except ValueError:
        return [0]


def _status_icon(status: TaskStatus) -> str:
    """Return a checkbox marker for the given status."""
    return "[x]" if status in _CHECKED_STATUSES else "[ ]"


def _compute_stats(tasks: List[RetroTask]) -> Dict[str, int]:
    """Compute summary statistics for a list of tasks."""
    total = len(tasks)
    completed = sum(1 for t in tasks if t.status in _CHECKED_STATUSES)
    in_progress = sum(1 for t in tasks if t.status == TaskStatus.IN_PROGRESS)
    pending = sum(1 for t in tasks if t.status == TaskStatus.PENDING)
    failed = sum(1 for t in tasks if t.status == TaskStatus.FAILED)
    max_depth = max((t.depth for t in tasks), default=0)
    terminal_precursors = sum(
        1 for t in tasks if t.status == TaskStatus.SKIPPED
    )
    return {
        "total": total,
        "completed": completed,
        "in_progress": in_progress,
        "pending": pending,
        "failed": failed,
        "max_depth": max_depth,
        "terminal_precursors": terminal_precursors,
    }


def _render_task_row(task: RetroTask, indent_level: int = 0) -> str:
    """Render a single task as a Markdown checkbox line."""
    indent = "  " * indent_level
    checkbox = _status_icon(task.status)
    type_tag = task.task_type.value
    smiles_part = f" `{task.target_smiles}`" if task.target_smiles else ""
    return f"{indent}- {checkbox} **{task.task_id} [{type_tag}]**{smiles_part}"


def _render_task_tree(
    task: RetroTask,
    all_tasks: List[RetroTask],
    phase_task_ids: set,
    lines: List[str],
) -> None:
    """Render a task and its sub-tasks recursively."""
    lines.append(_render_task_row(task, indent_level=task.depth))

    # Render children within the same phase
    children = [
        t
        for t in all_tasks
        if t.parent_task_id == task.task_id and t.task_id in phase_task_ids
    ]
    children.sort(key=lambda t: _task_id_sort_key(t.task_id))
    for child in children:
        _render_task_tree(child, all_tasks, phase_task_ids, lines)


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def render_status_table(task_list: RetroTaskList) -> str:
    """Render a Markdown status table for a :class:`RetroTaskList`.

    The output groups tasks by phase (侦察 → 断裂 → 可用性 → 收尾),
    renders each task as a checkbox line with hierarchical indentation,
    and appends a statistics summary table.

    Args:
        task_list: The task list to render.

    Returns:
        A Markdown-formatted string.
    """
    lines: List[str] = []
    tasks = task_list.tasks

    # Group tasks by phase
    for phase_num, phase_name, phase_types in _PHASE_GROUPS:
        phase_tasks = [t for t in tasks if t.task_type in phase_types]
        if not phase_tasks:
            continue

        lines.append(f"## Phase {phase_num}: {phase_name}")
        lines.append("")

        phase_task_ids = {t.task_id for t in phase_tasks}

        # Top-level tasks: no parent, or parent outside this phase
        top_level = [
            t
            for t in phase_tasks
            if t.parent_task_id is None
            or t.parent_task_id not in phase_task_ids
        ]
        top_level.sort(key=lambda t: _task_id_sort_key(t.task_id))

        for task in top_level:
            _render_task_tree(task, tasks, phase_task_ids, lines)

        lines.append("")

    # Statistics summary
    stats = _compute_stats(tasks)
    lines.append("## 统计摘要")
    lines.append("")
    lines.append("| 指标 | 值 |")
    lines.append("|------|-----|")
    lines.append(f"| 总任务数 | {stats['total']} |")
    lines.append(f"| 已完成 | {stats['completed']} |")
    lines.append(f"| 进行中 | {stats['in_progress']} |")
    lines.append(f"| 待执行 | {stats['pending']} |")
    lines.append(f"| 失败 | {stats['failed']} |")
    lines.append(f"| 最大深度 | {stats['max_depth']} |")
    lines.append(f"| 终端前体数 | {stats['terminal_precursors']} |")
    lines.append("")

    return "\n".join(lines)


def render_route_status(route: RetroRoute) -> str:
    """Render a full Markdown status report for a :class:`RetroRoute`.

    Includes route metadata (route ID, target SMILES), the task status
    table, and summary statistics.

    Args:
        route: The route to render.

    Returns:
        A Markdown-formatted string.
    """
    lines: List[str] = []

    # Title
    lines.append(f"# 逆合成路线状态: {route.route_id}")
    lines.append("")

    # Metadata
    lines.append("## 元数据")
    lines.append("")
    lines.append(f"- **Target SMILES**: `{route.target_smiles}`")
    lines.append(f"- **Route ID**: {route.route_id}")
    if route.metadata:
        for key, value in route.metadata.items():
            lines.append(f"- **{key}**: {value}")
    lines.append("")

    # Task status table
    task_table = render_status_table(route.tasks)
    lines.append(task_table)

    return "\n".join(lines)
