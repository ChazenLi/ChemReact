"""
Task state management — operational logic for RetroTask, RetroTaskList, RetroRoute.

This module provides the operational layer on top of the pure data models
defined in ``tools.models.workflow_models``.  It re-exports the data classes
for convenience and adds higher-level helpers used by the orchestrator and
other workflow tools (filtering, status transitions, subtask queries, etc.).

Key differences from the original ``core/taskflow/task_state.py``:

* Uses ``TaskStatus`` / ``TaskType`` from ``tools.common.status`` (unified enums).
* All task.md rendering / writing logic has been removed — that responsibility
  belongs to ``tools/output/status_renderer.py``.
* Data definitions live in ``tools.models.workflow_models``; this module
  contains only operational helpers.
"""

from __future__ import annotations

import logging
from datetime import datetime
from typing import Any, Dict, List, Optional, Set

from tools.common.status import TaskStatus, TaskType
from tools.models.workflow_models import (
    RetroRoute,
    RetroTask,
    RetroTaskList,
)

__all__ = [
    "RetroTask",
    "RetroTaskList",
    "RetroRoute",
    "TaskStatus",
    "TaskType",
    "get_next_pending",
    "get_next_pending_excluding",
    "get_subtasks",
    "get_tasks_by_status",
    "get_tasks_by_type",
    "transition_status",
    "complete_task",
    "fail_task",
    "start_task",
    "create_task",
    "next_subtask_id",
]

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Query helpers
# ---------------------------------------------------------------------------

def get_next_pending(task_list: RetroTaskList) -> Optional[RetroTask]:
    """Return the first task whose status is PENDING, or ``None``."""
    for task in task_list.tasks:
        if task.status == TaskStatus.PENDING:
            return task
    return None


def get_next_pending_excluding(
    task_list: RetroTaskList,
    excluded_types: Set[TaskType],
) -> Optional[RetroTask]:
    """Return the first PENDING task whose type is *not* in *excluded_types*."""
    for task in task_list.tasks:
        if task.status == TaskStatus.PENDING and task.task_type not in excluded_types:
            return task
    return None


def get_subtasks(task_list: RetroTaskList, parent_id: str) -> List[RetroTask]:
    """Return all tasks whose ``parent_task_id`` matches *parent_id*."""
    return [t for t in task_list.tasks if t.parent_task_id == parent_id]


def get_tasks_by_status(
    task_list: RetroTaskList, status: TaskStatus
) -> List[RetroTask]:
    """Return all tasks with the given *status*."""
    return [t for t in task_list.tasks if t.status == status]


def get_tasks_by_type(
    task_list: RetroTaskList, task_type: TaskType
) -> List[RetroTask]:
    """Return all tasks with the given *task_type*."""
    return [t for t in task_list.tasks if t.task_type == task_type]


# ---------------------------------------------------------------------------
# Status transition helpers
# ---------------------------------------------------------------------------

# Valid transitions — guards against illegal state changes.
_VALID_TRANSITIONS: Dict[TaskStatus, Set[TaskStatus]] = {
    TaskStatus.PENDING: {
        TaskStatus.IN_PROGRESS,
        TaskStatus.SKIPPED,
        TaskStatus.BLOCKED,
    },
    TaskStatus.IN_PROGRESS: {
        TaskStatus.COMPLETED,
        TaskStatus.FAILED,
        TaskStatus.AWAITING_DECISION,
        TaskStatus.VALIDATED,
    },
    TaskStatus.AWAITING_DECISION: {
        TaskStatus.IN_PROGRESS,
        TaskStatus.COMPLETED,
        TaskStatus.FAILED,
        TaskStatus.SKIPPED,
    },
    TaskStatus.VALIDATED: {
        TaskStatus.COMPLETED,
        TaskStatus.IN_PROGRESS,  # re-validation after repair
    },
    TaskStatus.FAILED: {
        TaskStatus.PENDING,  # retry
        TaskStatus.SKIPPED,
    },
    TaskStatus.BLOCKED: {
        TaskStatus.PENDING,
    },
    TaskStatus.COMPLETED: set(),  # terminal
    TaskStatus.SKIPPED: set(),    # terminal
}


def transition_status(
    task: RetroTask,
    new_status: TaskStatus,
    *,
    force: bool = False,
) -> None:
    """Transition *task* to *new_status*.

    Raises ``ValueError`` if the transition is not allowed unless *force*
    is ``True``.
    """
    if not force:
        allowed = _VALID_TRANSITIONS.get(task.status, set())
        if new_status not in allowed:
            raise ValueError(
                f"Invalid transition: {task.status.value!r} → {new_status.value!r} "
                f"for task {task.task_id!r}"
            )
    old = task.status
    task.status = new_status
    logger.debug(
        "Task %s: %s → %s", task.task_id, old.value, new_status.value
    )


def start_task(task: RetroTask) -> None:
    """Mark *task* as IN_PROGRESS."""
    transition_status(task, TaskStatus.IN_PROGRESS)


def complete_task(task: RetroTask) -> None:
    """Mark *task* as COMPLETED."""
    transition_status(task, TaskStatus.COMPLETED)


def fail_task(task: RetroTask) -> None:
    """Mark *task* as FAILED and bump retry_count."""
    transition_status(task, TaskStatus.FAILED)
    task.retry_count += 1


# ---------------------------------------------------------------------------
# Task creation helpers
# ---------------------------------------------------------------------------

def create_task(
    task_id: str,
    task_type: TaskType,
    target_smiles: str,
    *,
    parent_task_id: Optional[str] = None,
    depth: int = 0,
) -> RetroTask:
    """Factory for creating a new ``RetroTask`` with sensible defaults."""
    return RetroTask(
        task_id=task_id,
        task_type=task_type,
        status=TaskStatus.PENDING,
        target_smiles=target_smiles,
        parent_task_id=parent_task_id,
        depth=depth,
    )


def next_subtask_id(task_list: RetroTaskList, parent_id: str) -> str:
    """Generate the next hierarchical subtask id under *parent_id*.

    E.g. if parent_id is ``"1"`` and there are already subtasks ``"1.1"``
    and ``"1.2"``, this returns ``"1.3"``.
    """
    existing = get_subtasks(task_list, parent_id)
    next_idx = len(existing) + 1
    return f"{parent_id}.{next_idx}"
