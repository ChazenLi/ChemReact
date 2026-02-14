"""
Unified Status Enums
====================

Single source of truth for task status and task type enumerations.
Replaces string-based status comparisons throughout the codebase and
the duplicated enums previously in core/taskflow/task_state.py.

Usage::

    from tools.common.status import TaskStatus, TaskType

    if task.status == TaskStatus.PENDING:
        ...
    if task.task_type == TaskType.VALIDATE:
        ...
"""

from enum import Enum


class TaskStatus(Enum):
    """Status of a retrosynthesis task.

    Lifecycle: PENDING → IN_PROGRESS → COMPLETED / FAILED
    Special states:
        AWAITING_DECISION — paused at a decision point, waiting for LLM input.
        VALIDATED — reaction passed validation, ready for next phase.
        SKIPPED — intentionally bypassed (e.g. unnecessary repair).
        BLOCKED — cannot proceed due to unresolved dependency.
    """

    PENDING = "pending"
    IN_PROGRESS = "in_progress"
    AWAITING_DECISION = "awaiting_decision"
    VALIDATED = "validated"
    COMPLETED = "completed"
    FAILED = "failed"
    SKIPPED = "skipped"
    BLOCKED = "blocked"


class RouteStatus(Enum):
    """Status of a synthesis route.

    Lifecycle: PLANNING → COMPLETED / ABANDONED
    """

    PLANNING = "planning"
    COMPLETED = "completed"
    ABANDONED = "abandoned"
    PARTIAL = "partial"
    FAILED = "failed"


class TaskType(Enum):
    """Type of a retrosynthesis task.

    Each value maps to a specific skill or tool invocation in the
    orchestration pipeline:

        STRATEGY     — decide_strategy / get_global_strategy
        ANALYZE      — analyze_molecule
        DISCONNECT   — propose_disconnection / break_bond
        VALIDATE     — validate_reaction
        REPAIR       — repair_reaction
        AVAILABILITY — check_availability
        REPORT       — render_report / generate_synthesis_report
    """

    STRATEGY = "strategy"
    ANALYZE = "analyze"
    DISCONNECT = "disconnect"
    VALIDATE = "validate"
    REPAIR = "repair"
    AVAILABILITY = "availability"
    REPORT = "report"
