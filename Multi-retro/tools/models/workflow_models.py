"""
Workflow data models — task state, routes, and decision protocol.

Migrated from ``core/taskflow/task_state.py`` and ``decision_models.py``.
Key changes vs. the originals:

* Uses ``TaskStatus`` / ``TaskType`` from ``tools.common.status``
  (single source of truth).
* Removes task.md write-related fields (``title``, ``description``,
  ``validation_status``, ``subtasks``, ``images``, ``route_name``,
  ``strategy`` as a route-level field, etc.).
* ``DecisionContext.decision_type`` stored as plain ``str`` (not a
  separate enum) for simplicity — the 5 decision-type strings are
  still the same.
* Every class implements ``to_dict()`` / ``from_dict()`` with
  round-trip consistency and missing-field tolerance.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional

from tools.common.status import TaskStatus, TaskType


# ---------------------------------------------------------------------------
# RetroTask
# ---------------------------------------------------------------------------

@dataclass
class RetroTask:
    """Atomic unit of the retrosynthesis task flow."""

    task_id: str = ""
    task_type: TaskType = TaskType.ANALYZE
    status: TaskStatus = TaskStatus.PENDING
    target_smiles: str = ""
    parent_task_id: Optional[str] = None
    depth: int = 0
    result: Dict[str, Any] = field(default_factory=dict)
    retry_count: int = 0
    max_retries: int = 3
    selected_bond: Optional[Dict[str, Any]] = None
    reaction_smiles: str = ""
    precursors: List[str] = field(default_factory=list)

    def to_dict(self) -> Dict[str, Any]:
        return {
            "task_id": self.task_id,
            "task_type": self.task_type.value,
            "status": self.status.value,
            "target_smiles": self.target_smiles,
            "parent_task_id": self.parent_task_id,
            "depth": self.depth,
            "result": dict(self.result),
            "retry_count": self.retry_count,
            "max_retries": self.max_retries,
            "selected_bond": self.selected_bond,
            "reaction_smiles": self.reaction_smiles,
            "precursors": list(self.precursors),
        }

    @classmethod
    def from_dict(cls, d: Dict[str, Any]) -> "RetroTask":
        # Enum deserialization with fallback
        try:
            task_type = TaskType(d["task_type"]) if "task_type" in d else TaskType.ANALYZE
        except ValueError:
            task_type = TaskType.ANALYZE
        try:
            status = TaskStatus(d["status"]) if "status" in d else TaskStatus.PENDING
        except ValueError:
            status = TaskStatus.PENDING

        return cls(
            task_id=d.get("task_id", ""),
            task_type=task_type,
            status=status,
            target_smiles=d.get("target_smiles", ""),
            parent_task_id=d.get("parent_task_id"),
            depth=d.get("depth", 0),
            result=dict(d.get("result", {})) if d.get("result") is not None else {},
            retry_count=d.get("retry_count", 0),
            max_retries=d.get("max_retries", 3),
            selected_bond=d.get("selected_bond"),
            reaction_smiles=d.get("reaction_smiles", ""),
            precursors=list(d.get("precursors", [])),
        )


# ---------------------------------------------------------------------------
# RetroTaskList
# ---------------------------------------------------------------------------

@dataclass
class RetroTaskList:
    """Ordered collection of RetroTask objects with lookup helpers."""

    tasks: List[RetroTask] = field(default_factory=list)

    def get_task(self, task_id: str) -> Optional[RetroTask]:
        """Return the task with the given id, or ``None``."""
        for task in self.tasks:
            if task.task_id == task_id:
                return task
        return None

    def add_task(self, task: RetroTask) -> None:
        """Append a task to the list."""
        self.tasks.append(task)

    def get_pending_tasks(self) -> List[RetroTask]:
        """Return all tasks whose status is PENDING."""
        return [t for t in self.tasks if t.status == TaskStatus.PENDING]

    def to_dict(self) -> Dict[str, Any]:
        return {"tasks": [t.to_dict() for t in self.tasks]}

    @classmethod
    def from_dict(cls, d: Dict[str, Any]) -> "RetroTaskList":
        return cls(tasks=[RetroTask.from_dict(t) for t in d.get("tasks", [])])


# ---------------------------------------------------------------------------
# RetroRoute
# ---------------------------------------------------------------------------

@dataclass
class RetroRoute:
    """A complete synthesis route with its own task list and metadata."""

    route_id: str = ""
    target_smiles: str = ""
    tasks: RetroTaskList = field(default_factory=RetroTaskList)
    metadata: Dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> Dict[str, Any]:
        return {
            "route_id": self.route_id,
            "target_smiles": self.target_smiles,
            "tasks": self.tasks.to_dict(),
            "metadata": dict(self.metadata),
        }

    @classmethod
    def from_dict(cls, d: Dict[str, Any]) -> "RetroRoute":
        return cls(
            route_id=d.get("route_id", ""),
            target_smiles=d.get("target_smiles", ""),
            tasks=RetroTaskList.from_dict(d.get("tasks", {"tasks": []})),
            metadata=dict(d.get("metadata", {})),
        )


# ---------------------------------------------------------------------------
# DecisionContext
# ---------------------------------------------------------------------------

@dataclass
class DecisionContext:
    """Structured context output at a decision point for the LLM host."""

    decision_type: str = ""
    task_id: str = ""
    context: Dict[str, Any] = field(default_factory=dict)
    available_actions: List[Dict[str, Any]] = field(default_factory=list)
    decision_history: List[Dict[str, Any]] = field(default_factory=list)
    exploration_tools: List[Dict[str, Any]] = field(default_factory=list)
    exploration_budget: int = 5

    def to_dict(self) -> Dict[str, Any]:
        return {
            "decision_type": self.decision_type,
            "task_id": self.task_id,
            "context": self.context,
            "available_actions": list(self.available_actions),
            "decision_history": list(self.decision_history),
            "exploration_tools": list(self.exploration_tools),
            "exploration_budget": self.exploration_budget,
        }

    @classmethod
    def from_dict(cls, d: Dict[str, Any]) -> "DecisionContext":
        return cls(
            decision_type=d.get("decision_type", ""),
            task_id=d.get("task_id", ""),
            context=dict(d.get("context", {})),
            available_actions=list(d.get("available_actions", [])),
            decision_history=list(d.get("decision_history", [])),
            exploration_tools=list(d.get("exploration_tools", [])),
            exploration_budget=d.get("exploration_budget", 5),
        )


# ---------------------------------------------------------------------------
# DecisionInstruction
# ---------------------------------------------------------------------------

@dataclass
class DecisionInstruction:
    """Decision instruction submitted by the LLM host back to TaskFlow."""

    task_id: str = ""
    action: str = ""
    params: Dict[str, Any] = field(default_factory=dict)
    reasoning: Optional[str] = None
    exploration_log: Optional[List[str]] = None
    reaction_conditions: Optional[Dict[str, Any]] = None

    def to_dict(self) -> Dict[str, Any]:
        d: Dict[str, Any] = {
            "task_id": self.task_id,
            "action": self.action,
            "params": dict(self.params),
            "reasoning": self.reasoning,
            "exploration_log": list(self.exploration_log) if self.exploration_log is not None else None,
            "reaction_conditions": self.reaction_conditions,
        }
        return d

    @classmethod
    def from_dict(cls, d: Dict[str, Any]) -> "DecisionInstruction":
        exp_log = d.get("exploration_log")
        return cls(
            task_id=d.get("task_id", ""),
            action=d.get("action", ""),
            params=dict(d.get("params", {})),
            reasoning=d.get("reasoning"),
            exploration_log=list(exp_log) if exp_log is not None else None,
            reaction_conditions=d.get("reaction_conditions"),
        )
