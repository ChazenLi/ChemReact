"""
Skill dispatch — route a RetroTask to the correct Skill and manage
argument construction / result extraction.

Extracted from ``core/taskflow/executor.py`` (_dispatch_skill,
_build_skill_args, _extract_result).

Key improvements over the original:

* Uses typed dataclass arguments (``AnalyzeArgs``, ``BreakBondArgs``,
  ``ValidateArgs``, ``RepairArgs``) instead of ``SimpleNamespace``.
* Pure functions — no dependency on ``TaskExecutor`` instance state.
* ``dispatch_skill`` accepts an explicit *skill_registry* mapping so
  callers can inject / override skills freely.
"""

from __future__ import annotations

import logging
from types import SimpleNamespace
from typing import Any, Callable, Dict, Optional, Union

from tools.common.status import TaskType
from tools.models.skill_models import (
    AnalyzeArgs,
    BreakBondArgs,
    RepairArgs,
    ValidateArgs,
)
from tools.models.workflow_models import RetroTask

__all__ = [
    "dispatch_skill",
    "build_skill_args",
    "extract_result",
]

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _parse_selected_bond(selected_bond: Any) -> Optional[tuple]:
    """Parse *selected_bond* into ``(atom1_idx, atom2_idx)`` or ``None``.

    Supports both the legacy string format ``"3-7"`` and the new dict
    format ``{"atom1_idx": 3, "atom2_idx": 7}``.
    """
    if selected_bond is None:
        return None

    if isinstance(selected_bond, dict):
        try:
            return (int(selected_bond["atom1_idx"]), int(selected_bond["atom2_idx"]))
        except (KeyError, ValueError, TypeError):
            return None

    if isinstance(selected_bond, str):
        parts = selected_bond.split("-")
        if len(parts) == 2:
            try:
                return (int(parts[0]), int(parts[1]))
            except (ValueError, TypeError):
                return None

    return None


# ---------------------------------------------------------------------------
# build_skill_args
# ---------------------------------------------------------------------------

def build_skill_args(
    task: RetroTask,
    *,
    fallback_target: str = "",
) -> Union[AnalyzeArgs, BreakBondArgs, ValidateArgs, RepairArgs, SimpleNamespace]:
    """Build typed Skill arguments from a *task*.

    Returns a typed dataclass for the four core task types (ANALYZE,
    DISCONNECT with selected_bond, VALIDATE, REPAIR).  For other types
    or fallback scenarios a ``SimpleNamespace`` is returned to preserve
    backward compatibility with skills that still expect namespace args.

    Parameters
    ----------
    task:
        The task to build arguments for.
    fallback_target:
        Target SMILES to use when ``task.target_smiles`` is empty
        (typically the route-level target).
    """
    target = task.target_smiles or fallback_target

    # -- STRATEGY / ANALYZE → AnalyzeArgs
    if task.task_type in (TaskType.STRATEGY, TaskType.ANALYZE):
        return AnalyzeArgs(smiles=target)

    # -- DISCONNECT
    if task.task_type == TaskType.DISCONNECT:
        bond = _parse_selected_bond(task.selected_bond)
        if bond is not None:
            return BreakBondArgs(smiles=target, atom1_idx=bond[0], atom2_idx=bond[1])
        # No selected_bond or parse failure → ProposeDisconnectionSkill path
        if task.selected_bond is not None:
            logger.warning(
                "Failed to parse selected_bond %r for task %s; "
                "falling back to ProposeDisconnectionSkill path",
                task.selected_bond,
                task.task_id,
            )
        return SimpleNamespace(target_smiles=target)

    # -- VALIDATE → ValidateArgs
    if task.task_type == TaskType.VALIDATE:
        return ValidateArgs(reaction_smiles=task.reaction_smiles or "")

    # -- REPAIR → RepairArgs
    if task.task_type == TaskType.REPAIR:
        issues: list = []
        if task.result and isinstance(task.result, dict):
            issues = task.result.get("issues", [])
            if not issues:
                issues = task.result.get("verdict_details", [])
        return RepairArgs(
            reaction_smiles=task.reaction_smiles or "",
            issues=issues,
        )

    # -- AVAILABILITY
    if task.task_type == TaskType.AVAILABILITY:
        return {"smiles_list": task.precursors or []}

    # -- REPORT
    if task.task_type == TaskType.REPORT:
        return {
            "retro_graph": {},
            "output_dir": "",
        }

    # -- Fallback
    return SimpleNamespace(smiles=target, target_smiles=target)


# ---------------------------------------------------------------------------
# dispatch_skill
# ---------------------------------------------------------------------------

def dispatch_skill(
    task: RetroTask,
    skill_registry: Dict[str, Any],
    *,
    fallback_target: str = "",
) -> Dict[str, Any]:
    """Dispatch *task* to the appropriate Skill and return the raw result.

    Parameters
    ----------
    task:
        The task to execute.
    skill_registry:
        Mapping of ``TaskType.value`` (or special keys like
        ``"disconnect_break_bond"`` / ``"disconnect_propose"``) to Skill
        instances or callables.  Each entry must support
        ``skill.execute(args)`` returning a ``Dict[str, Any]``.
    fallback_target:
        Fallback target SMILES (see ``build_skill_args``).

    Returns
    -------
    Dict[str, Any]
        The raw result dictionary from the Skill's ``execute()`` method.

    Raises
    ------
    ValueError
        If no skill is registered for the task's type.
    """
    # Resolve the skill
    skill = _resolve_skill(task, skill_registry)
    if skill is None:
        raise ValueError(
            f"No skill registered for task type {task.task_type.value!r}"
        )

    args = build_skill_args(task, fallback_target=fallback_target)
    return skill.execute(args)


def _resolve_skill(task: RetroTask, registry: Dict[str, Any]) -> Any:
    """Pick the right skill from *registry* for *task*.

    For DISCONNECT tasks the registry may contain two entries:
    * ``"disconnect_break_bond"`` — used when ``selected_bond`` is present.
    * ``"disconnect_propose"``    — used otherwise.
    If neither is found, falls back to the generic ``"disconnect"`` key.
    """
    if task.task_type == TaskType.DISCONNECT:
        bond = _parse_selected_bond(task.selected_bond)
        if bond is not None:
            skill = registry.get("disconnect_break_bond")
            if skill is not None:
                return skill
        else:
            skill = registry.get("disconnect_propose")
            if skill is not None:
                return skill
        # Generic fallback
        return registry.get(TaskType.DISCONNECT.value)

    return registry.get(task.task_type.value)


# ---------------------------------------------------------------------------
# extract_result
# ---------------------------------------------------------------------------

def extract_result(task: RetroTask, raw_result: Dict[str, Any]) -> None:
    """Extract key fields from a Skill's raw result and write them to *task*.

    Mutates *task* in place:

    * **DISCONNECT** (success): populates ``task.precursors``,
      ``task.reaction_smiles``, and optionally ``task.selected_bond``.
    * **VALIDATE**: sets ``task.result["validation_status"]`` to
      ``"PASS"`` or ``"FAIL"``.

    Parameters
    ----------
    task:
        The task to update.
    raw_result:
        The dictionary returned by the Skill's ``execute()`` method.
    """
    if not isinstance(raw_result, dict):
        return

    if task.task_type == TaskType.DISCONNECT and raw_result.get("success"):
        _extract_disconnect_result(task, raw_result)

    elif task.task_type == TaskType.VALIDATE:
        validation_status = "PASS" if raw_result.get("is_valid") else "FAIL"
        task.result = {**task.result, "validation_status": validation_status}


def _extract_disconnect_result(task: RetroTask, raw_result: Dict[str, Any]) -> None:
    """Handle DISCONNECT result extraction (both propose and break-bond paths)."""
    # ProposeDisconnectionSkill returns "disconnections" or legacy "proposals"
    proposals = raw_result.get("disconnections") or raw_result.get("proposals")
    if proposals and isinstance(proposals, list):
        best = proposals[0]
        task.precursors = best.get("fragments") or best.get("precursors", [])
        task.reaction_smiles = best.get("reaction_smiles", "")
        if "atom1_idx" in best:
            task.selected_bond = {
                "atom1_idx": best["atom1_idx"],
                "atom2_idx": best["atom2_idx"],
            }
        return

    # BreakBondSkill standard path — suggested_fragments preferred
    if "fragments" in raw_result or "suggested_fragments" in raw_result:
        task.precursors = (
            raw_result.get("suggested_fragments")
            or raw_result.get("fragments", [])
        )
        task.reaction_smiles = (
            raw_result.get("suggested_reaction_smiles")
            or raw_result.get("reaction_smiles", "")
        )
