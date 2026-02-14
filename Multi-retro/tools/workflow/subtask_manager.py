"""
Subtask management — create and append child tasks for the retrosynthesis pipeline.

Extracted from ``core/taskflow/executor.py`` (_append_subtasks,
_create_validate_subtask, _create_repair_subtask, expand_precursor_subtasks).

Public API:
    create_validate_subtask(parent_task, reaction_smiles) -> RetroTask
    create_repair_subtask(parent_task, issues)            -> RetroTask
    expand_precursor_subtasks(parent_task, precursors, sa_scores) -> List[RetroTask]
    append_subtasks(task_list, subtasks)                  -> None
"""

from __future__ import annotations

import logging
import re
from typing import Dict, List, Optional

from tools.common.constants import RetroLimits, SAThresholds
from tools.common.status import TaskStatus, TaskType
from tools.models.workflow_models import RetroTask, RetroTaskList

__all__ = [
    "create_validate_subtask",
    "create_repair_subtask",
    "expand_precursor_subtasks",
    "append_subtasks",
]

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _clean_smiles(smiles: str) -> str:
    """Remove text annotations from a SMILES string (e.g. ' [apply: C=O → C-OH]')."""
    if not smiles:
        return smiles
    match = re.match(r"^([^\s\[]*(?:\[[^\]]*\])*[^\s\[]*)", smiles)
    return match.group(1) if match else smiles.split()[0]


def _estimate_mw(smiles: str) -> float:
    """Estimate molecular weight from SMILES.

    Uses RDKit when available, otherwise falls back to a rough heuristic
    of ~14 Da per heavy-atom character.
    """
    try:
        from rdkit import Chem
        from rdkit.Chem import Descriptors

        mol = Chem.MolFromSmiles(smiles)
        if mol is not None:
            return Descriptors.MolWt(mol)
    except Exception:
        pass
    heavy = sum(1 for c in smiles if c.isalpha() and c.isupper())
    return heavy * 14.0


def _normalize_reaction_smiles(
    reaction_smiles: Optional[str],
    precursors: List[str],
    target_smiles: str,
) -> str:
    """Normalize reaction SMILES to ``reactant1.reactant2>>product`` format.

    If *reaction_smiles* already contains ``>>``, return as-is.
    Otherwise rebuild from *precursors* and *target_smiles*.
    """
    if reaction_smiles and ">>" in reaction_smiles:
        return reaction_smiles
    if precursors and target_smiles:
        reactants = ".".join(precursors)
        return f"{reactants}>>{target_smiles}"
    return reaction_smiles or ""


def _get_sa_score(smiles: str, sa_scores: Optional[Dict[str, float]] = None) -> float:
    """Return SA score from the provided mapping, or estimate it."""
    if sa_scores and smiles in sa_scores:
        return sa_scores[smiles]
    try:
        from tools.chem.sa_scorer import estimate_sa_score
        return estimate_sa_score(smiles)
    except Exception:
        return 3.0  # conservative default


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def create_validate_subtask(
    parent_task: RetroTask,
    reaction_smiles: str,
) -> RetroTask:
    """Create a VALIDATE subtask for a completed DISCONNECT task.

    Args:
        parent_task: The parent DISCONNECT task.
        reaction_smiles: The reaction SMILES to validate.  If it does not
            contain ``>>``, it is rebuilt from the parent's precursors and
            target SMILES.

    Returns:
        A new ``RetroTask`` with type VALIDATE and status PENDING.
    """
    normalized_rxn = _normalize_reaction_smiles(
        reaction_smiles,
        parent_task.precursors,
        parent_task.target_smiles,
    )

    validate_task = RetroTask(
        task_id="",  # caller assigns via next_subtask_id
        task_type=TaskType.VALIDATE,
        status=TaskStatus.PENDING,
        target_smiles=parent_task.target_smiles,
        parent_task_id=parent_task.task_id,
        depth=parent_task.depth + 1,
        precursors=list(parent_task.precursors),
        reaction_smiles=normalized_rxn,
    )
    return validate_task


def create_repair_subtask(
    parent_task: RetroTask,
    issues: List[str],
) -> RetroTask:
    """Create a REPAIR subtask for a failed validation.

    Args:
        parent_task: The parent task (typically VALIDATE or DISCONNECT)
            whose reaction needs repair.
        issues: List of issue descriptions from the validation result.

    Returns:
        A new ``RetroTask`` with type REPAIR and status PENDING.
        The ``result`` field stores the issues for downstream consumption.
    """
    repair_task = RetroTask(
        task_id="",  # caller assigns via next_subtask_id
        task_type=TaskType.REPAIR,
        status=TaskStatus.PENDING,
        target_smiles=parent_task.target_smiles,
        parent_task_id=parent_task.task_id,
        depth=parent_task.depth + 1,
        reaction_smiles=parent_task.reaction_smiles,
        result={"issues": list(issues)},
    )
    return repair_task


def expand_precursor_subtasks(
    parent_task: RetroTask,
    precursors: List[str],
    sa_scores: Optional[Dict[str, float]] = None,
) -> List[RetroTask]:
    """Expand precursors into child subtasks based on SA score and MW.

    For each precursor SMILES:
    - SA < PURCHASABLE_MAX or MW < MW_TERMINAL → terminal AVAILABILITY task
      (status SKIPPED).
    - SA in [PURCHASABLE_MAX, COMPLEX_MIN) → optional STRATEGY → ANALYZE →
      DISCONNECT chain.
    - SA >= COMPLEX_MIN → required STRATEGY → ANALYZE → DISCONNECT chain.
    - Precursors identical to the parent target (after cleaning) are skipped
      to avoid no-op disconnection cycles.
    - Depth limit (``RetroLimits.MAX_DEPTH``) is enforced.

    Args:
        parent_task: The parent task whose precursors are being expanded.
        precursors: List of precursor SMILES strings.
        sa_scores: Optional pre-computed SA scores keyed by SMILES.
            If a precursor is missing from the map, its score is estimated.

    Returns:
        List of newly created ``RetroTask`` objects.  The caller is
        responsible for assigning ``task_id`` values and adding them to
        the task list via :func:`append_subtasks`.
    """
    child_depth = parent_task.depth + 1

    if child_depth >= RetroLimits.MAX_DEPTH:
        logger.warning(
            "Depth limit reached (%d >= %d) for task %s; not expanding precursors",
            child_depth,
            RetroLimits.MAX_DEPTH,
            parent_task.task_id,
        )
        return []

    parent_target_clean = _clean_smiles(parent_task.target_smiles)
    new_tasks: List[RetroTask] = []

    for smiles in precursors:
        # Skip precursors identical to parent target (no-op disconnect)
        if _clean_smiles(smiles) == parent_target_clean:
            logger.info(
                "Skipping precursor identical to parent target for task %s: %s",
                parent_task.task_id,
                smiles[:60],
            )
            continue

        sa_score = _get_sa_score(smiles, sa_scores)
        mw = _estimate_mw(smiles)

        if mw < RetroLimits.MW_TERMINAL or sa_score < SAThresholds.PURCHASABLE_MAX:
            # Terminal precursor
            terminal = RetroTask(
                task_id="",  # caller assigns
                task_type=TaskType.AVAILABILITY,
                status=TaskStatus.SKIPPED,
                target_smiles=smiles,
                parent_task_id=parent_task.task_id,
                depth=child_depth,
                result={
                    "sa_score": round(sa_score, 2),
                    "mw": round(mw, 1),
                    "terminal": True,
                },
            )
            new_tasks.append(terminal)
        else:
            # Non-terminal: STRATEGY → ANALYZE → DISCONNECT chain
            is_optional = sa_score < SAThresholds.COMPLEX_MIN
            sa_desc = (
                f"SA={sa_score:.2f}, MW={mw:.1f} "
                f"({'optional' if is_optional else 'required'})"
            )

            strategy_task = RetroTask(
                task_id="",
                task_type=TaskType.STRATEGY,
                status=TaskStatus.PENDING,
                target_smiles=smiles,
                parent_task_id=parent_task.task_id,
                depth=child_depth,
                result={"sa_desc": sa_desc},
            )

            analyze_task = RetroTask(
                task_id="",
                task_type=TaskType.ANALYZE,
                status=TaskStatus.PENDING,
                target_smiles=smiles,
                parent_task_id=parent_task.task_id,
                depth=child_depth,
            )

            disconnect_task = RetroTask(
                task_id="",
                task_type=TaskType.DISCONNECT,
                status=TaskStatus.PENDING,
                target_smiles=smiles,
                parent_task_id=parent_task.task_id,
                depth=child_depth,
            )

            # Order matters: STRATEGY → ANALYZE → DISCONNECT
            new_tasks.extend([strategy_task, analyze_task, disconnect_task])

    return new_tasks


def append_subtasks(
    task_list: RetroTaskList,
    subtasks: List[RetroTask],
) -> None:
    """Add a list of subtasks to a task list, auto-assigning task IDs.

    Each subtask must have ``parent_task_id`` set.  The ``task_id`` is
    generated as ``"{parent_task_id}.{n}"`` where *n* is the next
    available index under that parent.

    Args:
        task_list: The task list to append to.
        subtasks: Subtasks to add.  Their ``task_id`` fields will be
            overwritten with auto-generated hierarchical IDs.
    """
    if not subtasks:
        return

    # Pre-compute next index per parent to avoid repeated scans
    parent_counts: Dict[str, int] = {}

    for subtask in subtasks:
        parent_id = subtask.parent_task_id or ""
        if parent_id not in parent_counts:
            existing = [
                t for t in task_list.tasks if t.parent_task_id == parent_id
            ]
            parent_counts[parent_id] = len(existing)

        parent_counts[parent_id] += 1
        subtask.task_id = f"{parent_id}.{parent_counts[parent_id]}"
        task_list.add_task(subtask)

    logger.debug(
        "Appended %d subtask(s) to task list (total: %d)",
        len(subtasks),
        len(task_list.tasks),
    )
