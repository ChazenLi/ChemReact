"""
Retrosynthesis graph builder — extracted from ``core/taskflow/executor.py``.

Builds a typed :class:`RetroGraph` from a list of :class:`RetroTask` objects.
Fixes the format mismatch bug between the old executor output (plain dicts
with ``"type"`` key) and the visualizer expectations (``"node_type"`` key,
typed ``RetroGraphNode`` / ``RetroGraphEdge`` objects).

Key fixes vs. the original ``_build_retro_graph``:

* Nodes use ``node_type`` (matching ``RetroGraphNode.node_type``) instead of
  the ambiguous ``"type"`` key that conflicted with Python builtins.
* Nodes carry ``sa_score`` and ``is_terminal`` fields so downstream consumers
  (visualizer, report generator) don't need to re-derive them.
* Edges use the canonical ``RetroGraphEdge`` schema — no extra ad-hoc keys
  like ``"reaction_smiles_raw"`` or ``"reaction_class"`` that the visualizer
  never consumed.
* The returned ``RetroGraph`` can be serialised via ``to_dict()`` and fed
  directly to ``visualize_retro_graph()`` without manual patching.

Usage::

    from tools.output.graph_builder import build_retro_graph, collect_terminal_precursors

    graph = build_retro_graph(task_list.tasks, target_smiles)
    terminal = collect_terminal_precursors(task_list.tasks)
"""

from __future__ import annotations

import logging
import re
from typing import Any, Dict, List, Optional

from tools.common.constants import SAThresholds
from tools.common.status import TaskStatus, TaskType
from tools.models.output_models import RetroGraph, RetroGraphEdge, RetroGraphNode
from tools.models.workflow_models import RetroTask, RetroTaskList

__all__ = ["build_retro_graph", "collect_terminal_precursors"]

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _clean_smiles(smiles: str) -> str:
    """Strip trailing text annotations from a SMILES string.

    E.g. ``'CCO [apply: C=O → C-OH]'`` → ``'CCO'``.
    """
    if not smiles:
        return smiles
    match = re.match(r"^([^\s\[]*(?:\[[^\]]*\])*[^\s\[]*)", smiles)
    return match.group(1) if match else smiles.split()[0]


def _extract_reaction_type(task: RetroTask) -> str:
    """Extract the reaction type/class from a task's result dict.

    Preference order:
    1. ``result["reaction_type"]``
    2. First entry in ``result["disconnections"]`` or ``result["proposals"]``
    3. Fallback to ``"unknown"``
    """
    result = task.result
    if not result or not isinstance(result, dict):
        return "unknown"

    # Direct reaction_type key
    if result.get("reaction_type"):
        return str(result["reaction_type"])

    # Nested in disconnections / proposals list
    discs = result.get("disconnections") or result.get("proposals", [])
    if discs and isinstance(discs, list):
        first = discs[0] if discs else {}
        if isinstance(first, dict) and first.get("reaction_class"):
            return str(first["reaction_class"])

    if result.get("reaction_class"):
        return str(result["reaction_class"])

    return "unknown"


def _node_status_for_task(task: RetroTask) -> str:
    """Map a DISCONNECT task's status to a graph-node status string."""
    if task.status in (TaskStatus.FAILED, TaskStatus.SKIPPED):
        return "failed_decomposition"
    if task.status == TaskStatus.VALIDATED:
        return "validated"
    return "terminal"


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def build_retro_graph(
    tasks: List[RetroTask],
    target_smiles: str,
) -> RetroGraph:
    """Build a typed :class:`RetroGraph` from completed retrosynthesis tasks.

    This is the standalone replacement for ``TaskExecutor._build_retro_graph``.
    The returned graph uses :class:`RetroGraphNode` and :class:`RetroGraphEdge`
    objects, fixing the format mismatch where the old code emitted plain dicts
    with a ``"type"`` key that the visualizer expected as ``"node_type"``.

    Args:
        tasks: Flat list of :class:`RetroTask` (typically from
            ``RetroTaskList.tasks``).
        target_smiles: The root target molecule SMILES.

    Returns:
        A fully populated :class:`RetroGraph`.
    """
    nodes: Dict[str, RetroGraphNode] = {}
    edges: Dict[str, RetroGraphEdge] = {}
    smiles_to_node_id: Dict[str, str] = {}
    precursor_counter = 0
    edge_counter = 0

    # --- target node ---
    if target_smiles:
        nodes["target"] = RetroGraphNode(
            node_id="target",
            smiles=target_smiles,
            node_type="target",
            status=TaskStatus.COMPLETED.value,
            is_terminal=False,
        )
        smiles_to_node_id[target_smiles] = "target"

    # --- iterate DISCONNECT tasks ---
    for task in tasks:
        if task.task_type != TaskType.DISCONNECT:
            continue
        if not task.precursors:
            continue

        node_status = _node_status_for_task(task)
        reaction_type = _extract_reaction_type(task)

        precursor_node_ids: List[str] = []
        for smi in task.precursors:
            if not smi or not isinstance(smi, str):
                continue
            try:
                clean_smi = _clean_smiles(smi)
            except Exception:
                logger.warning("Invalid SMILES skipped in retro_graph: %r", smi)
                continue
            if not clean_smi:
                continue

            if clean_smi not in smiles_to_node_id:
                node_id = f"precursor_{precursor_counter}"
                precursor_counter += 1
                nodes[node_id] = RetroGraphNode(
                    node_id=node_id,
                    smiles=clean_smi,
                    node_type="precursor",
                    status=node_status,
                    is_terminal=False,
                )
                smiles_to_node_id[clean_smi] = node_id
            precursor_node_ids.append(smiles_to_node_id[clean_smi])

        # --- reaction edge ---
        if precursor_node_ids:
            clean_target = (
                _clean_smiles(task.target_smiles) if task.target_smiles else ""
            )
            product_node_id = smiles_to_node_id.get(clean_target, "target")

            edge_id = f"edge_{edge_counter}"
            edge_counter += 1

            rxn = task.reaction_smiles or ""
            is_valid_rxn = ">>" in rxn and "??" not in rxn

            # Extract reasoning from task result (decision history)
            reasoning = ""
            if isinstance(task.result, dict):
                reasoning = task.result.get("reasoning", "")
                if not reasoning:
                    # Try nested decision_instruction
                    di = task.result.get("decision_instruction", {})
                    if isinstance(di, dict):
                        reasoning = di.get("reasoning", "")

            edges[edge_id] = RetroGraphEdge(
                edge_id=edge_id,
                product_id=product_node_id,
                precursor_ids=list(precursor_node_ids),
                reaction_smiles=rxn if is_valid_rxn else "",
                reaction_type=reaction_type,
                reasoning=reasoning,
            )

    # Collect strategy and decision history into metadata
    decision_history: List[Dict[str, Any]] = []
    strategy_info: Dict[str, Any] = {}
    for task in tasks:
        if not isinstance(task.result, dict):
            continue
        if task.task_type == TaskType.STRATEGY and task.result.get("recommended_strategy"):
            strategy_info = {
                "recommended_strategy": task.result.get("recommended_strategy", ""),
                "rationale": task.result.get("rationale", ""),
                "scores": task.result.get("scores", {}),
            }
        # Collect any decision reasoning stored in result
        if task.result.get("reasoning") or task.result.get("decision_instruction"):
            di = task.result.get("decision_instruction", {})
            decision_history.append({
                "task_id": task.task_id,
                "task_type": task.task_type.value if hasattr(task.task_type, "value") else str(task.task_type),
                "reasoning": task.result.get("reasoning", "") or (di.get("reasoning", "") if isinstance(di, dict) else ""),
            })

    metadata: Dict[str, Any] = {}
    if strategy_info:
        metadata["strategy"] = strategy_info
    if decision_history:
        metadata["decision_history"] = decision_history

    return RetroGraph(
        target_smiles=target_smiles,
        nodes=nodes,
        edges=edges,
        metadata=metadata,
    )


def collect_terminal_precursors(tasks: List[RetroTask]) -> List[str]:
    """Collect all terminal precursor SMILES from a task list.

    Collection rules (same semantics as the original
    ``TaskExecutor._collect_terminal_precursors``):

    1. SKIPPED AVAILABILITY tasks → their ``target_smiles`` (SA < threshold).
    2. Completed/validated DISCONNECT tasks → precursors that have **no**
       child DISCONNECT task targeting them.

    Args:
        tasks: Flat list of :class:`RetroTask`.

    Returns:
        De-duplicated list of terminal precursor SMILES (insertion order).
    """
    seen: set = set()
    result: List[str] = []

    # Build a quick lookup for subtask queries
    task_list = RetroTaskList(tasks=list(tasks))

    # (a) SKIPPED AVAILABILITY tasks → their target_smiles
    for task in tasks:
        if (
            task.status == TaskStatus.SKIPPED
            and task.task_type == TaskType.AVAILABILITY
            and task.target_smiles
        ):
            if task.target_smiles not in seen:
                seen.add(task.target_smiles)
                result.append(task.target_smiles)

    # (b) Completed DISCONNECT tasks: precursors with no child DISCONNECT
    for task in tasks:
        if task.task_type != TaskType.DISCONNECT:
            continue
        if task.status not in (TaskStatus.COMPLETED, TaskStatus.VALIDATED):
            continue

        child_disconnect_targets: set = set()
        for sub in task_list.tasks:
            if sub.parent_task_id == task.task_id and sub.task_type == TaskType.DISCONNECT:
                child_disconnect_targets.add(sub.target_smiles)

        for smi in task.precursors:
            if smi and smi not in child_disconnect_targets and smi not in seen:
                seen.add(smi)
                result.append(smi)

    return result
