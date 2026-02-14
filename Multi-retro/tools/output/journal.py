"""
Process journal — append-only process recorder for retrosynthesis task execution.

Migrated from ``core/taskflow/journal.py``.  Writes structured step records
to ``process.md``, including tool call summaries, LLM thinking, conclusions,
and optional discrepancy warnings.

Key changes vs. the original:

* Uses ``TaskStatus`` / ``TaskType`` from ``tools.common.status``.
* Uses ``RetroTask`` / ``RetroRoute`` from ``tools.models.workflow_models``
  (the migrated models no longer carry ``title`` / ``route_name`` /
  ``strategy`` as top-level fields — those live in ``metadata``).
* Raises :class:`ExecutionError` from ``tools.common.errors`` when the
  journal path is invalid or a write fails.
* Follows the same module conventions as ``graph_builder.py`` and
  ``report_generator.py``.

Usage::

    from tools.output.journal import ProcessJournal

    journal = ProcessJournal(Path("outputs/process.md"))
    journal.write_header(route)
    journal.append_step(task, tool_call_summary="...", ...)
"""

from __future__ import annotations

import logging
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional

from tools.common.errors import ErrorSeverity, ExecutionError
from tools.common.status import TaskType
from tools.models.workflow_models import RetroRoute, RetroTask

__all__ = ["ProcessJournal"]

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Label map for reaction-condition keys → Chinese labels
# ---------------------------------------------------------------------------

_CONDITION_LABEL_MAP: Dict[str, str] = {
    "solvent": "溶剂",
    "base": "碱",
    "temperature": "温度",
    "catalyst": "催化剂",
    "reagents": "试剂",
    "reaction_type": "反应类型",
    "mechanism": "机理",
    "time": "时间",
    "notes": "备注",
}


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _format_conditions_block(conditions: Dict[str, Any]) -> List[str]:
    """Format a reaction-conditions dict into Markdown list lines.

    Empty / ``None`` values are silently skipped.  List values are
    joined with ``", "``.

    Returns:
        A list of Markdown lines (each ending with ``"\\n"``).
    """
    lines: List[str] = []
    for key, val in conditions.items():
        if val is None or val == "" or val == []:
            continue
        label = _CONDITION_LABEL_MAP.get(key, key)
        if isinstance(val, list):
            val = ", ".join(str(v) for v in val)
        lines.append(f"- **{label}**: {val}\n")
    return lines


def _build_thinking_text(
    llm_thinking: str,
    is_llm_driven: bool,
) -> str:
    """Determine the thinking text for a journal step.

    When *is_llm_driven* is ``True`` and *llm_thinking* is non-empty the
    thinking text is used as-is (LLM's real reasoning).  When
    *is_llm_driven* is ``True`` but *llm_thinking* is empty a
    ``"[规则驱动]"`` prefix is added.  Otherwise the raw text is returned.
    """
    if is_llm_driven and llm_thinking:
        return llm_thinking
    if is_llm_driven:
        return f"[规则驱动] {llm_thinking}" if llm_thinking else "[规则驱动]"
    return llm_thinking


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

class ProcessJournal:
    """Manages append-only writing to a ``process.md`` file.

    Args:
        process_md_path: Filesystem path for the journal file.

    Raises:
        ExecutionError: If the parent directory cannot be created or a
            write operation fails.
    """

    def __init__(self, process_md_path: Path | str) -> None:
        self.path = Path(process_md_path)

    # -- header ------------------------------------------------------------

    def write_header(self, route: RetroRoute) -> None:
        """Write the ``process.md`` header (called once at creation).

        Uses write mode (``'w'``) — only call this when first creating
        the file.  Route name and strategy are read from
        ``route.metadata`` with sensible fallbacks.
        """
        route_name = route.metadata.get("route_name", route.route_id)
        strategy = route.metadata.get("strategy", "unknown")
        timestamp = datetime.now().isoformat(timespec="seconds")

        header = (
            f"# 逆合成过程记录: {route_name}\n"
            f"\n"
            f"**Target**: `{route.target_smiles}`\n"
            f"**Strategy**: {strategy}\n"
            f"**Created**: {timestamp}\n"
            f"\n"
            f"---\n"
        )

        try:
            self.path.parent.mkdir(parents=True, exist_ok=True)
            self.path.write_text(header, encoding="utf-8")
        except OSError as exc:
            raise ExecutionError(
                code="JOURNAL_WRITE_FAILED",
                message=f"Failed to write journal header to {self.path}: {exc}",
                severity=ErrorSeverity.HIGH,
            ) from exc

        logger.info("Journal header written to %s", self.path)

    # -- append step -------------------------------------------------------

    def append_step(
        self,
        task: RetroTask,
        tool_call_summary: str,
        tool_result_summary: str,
        llm_thinking: str,
        conclusion: str,
        *,
        discrepancy: Optional[str] = None,
        reaction_smiles: Optional[str] = None,
        exploration_log: Optional[List[str]] = None,
        is_llm_driven: bool = False,
        reaction_conditions: Optional[Dict[str, Any]] = None,
    ) -> None:
        """Append a step record to ``process.md``.

        Opens in append mode (``'a'``) so existing content is never
        overwritten.  Generates a unique anchor ``#step-{task_id}`` for
        cross-referencing.

        Args:
            task: The :class:`RetroTask` this step belongs to.
            tool_call_summary: Human-readable summary of the tool call.
            tool_result_summary: Human-readable summary of the tool result.
            llm_thinking: LLM reasoning text (or empty string).
            conclusion: Conclusion text for this step.
            discrepancy: Optional discrepancy warning text.
            reaction_smiles: Optional reaction SMILES string.
            exploration_log: Optional list of exploration log entries.
            is_llm_driven: Whether this step was driven by LLM reasoning.
            reaction_conditions: Optional dict of reaction conditions.

        Raises:
            ExecutionError: If the append operation fails.
        """
        timestamp = datetime.now().isoformat(timespec="seconds")
        anchor = f"step-{task.task_id}"
        task_label = task.task_type.value
        thinking_text = _build_thinking_text(llm_thinking, is_llm_driven)

        lines: List[str] = [
            f"\n## <a id=\"{anchor}\"></a>Step {task.task_id}: [{task_label}]\n",
            f"\n**时间**: {timestamp}\n",
            f"\n### 🧪 思考\n{thinking_text}\n",
        ]

        # Exploration log subsection
        if exploration_log:
            lines.append("\n#### 探索过程\n")
            for entry in exploration_log:
                lines.append(f"- {entry}\n")

        lines.append(
            f"\n### 🔧 工具\n{tool_call_summary}\n"
            f"**结果**: {tool_result_summary}\n"
        )

        if discrepancy is not None:
            lines.append(f"\n### ⚠️ 差异警告\n{discrepancy}\n")

        if reaction_smiles:
            lines.append(f"\n### ⚗️ 反应式\n`{reaction_smiles}`\n")

        # Reaction conditions
        rc = reaction_conditions or {}
        if rc and isinstance(rc, dict):
            cond_lines = _format_conditions_block(rc)
            if cond_lines:
                lines.append("\n### 🧪 反应条件\n")
                lines.extend(cond_lines)

        lines.extend([
            f"\n### 📋 结论\n{conclusion}\n",
            "\n---\n",
        ])

        try:
            self.path.parent.mkdir(parents=True, exist_ok=True)
            with self.path.open("a", encoding="utf-8") as f:
                f.writelines(lines)
        except OSError as exc:
            raise ExecutionError(
                code="JOURNAL_APPEND_FAILED",
                message=f"Failed to append step to journal {self.path}: {exc}",
                severity=ErrorSeverity.HIGH,
            ) from exc

        logger.debug("Journal step %s appended to %s", task.task_id, self.path)
