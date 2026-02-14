"""
Decision Manager — manages decision points in the retrosynthesis workflow.

Integrates the original ``core/taskflow/decision_manager.py`` and
``core/taskflow/decision_models.py`` into a single module.

Key capabilities:
* Detects 5 decision point types (strategy, disconnection, validation,
  repair, recursion).
* Builds ``DecisionContext`` for the LLM host with full context.
* Applies ``DecisionInstruction`` from the LLM (pause-resume protocol).
* Provides rule-driven fallback decisions for backward compatibility.
* Persists decision history and pending decision state to disk.
* Detects discrepancies between LLM decisions and Skill recommendations.

All data classes (``DecisionContext``, ``DecisionInstruction``, ``RetroTask``,
etc.) are imported from ``tools.models.workflow_models``.  Status enums come
from ``tools.common.status``.
"""

from __future__ import annotations

import json
import logging
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional

from tools.common.status import TaskStatus, TaskType
from tools.common.constants import SAThresholds
from tools.models.workflow_models import (
    DecisionContext,
    DecisionInstruction,
    RetroTask,
    RetroRoute,
)

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Decision type constants (the 5 decision point types)
# ---------------------------------------------------------------------------

STRATEGY_SELECTION = "strategy_selection"
DISCONNECTION_DECISION = "disconnection_decision"
VALIDATION_JUDGMENT = "validation_judgment"
REPAIR_JUDGMENT = "repair_judgment"
RECURSION_DECISION = "recursion_decision"

# Map TaskType → decision type string for automatic detection
_TASK_TYPE_TO_DECISION: Dict[TaskType, str] = {
    TaskType.STRATEGY: STRATEGY_SELECTION,
    TaskType.ANALYZE: DISCONNECTION_DECISION,
    TaskType.VALIDATE: VALIDATION_JUDGMENT,
    TaskType.REPAIR: REPAIR_JUDGMENT,
}


class DecisionManager:
    """Manages decision points in the retrosynthesis workflow.

    Supports 5 decision point types and the pause-resume protocol:
    1. Strategy selection — which synthesis strategy to use
    2. Disconnection decision — which bond to break
    3. Validation judgment — accept/repair/retry after reaction validation
    4. Repair judgment — accept/reject automatic repair
    5. Recursion decision — how to handle precursor expansion

    Parameters
    ----------
    route : RetroRoute
        The current synthesis route being executed.
    route_dir : Path
        Directory for persisting decision history and pending state.
    """

    def __init__(self, route: RetroRoute, route_dir: Path) -> None:
        self.route = route
        self.route_dir = route_dir
        self._history_path = route_dir / "decision_history.json"
        self._pending_path = route_dir / "pending_decision.json"
        self._history: List[Dict[str, Any]] = self._load_history()

    # ------------------------------------------------------------------
    # Decision history
    # ------------------------------------------------------------------

    @property
    def decision_history(self) -> List[Dict[str, Any]]:
        """Return a copy of the current decision history."""
        return list(self._history)

    def _load_history(self) -> List[Dict[str, Any]]:
        """Load decision history from disk."""
        if self._history_path.exists():
            try:
                data = json.loads(self._history_path.read_text(encoding="utf-8"))
                return data.get("decisions", [])
            except Exception as exc:
                logger.warning("Failed to load decision history: %s", exc)
        return []

    def _save_history(self) -> None:
        """Persist decision history to disk."""
        self._history_path.parent.mkdir(parents=True, exist_ok=True)
        self._history_path.write_text(
            json.dumps(
                {"decisions": self._history},
                ensure_ascii=False,
                indent=2,
                default=str,
            ),
            encoding="utf-8",
        )

    def record_decision(
        self,
        instruction: DecisionInstruction,
        decision_type: str,
        *,
        was_fallback: bool = False,
    ) -> None:
        """Record a completed decision to history.

        Parameters
        ----------
        instruction : DecisionInstruction
            The decision that was applied.
        decision_type : str
            One of the 5 decision type constants.
        was_fallback : bool
            Whether this was a rule-driven fallback decision.
        """
        entry = {
            "decision_type": decision_type,
            "task_id": instruction.task_id,
            "action": instruction.action,
            "params": instruction.params,
            "reasoning": instruction.reasoning or "",
            "exploration_log": instruction.exploration_log or [],
            "timestamp": datetime.now().isoformat(),
            "was_fallback": was_fallback,
        }
        self._history.append(entry)
        self._save_history()

    def _build_history_summary(self) -> List[Dict[str, Any]]:
        """Build a compact decision history summary for DecisionContext."""
        summary = []
        for h in self._history:
            reasoning = h.get("reasoning", "")
            summary.append({
                "decision_type": h["decision_type"],
                "task_id": h["task_id"],
                "action": h["action"],
                "reasoning_preview": reasoning[:200] if reasoning else "",
            })
        return summary

    # ------------------------------------------------------------------
    # Decision point detection
    # ------------------------------------------------------------------

    def check_decision_point(
        self, task: RetroTask, result: Dict[str, Any]
    ) -> Optional[DecisionContext]:
        """Check if the completed *task* triggers a decision point.

        Returns a ``DecisionContext`` if a decision is needed, ``None``
        otherwise.
        """
        decision_type = _TASK_TYPE_TO_DECISION.get(task.task_type)
        if decision_type is None:
            return None

        builders = {
            STRATEGY_SELECTION: self._build_strategy_context,
            DISCONNECTION_DECISION: self._build_disconnection_context,
            VALIDATION_JUDGMENT: self._build_validation_context,
            REPAIR_JUDGMENT: self._build_repair_context,
        }
        builder = builders.get(decision_type)
        if builder is None:
            return None
        return builder(task, result)

    def check_recursion_decision_point(
        self, task: RetroTask, precursors: List[str]
    ) -> Optional[DecisionContext]:
        """Check if precursor expansion triggers a recursion decision."""
        if not precursors:
            return None
        return self._build_recursion_context(task, precursors)

    # ------------------------------------------------------------------
    # Context builders
    # ------------------------------------------------------------------

    def _route_summary(self) -> Dict[str, Any]:
        """Build a compact route summary for inclusion in DecisionContext."""
        completed = sum(
            1
            for t in self.route.tasks.tasks
            if t.status in (TaskStatus.COMPLETED, TaskStatus.VALIDATED)
        )
        max_depth = max((t.depth for t in self.route.tasks.tasks), default=0)
        return {
            "route_id": self.route.route_id,
            "target_smiles": self.route.target_smiles,
            "completed_steps": completed,
            "total_tasks": len(self.route.tasks.tasks),
            "current_max_depth": max_depth,
        }

    def _build_strategy_context(
        self, task: RetroTask, result: Dict[str, Any]
    ) -> DecisionContext:
        """Build context for strategy selection decision."""
        scores = result.get("scores", {})
        context = {
            "target_smiles": task.target_smiles or self.route.target_smiles,
            "skill_result": result,
            "recommended_strategy": result.get("recommended_strategy"),
            "rationale": result.get("rationale"),
            "scores": scores,
            "convergent_feasibility": scores.get("convergent_feasibility"),
            "linear_feasibility": scores.get("linear_feasibility"),
            "scaffold_summary": result.get("scaffold_summary", {}),
            "selectivity_summary": result.get("selectivity_summary", {}),
            "route_summary": self._route_summary(),
        }
        actions = [
            {
                "action_id": "choose_strategy",
                "description": "选择合成策略",
                "params": ["chosen_strategy"],
            },
            {"action_id": "use_default", "description": "使用规则驱动默认选择"},
        ]
        return DecisionContext(
            decision_type=STRATEGY_SELECTION,
            task_id=task.task_id,
            context=context,
            available_actions=actions,
            decision_history=self._build_history_summary(),
            exploration_tools=_get_exploration_tools(STRATEGY_SELECTION),
        )

    def _build_disconnection_context(
        self, task: RetroTask, result: Dict[str, Any]
    ) -> DecisionContext:
        """Build context for disconnection decision."""
        atom_bond_map = result.get("atom_bond_map", {})
        strategic_bonds = atom_bond_map.get("strategic_bonds", [])
        context = {
            "target_smiles": task.target_smiles or self.route.target_smiles,
            "skill_result": result,
            "atom_bond_map": atom_bond_map,
            "strategic_bonds": strategic_bonds,
            "retro_guidance": result.get("retro_guidance", ""),
            "functional_groups": result.get("functional_groups", []),
            "summary": result.get("summary", {}),
            "route_summary": self._route_summary(),
        }
        actions = [
            {
                "action_id": "select_bond",
                "description": "选择断裂位点",
                "params": ["atom1_idx", "atom2_idx"],
            },
            {
                "action_id": "custom_reaction",
                "description": "提供自定义反应式",
                "params": ["reaction_smiles", "precursors"],
            },
            {"action_id": "use_default", "description": "使用规则驱动默认选择"},
        ]
        return DecisionContext(
            decision_type=DISCONNECTION_DECISION,
            task_id=task.task_id,
            context=context,
            available_actions=actions,
            decision_history=self._build_history_summary(),
            exploration_tools=_get_exploration_tools(DISCONNECTION_DECISION),
        )

    def _build_validation_context(
        self, task: RetroTask, result: Dict[str, Any]
    ) -> DecisionContext:
        """Build context for validation judgment decision."""
        context = {
            "target_smiles": task.target_smiles or self.route.target_smiles,
            "skill_result": result,
            "is_valid": result.get("is_valid", False),
            "verdict": result.get("verdict", ""),
            "issues": result.get("issues", []),
            "atom_balance": result.get("atom_balance", {}),
            "reaction_smiles": task.reaction_smiles or "",
            "route_summary": self._route_summary(),
            "feasibility_assessment": result.get("feasibility_assessment", {}),
        }
        actions = [
            {"action_id": "accept", "description": "接受验证结果"},
            {"action_id": "repair", "description": "创建修复子任务"},
            {
                "action_id": "retry_other_bond",
                "description": "放弃当前断键，回到 DISCONNECT 重新选择",
            },
            {
                "action_id": "re_analyze_disconnection",
                "description": "回到 ANALYZE 重新分析断裂方案",
            },
            {
                "action_id": "accept_with_note",
                "description": "带备注接受 (附加 reaction_conditions)",
            },
            {"action_id": "use_default", "description": "使用规则驱动默认选择"},
        ]
        return DecisionContext(
            decision_type=VALIDATION_JUDGMENT,
            task_id=task.task_id,
            context=context,
            available_actions=actions,
            decision_history=self._build_history_summary(),
            exploration_tools=_get_exploration_tools(VALIDATION_JUDGMENT),
        )

    def _build_repair_context(
        self, task: RetroTask, result: Dict[str, Any]
    ) -> DecisionContext:
        """Build context for repair judgment decision."""
        context = {
            "target_smiles": task.target_smiles or self.route.target_smiles,
            "skill_result": result,
            "repaired_smiles": result.get("repaired_smiles", ""),
            "repaired_precursors": result.get("repaired_precursors", []),
            "original_issues": result.get("issues", []),
            "repair_strategy": result.get("repair_strategy", ""),
            "changes_made": result.get("changes_made", ""),
            "original_reaction_smiles": task.reaction_smiles or "",
            "route_summary": self._route_summary(),
        }
        actions = [
            {"action_id": "accept_repair", "description": "接受自动修复结果"},
            {
                "action_id": "use_custom_repair",
                "description": "使用自定义修复",
                "params": ["repaired_reaction_smiles", "repaired_precursors"],
            },
            {"action_id": "reject_and_retry", "description": "拒绝修复，重试"},
            {"action_id": "abandon_bond", "description": "放弃该断裂路径"},
            {"action_id": "use_default", "description": "使用规则驱动默认选择"},
        ]
        return DecisionContext(
            decision_type=REPAIR_JUDGMENT,
            task_id=task.task_id,
            context=context,
            available_actions=actions,
            decision_history=self._build_history_summary(),
            exploration_tools=_get_exploration_tools(REPAIR_JUDGMENT),
        )

    def _build_recursion_context(
        self, task: RetroTask, precursors: List[str]
    ) -> DecisionContext:
        """Build context for recursion decision (precursor expansion)."""
        precursor_info = []
        for smi in precursors:
            sa = _estimate_sa_score(smi)
            mw = _estimate_mw(smi)
            precursor_info.append({
                "smiles": smi,
                "sa_score": round(sa, 2),
                "molecular_weight": round(mw, 1),
                "suggested_action": (
                    "mark_terminal"
                    if mw < 120 or sa < SAThresholds.PURCHASABLE_MAX
                    else "mark_optional"
                    if sa < SAThresholds.COMPLEX_MIN
                    else "decompose"
                ),
            })
        context = {
            "target_smiles": task.target_smiles or self.route.target_smiles,
            "precursors": precursor_info,
            "current_depth": task.depth,
            "route_summary": self._route_summary(),
        }
        actions = [
            {
                "action_id": "recursion_decision",
                "description": "为每个前体指定操作",
                "params": ["precursor_decisions"],
            },
            {"action_id": "use_default", "description": "使用规则驱动默认选择"},
        ]
        return DecisionContext(
            decision_type=RECURSION_DECISION,
            task_id=task.task_id,
            context=context,
            available_actions=actions,
            decision_history=self._build_history_summary(),
            exploration_tools=_get_exploration_tools(RECURSION_DECISION),
        )


    # ------------------------------------------------------------------
    # Decision application (pause-resume protocol)
    # ------------------------------------------------------------------

    def apply_decision(
        self,
        instruction: DecisionInstruction,
        task: RetroTask,
        pending_context: Optional[DecisionContext] = None,
    ) -> Dict[str, Any]:
        """Apply a decision instruction to the task.

        Parameters
        ----------
        instruction : DecisionInstruction
            The LLM's chosen action.
        task : RetroTask
            The task in AWAITING_DECISION state.
        pending_context : DecisionContext, optional
            The original context that was presented to the LLM.

        Returns
        -------
        dict
            Result of applying the decision, including any discrepancy info.

        Raises
        ------
        ValueError
            If the task is not in AWAITING_DECISION state or the decision
            type / action is unknown.
        """
        if task.status != TaskStatus.AWAITING_DECISION:
            raise ValueError(
                f"Task {task.task_id} is not in AWAITING_DECISION state "
                f"(current: {task.status.value})"
            )

        # Determine decision type
        decision_type = (
            pending_context.decision_type
            if pending_context
            else _TASK_TYPE_TO_DECISION.get(task.task_type)
        )
        if decision_type is None:
            raise ValueError(
                f"Cannot determine decision type for task {task.task_id}"
            )

        # Handle "use_default" — delegate to fallback
        if instruction.action == "use_default":
            fallback = self.get_fallback_decision(
                decision_type, task, task.result or {}
            )
            if instruction.reasoning:
                fallback.reasoning = instruction.reasoning
            instruction = fallback

        # Dispatch to type-specific handler
        handlers = {
            STRATEGY_SELECTION: self._apply_strategy_decision,
            DISCONNECTION_DECISION: self._apply_disconnection_decision,
            VALIDATION_JUDGMENT: self._apply_validation_decision,
            REPAIR_JUDGMENT: self._apply_repair_decision,
            RECURSION_DECISION: self._apply_recursion_decision,
        }
        handler = handlers.get(decision_type)
        if handler is None:
            raise ValueError(
                f"No handler for decision type {decision_type!r}"
            )

        result = handler(instruction, task)

        # Detect discrepancy between LLM and Skill recommendation
        discrepancy = self._detect_discrepancy(
            instruction, task, pending_context
        )
        if discrepancy:
            result["discrepancy"] = discrepancy

        # Record and clean up
        self.record_decision(instruction, decision_type)
        self._remove_pending()

        return result

    # -- Type-specific apply handlers --

    def _apply_strategy_decision(
        self, instruction: DecisionInstruction, task: RetroTask
    ) -> Dict[str, Any]:
        chosen = instruction.params.get("chosen_strategy", "linear")
        if task.depth == 0:
            self.route.metadata["strategy"] = chosen
        task.status = TaskStatus.COMPLETED
        return {"applied": "strategy_selection", "chosen_strategy": chosen}

    def _apply_disconnection_decision(
        self, instruction: DecisionInstruction, task: RetroTask
    ) -> Dict[str, Any]:
        if instruction.action == "select_bond":
            a1 = instruction.params.get("atom1_idx")
            a2 = instruction.params.get("atom2_idx")
            if a1 is not None and a2 is not None:
                task.selected_bond = {"atom1_idx": a1, "atom2_idx": a2}
            task.status = TaskStatus.COMPLETED
            return {
                "applied": "select_bond",
                "selected_bond": task.selected_bond,
            }
        elif instruction.action == "custom_reaction":
            rxn = instruction.params.get("reaction_smiles", "")
            precursors = instruction.params.get("precursors", [])
            task.reaction_smiles = rxn
            task.precursors = precursors
            task.status = TaskStatus.COMPLETED
            return {
                "applied": "custom_reaction",
                "reaction_smiles": rxn,
                "precursors": precursors,
            }
        else:
            raise ValueError(
                f"Unknown disconnection action: {instruction.action}"
            )

    def _apply_validation_decision(
        self, instruction: DecisionInstruction, task: RetroTask
    ) -> Dict[str, Any]:
        action = instruction.action
        if action in ("accept", "accept_with_note"):
            task.status = TaskStatus.COMPLETED
            return {"applied": action, "validation_status": "PASS"}
        elif action == "repair":
            task.status = TaskStatus.COMPLETED
            return {"applied": "repair", "needs_repair": True}
        elif action == "retry_other_bond":
            task.status = TaskStatus.FAILED
            return {"applied": "retry_other_bond", "reset_disconnect": True}
        elif action == "re_analyze_disconnection":
            task.status = TaskStatus.FAILED
            return {
                "applied": "re_analyze_disconnection",
                "reset_analyze": True,
            }
        else:
            raise ValueError(f"Unknown validation action: {action}")

    def _apply_repair_decision(
        self, instruction: DecisionInstruction, task: RetroTask
    ) -> Dict[str, Any]:
        action = instruction.action
        if action == "accept_repair":
            task.status = TaskStatus.COMPLETED
            return {"applied": "accept_repair"}
        elif action == "use_custom_repair":
            rxn = instruction.params.get("repaired_reaction_smiles", "")
            precursors = instruction.params.get("repaired_precursors", [])
            task.reaction_smiles = rxn
            task.precursors = precursors
            task.status = TaskStatus.COMPLETED
            return {"applied": "use_custom_repair", "reaction_smiles": rxn}
        elif action == "reject_and_retry":
            task.status = TaskStatus.FAILED
            return {"applied": "reject_and_retry"}
        elif action == "abandon_bond":
            task.status = TaskStatus.FAILED
            return {"applied": "abandon_bond"}
        else:
            raise ValueError(f"Unknown repair action: {action}")

    def _apply_recursion_decision(
        self, instruction: DecisionInstruction, task: RetroTask
    ) -> Dict[str, Any]:
        decisions = instruction.params.get("precursor_decisions", [])
        task.status = TaskStatus.COMPLETED
        return {
            "applied": "recursion_decision",
            "precursor_decisions": decisions,
        }

    # ------------------------------------------------------------------
    # Fallback (rule-driven) decisions
    # ------------------------------------------------------------------

    def get_fallback_decision(
        self, decision_type: str, task: RetroTask, result: Dict[str, Any]
    ) -> DecisionInstruction:
        """Generate a rule-driven fallback decision.

        Replicates the existing automatic decision logic so the system
        can operate without LLM intervention when needed.
        """
        fallbacks = {
            STRATEGY_SELECTION: self._fallback_strategy,
            DISCONNECTION_DECISION: self._fallback_disconnection,
            VALIDATION_JUDGMENT: self._fallback_validation,
            REPAIR_JUDGMENT: self._fallback_repair,
            RECURSION_DECISION: self._fallback_recursion,
        }
        fb = fallbacks.get(decision_type)
        if fb is None:
            logger.warning(
                "No fallback handler for decision type %s on task %s; "
                "marking task completed",
                decision_type,
                task.task_id,
            )
            return DecisionInstruction(
                task_id=task.task_id,
                action="accept",
                reasoning="[规则驱动] 无对应决策处理器，自动接受",
            )
        return fb(task, result)

    def _fallback_strategy(
        self, task: RetroTask, result: Dict[str, Any]
    ) -> DecisionInstruction:
        strategy = result.get("recommended_strategy", "linear")
        return DecisionInstruction(
            task_id=task.task_id,
            action="choose_strategy",
            params={"chosen_strategy": strategy},
            reasoning=f"[规则驱动] 使用 Skill 推荐策略: {strategy}",
        )

    def _fallback_disconnection(
        self, task: RetroTask, result: Dict[str, Any]
    ) -> DecisionInstruction:
        atom_bond_map = result.get("atom_bond_map", {})
        bonds = atom_bond_map.get("strategic_bonds", [])
        if bonds:
            sorted_bonds = sorted(
                bonds,
                key=lambda b: b.get("strategic_priority", 0),
                reverse=True,
            )
            best = sorted_bonds[0]
            return DecisionInstruction(
                task_id=task.task_id,
                action="select_bond",
                params={
                    "atom1_idx": best.get("atom1_idx"),
                    "atom2_idx": best.get("atom2_idx"),
                },
                reasoning=(
                    f"[规则驱动] 按优先级选择第一个战略键: "
                    f"{best.get('chemical_label', 'N/A')}"
                ),
            )
        return DecisionInstruction(
            task_id=task.task_id,
            action="select_bond",
            params={"atom1_idx": None, "atom2_idx": None},
            reasoning="[规则驱动] 无战略键信息，由 ProposeDisconnectionSkill 自动选择",
        )

    def _fallback_validation(
        self, task: RetroTask, result: Dict[str, Any]
    ) -> DecisionInstruction:
        is_valid = result.get("is_valid", False)
        if is_valid:
            return DecisionInstruction(
                task_id=task.task_id,
                action="accept",
                reasoning="[规则驱动] 验证通过",
            )
        return DecisionInstruction(
            task_id=task.task_id,
            action="repair",
            reasoning="[规则驱动] 验证失败，创建修复任务",
        )

    def _fallback_repair(
        self, task: RetroTask, result: Dict[str, Any]
    ) -> DecisionInstruction:
        if result.get("repaired_smiles"):
            return DecisionInstruction(
                task_id=task.task_id,
                action="accept_repair",
                reasoning="[规则驱动] 接受自动修复结果",
            )
        return DecisionInstruction(
            task_id=task.task_id,
            action="abandon_bond",
            reasoning="[规则驱动] 修复失败，放弃该路径",
        )

    def _fallback_recursion(
        self, task: RetroTask, result: Dict[str, Any]
    ) -> DecisionInstruction:
        precursor_decisions = []
        for smi in task.precursors or []:
            sa = _estimate_sa_score(smi)
            mw = _estimate_mw(smi)
            if mw < 120 or sa < SAThresholds.PURCHASABLE_MAX:
                action = "mark_terminal"
                reason = f"SA={sa:.2f}, MW={mw:.1f}; terminal"
            elif sa < SAThresholds.COMPLEX_MIN:
                action = "mark_optional"
                reason = f"SA={sa:.2f}, MW={mw:.1f}; optional"
            else:
                action = "decompose"
                reason = f"SA={sa:.2f}, MW={mw:.1f}; requires decomposition"
            precursor_decisions.append(
                {"smiles": smi, "action": action, "reason": reason}
            )
        return DecisionInstruction(
            task_id=task.task_id,
            action="recursion_decision",
            params={"precursor_decisions": precursor_decisions},
            reasoning="[规则驱动] 基于 SA_Score 阈值自动决策",
        )

    # ------------------------------------------------------------------
    # Pending decision persistence (pause-resume)
    # ------------------------------------------------------------------

    def save_pending_decision(
        self,
        context: DecisionContext,
        task: RetroTask,
        executor_state: Optional[Dict[str, Any]] = None,
    ) -> None:
        """Persist pending decision state to ``pending_decision.json``."""
        data = {
            "decision_context": context.to_dict(),
            "task_snapshot": task.to_dict(),
            "executor_state": executor_state or {},
            "timestamp": datetime.now().isoformat(),
        }
        self._pending_path.parent.mkdir(parents=True, exist_ok=True)
        self._pending_path.write_text(
            json.dumps(data, ensure_ascii=False, indent=2, default=str),
            encoding="utf-8",
        )

    def load_pending_decision(self) -> Optional[Dict[str, Any]]:
        """Load pending decision from ``pending_decision.json``.

        Returns dict with ``decision_context``, ``task_snapshot``,
        ``executor_state``, ``timestamp``, or ``None`` if no pending
        decision exists.
        """
        if not self._pending_path.exists():
            return None
        try:
            data = json.loads(self._pending_path.read_text(encoding="utf-8"))
            data["decision_context"] = DecisionContext.from_dict(
                data["decision_context"]
            )
            return data
        except Exception as exc:
            logger.warning("Failed to load pending decision: %s", exc)
            return None

    def _remove_pending(self) -> None:
        """Remove ``pending_decision.json`` after decision is applied."""
        try:
            if self._pending_path.exists():
                self._pending_path.unlink()
        except Exception as exc:
            logger.warning("Failed to remove pending decision file: %s", exc)

    # ------------------------------------------------------------------
    # Discrepancy detection
    # ------------------------------------------------------------------

    def _detect_discrepancy(
        self,
        instruction: DecisionInstruction,
        task: RetroTask,
        pending_context: Optional[DecisionContext] = None,
    ) -> Optional[str]:
        """Detect if LLM decision differs from Skill recommendation."""
        if pending_context is None:
            return None

        ctx = pending_context.context
        dt = pending_context.decision_type

        if dt == STRATEGY_SELECTION:
            recommended = ctx.get("recommended_strategy")
            chosen = instruction.params.get("chosen_strategy")
            if recommended and chosen and recommended != chosen:
                return (
                    f"LLM 选择策略 '{chosen}' 与 Skill 推荐 "
                    f"'{recommended}' 不同"
                )

        elif dt == VALIDATION_JUDGMENT:
            is_valid = ctx.get("is_valid", False)
            if is_valid and instruction.action == "repair":
                return "Skill 判定验证通过，但 LLM 选择修复"
            elif not is_valid and instruction.action in (
                "accept",
                "accept_with_note",
            ):
                return "Skill 判定验证失败，但 LLM 选择接受"

        elif dt == REPAIR_JUDGMENT:
            has_repair = bool(ctx.get("repaired_smiles"))
            if has_repair and instruction.action == "abandon_bond":
                return "Skill 提供了修复结果，但 LLM 选择放弃"
            elif not has_repair and instruction.action == "accept_repair":
                return "Skill 未提供修复结果，但 LLM 选择接受修复"

        return None


# ---------------------------------------------------------------------------
# Module-level helpers
# ---------------------------------------------------------------------------

def _get_exploration_tools(decision_type: str) -> List[Dict[str, Any]]:
    """Return exploration tool definitions applicable to *decision_type*.

    This is a lightweight inline version.  When ``tools.workflow.exploration``
    is available, it should be used instead.
    """
    try:
        from tools.workflow.exploration import get_tools_for_decision  # type: ignore[import-not-found]
        return get_tools_for_decision(decision_type)
    except ImportError:
        pass

    # Fallback: minimal built-in tool list
    _ALL_TOOLS: List[Dict[str, Any]] = [
        {
            "tool_id": "assess_feasibility",
            "description": "正向合成可行性评估",
            "applicable_decisions": [
                DISCONNECTION_DECISION,
                VALIDATION_JUDGMENT,
                REPAIR_JUDGMENT,
            ],
        },
        {
            "tool_id": "analyze_selectivity",
            "description": "化学选择性分析",
            "applicable_decisions": [
                STRATEGY_SELECTION,
                DISCONNECTION_DECISION,
                VALIDATION_JUDGMENT,
                RECURSION_DECISION,
            ],
        },
        {
            "tool_id": "analyze_molecule",
            "description": "分子分析",
            "applicable_decisions": [DISCONNECTION_DECISION, RECURSION_DECISION],
        },
        {
            "tool_id": "validate_reaction",
            "description": "反应验证",
            "applicable_decisions": [DISCONNECTION_DECISION, REPAIR_JUDGMENT],
        },
        {
            "tool_id": "check_availability",
            "description": "前体可用性检查",
            "applicable_decisions": [RECURSION_DECISION, VALIDATION_JUDGMENT],
        },
    ]
    return [
        t for t in _ALL_TOOLS if decision_type in t["applicable_decisions"]
    ]


def _estimate_sa_score(smiles: str) -> float:
    """Estimate SA Score — delegates to ``tools.chem.sa_scorer``."""
    try:
        from tools.chem.sa_scorer import estimate_sa_score
        return estimate_sa_score(smiles)
    except Exception:
        return 3.0


def _estimate_mw(smiles: str) -> float:
    """Estimate molecular weight from SMILES."""
    try:
        from tools.chem.mol_parser import parse_smiles
        mol = parse_smiles(smiles)
        if mol is not None:
            from rdkit.Chem import Descriptors
            return Descriptors.MolWt(mol)
    except Exception:
        pass
    # Rough heuristic fallback
    heavy = sum(1 for c in smiles if c.isalpha() and c.isupper())
    return heavy * 14.0
