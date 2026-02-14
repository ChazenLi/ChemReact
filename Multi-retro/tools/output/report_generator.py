"""
Synthesis report generator — integrated from ``core/taskflow/report_builder.py``
and ``core/skills/render.py``.

Generates a ``SYNTHESIS_REPORT.md`` from a typed :class:`RetroGraph`.  The
report includes four sections:

1. Target molecule
2. Synthesis steps (reaction edges with conditions)
3. Terminal precursors
4. Reaction conditions summary

Usage::

    from tools.output.report_generator import generate_synthesis_report

    path = generate_synthesis_report(retro_graph, Path("outputs/report"))
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any, Dict, List

from tools.models.output_models import RetroGraph, RetroGraphEdge, RetroGraphNode

__all__ = ["generate_synthesis_report"]

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _collect_terminal_nodes(graph: RetroGraph) -> List[RetroGraphNode]:
    """Return nodes that are terminal precursors (leaf nodes)."""
    # A node is terminal if it is flagged is_terminal, or if it is a
    # precursor that never appears as a product_id in any edge.
    product_ids = {edge.product_id for edge in graph.edges.values()}
    terminals: List[RetroGraphNode] = []
    for node in graph.nodes.values():
        if node.is_terminal:
            terminals.append(node)
        elif node.node_type == "precursor" and node.node_id not in product_ids:
            terminals.append(node)
    return terminals


def _format_conditions(conditions: Dict[str, Any] | None) -> str:
    """Format reaction conditions dict into a readable string."""
    if not conditions:
        return "N/A"
    parts: List[str] = []
    for key, value in conditions.items():
        parts.append(f"{key}: {value}")
    return "; ".join(parts) if parts else "N/A"


def _build_step_number_map(graph: RetroGraph) -> Dict[str, int]:
    """Assign a 1-based step number to each edge, ordered by edge_id."""
    sorted_edge_ids = sorted(graph.edges.keys())
    return {eid: idx + 1 for idx, eid in enumerate(sorted_edge_ids)}


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def generate_synthesis_report(
    retro_graph: RetroGraph,
    output_dir: Path | str,
) -> str:
    """Generate a ``SYNTHESIS_REPORT.md`` from a :class:`RetroGraph`.

    The report contains:

    1. **Target Molecule** — the root SMILES being synthesised.
    2. **Synthesis Steps** — each reaction edge with product, precursors,
       reaction SMILES, reaction type, and conditions.
    3. **Terminal Precursors** — leaf-node molecules (purchasable starting
       materials).
    4. **Reaction Conditions** — a summary table of conditions per step.

    Args:
        retro_graph: A fully populated :class:`RetroGraph`.
        output_dir: Directory where ``SYNTHESIS_REPORT.md`` will be written.
            Created automatically if it does not exist.

    Returns:
        The full Markdown content of the generated report.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    lines: List[str] = []
    step_map = _build_step_number_map(retro_graph)
    terminals = _collect_terminal_nodes(retro_graph)

    # --- Header ---
    lines.append("# Synthesis Report")
    lines.append("")

    # --- 1. Target Molecule ---
    lines.append("## Target Molecule")
    lines.append("")
    lines.append(f"**SMILES:** `{retro_graph.target_smiles}`")
    lines.append("")

    # --- 2. Synthesis Steps ---
    lines.append("## Synthesis Steps")
    lines.append("")
    if not retro_graph.edges:
        lines.append("No synthesis steps recorded.")
        lines.append("")
    else:
        lines.append("| Step | Product | Precursors | Reaction SMILES | Reaction Type |")
        lines.append("|------|---------|------------|-----------------|---------------|")
        for edge_id in sorted(retro_graph.edges.keys()):
            edge = retro_graph.edges[edge_id]
            step_num = step_map[edge_id]

            product_node = retro_graph.nodes.get(edge.product_id)
            product_smi = product_node.smiles if product_node else edge.product_id

            precursor_smiles_list: List[str] = []
            for pid in edge.precursor_ids:
                pnode = retro_graph.nodes.get(pid)
                precursor_smiles_list.append(pnode.smiles if pnode else pid)
            precursors_str = ", ".join(precursor_smiles_list) if precursor_smiles_list else "N/A"

            rxn_smi = edge.reaction_smiles or "N/A"
            rxn_type = edge.reaction_type or "unknown"

            lines.append(
                f"| {step_num} | `{product_smi}` | `{precursors_str}` | `{rxn_smi}` | {rxn_type} |"
            )
        lines.append("")

    # --- 3. Terminal Precursors ---
    lines.append("## Terminal Precursors")
    lines.append("")
    if not terminals:
        lines.append("No terminal precursors identified.")
        lines.append("")
    else:
        lines.append("| # | SMILES | SA Score | Status |")
        lines.append("|---|--------|----------|--------|")
        for idx, node in enumerate(terminals, 1):
            sa_str = f"{node.sa_score:.2f}" if node.sa_score is not None else "N/A"
            lines.append(f"| {idx} | `{node.smiles}` | {sa_str} | {node.status} |")
        lines.append("")

    # --- 4. Reaction Conditions ---
    lines.append("## Reaction Conditions")
    lines.append("")
    edges_with_conditions = [
        (eid, e) for eid, e in retro_graph.edges.items() if e.conditions
    ]
    if not edges_with_conditions:
        lines.append("No reaction conditions recorded.")
        lines.append("")
    else:
        lines.append("| Step | Conditions |")
        lines.append("|------|------------|")
        for edge_id, edge in sorted(edges_with_conditions, key=lambda x: x[0]):
            step_num = step_map[edge_id]
            cond_str = _format_conditions(edge.conditions)
            lines.append(f"| {step_num} | {cond_str} |")
        lines.append("")

    # --- 5. Strategy & Decision Rationale ---
    has_rationale = False

    # Global strategy from metadata
    strategy = retro_graph.metadata.get("strategy", {})
    if isinstance(strategy, dict) and strategy.get("recommended_strategy"):
        has_rationale = True
        lines.append("## Strategy & Decision Rationale")
        lines.append("")
        lines.append(f"**Strategy:** {strategy['recommended_strategy']}")
        if strategy.get("rationale"):
            lines.append(f"**Rationale:** {strategy['rationale']}")
        scores = strategy.get("scores", {})
        if scores:
            parts = [f"{k}: {v}" for k, v in scores.items()]
            lines.append(f"**Scores:** {'; '.join(parts)}")
        lines.append("")

    # Per-step reasoning from edges
    edges_with_reasoning = [
        (eid, e) for eid, e in retro_graph.edges.items() if e.reasoning
    ]
    if edges_with_reasoning:
        if not has_rationale:
            lines.append("## Strategy & Decision Rationale")
            lines.append("")
        lines.append("### Step Reasoning")
        lines.append("")
        lines.append("| Step | Reaction Type | Reasoning |")
        lines.append("|------|---------------|-----------|")
        for edge_id, edge in sorted(edges_with_reasoning, key=lambda x: x[0]):
            step_num = step_map[edge_id]
            rxn_type = edge.reaction_type or "unknown"
            lines.append(f"| {step_num} | {rxn_type} | {edge.reasoning} |")
        lines.append("")

    # Decision history from metadata
    decision_history = retro_graph.metadata.get("decision_history", [])
    if decision_history and isinstance(decision_history, list):
        if not has_rationale and not edges_with_reasoning:
            lines.append("## Strategy & Decision Rationale")
            lines.append("")
        lines.append("### Decision History")
        lines.append("")
        lines.append("| Task | Type | Reasoning |")
        lines.append("|------|------|-----------|")
        for entry in decision_history:
            if isinstance(entry, dict):
                tid = entry.get("task_id", "")
                ttype = entry.get("task_type", "")
                reason = entry.get("reasoning", "")
                if reason:
                    lines.append(f"| {tid} | {ttype} | {reason} |")
        lines.append("")

    report_content = "\n".join(lines)

    report_path = output_dir / "SYNTHESIS_REPORT.md"
    report_path.write_text(report_content, encoding="utf-8")
    logger.info("Synthesis report written to %s", report_path)

    # --- 5. Synthesis Tree Visualization ---
    try:
        from tools.output.visualizer import draw_synthesis_tree

        tree_path = output_dir / "synthesis_tree.png"
        drawn = draw_synthesis_tree(retro_graph, tree_path)
        if drawn:
            logger.info("Synthesis tree image written to %s", tree_path)
        else:
            logger.info("Synthesis tree image skipped (PIL or data unavailable)")
    except Exception as exc:
        logger.warning("Synthesis tree visualization failed: %s", exc)

    return report_content
