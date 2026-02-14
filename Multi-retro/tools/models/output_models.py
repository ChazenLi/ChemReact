"""
Output data models — retrosynthesis graph representation.

Provides typed dataclasses for the retrosynthesis graph (nodes, edges,
and the full graph).  Every class implements ``to_dict()`` /
``from_dict()`` with round-trip consistency and missing-field tolerance.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional


# ---------------------------------------------------------------------------
# RetroGraphNode
# ---------------------------------------------------------------------------

@dataclass
class RetroGraphNode:
    """A node in the retrosynthesis graph."""

    node_id: str = ""
    smiles: str = ""
    node_type: str = ""  # "target", "precursor", "intermediate"
    status: str = ""
    sa_score: Optional[float] = None
    is_terminal: bool = False

    def to_dict(self) -> Dict[str, Any]:
        return {
            "node_id": self.node_id,
            "smiles": self.smiles,
            "node_type": self.node_type,
            "status": self.status,
            "sa_score": self.sa_score,
            "is_terminal": self.is_terminal,
        }

    @classmethod
    def from_dict(cls, d: Dict[str, Any]) -> "RetroGraphNode":
        return cls(
            node_id=d.get("node_id", ""),
            smiles=d.get("smiles", ""),
            node_type=d.get("node_type", ""),
            status=d.get("status", ""),
            sa_score=d.get("sa_score"),
            is_terminal=d.get("is_terminal", False),
        )


# ---------------------------------------------------------------------------
# RetroGraphEdge
# ---------------------------------------------------------------------------

@dataclass
class RetroGraphEdge:
    """An edge in the retrosynthesis graph (reaction step)."""

    edge_id: str = ""
    product_id: str = ""
    precursor_ids: List[str] = field(default_factory=list)
    reaction_smiles: str = ""
    reaction_type: str = ""
    conditions: Optional[Dict[str, Any]] = None
    reasoning: str = ""

    def to_dict(self) -> Dict[str, Any]:
        d: Dict[str, Any] = {
            "edge_id": self.edge_id,
            "product_id": self.product_id,
            "precursor_ids": list(self.precursor_ids),
            "reaction_smiles": self.reaction_smiles,
            "reaction_type": self.reaction_type,
            "conditions": self.conditions,
        }
        if self.reasoning:
            d["reasoning"] = self.reasoning
        return d

    @classmethod
    def from_dict(cls, d: Dict[str, Any]) -> "RetroGraphEdge":
        return cls(
            edge_id=d.get("edge_id", ""),
            product_id=d.get("product_id", ""),
            precursor_ids=list(d.get("precursor_ids", [])),
            reaction_smiles=d.get("reaction_smiles", ""),
            reaction_type=d.get("reaction_type", ""),
            conditions=d.get("conditions"),
            reasoning=d.get("reasoning", ""),
        )


# ---------------------------------------------------------------------------
# RetroGraph
# ---------------------------------------------------------------------------

@dataclass
class RetroGraph:
    """Complete retrosynthesis graph."""

    target_smiles: str = ""
    nodes: Dict[str, RetroGraphNode] = field(default_factory=dict)
    edges: Dict[str, RetroGraphEdge] = field(default_factory=dict)
    metadata: Dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> Dict[str, Any]:
        return {
            "target_smiles": self.target_smiles,
            "nodes": {k: v.to_dict() for k, v in self.nodes.items()},
            "edges": {k: v.to_dict() for k, v in self.edges.items()},
            "metadata": dict(self.metadata),
        }

    @classmethod
    def from_dict(cls, d: Dict[str, Any]) -> "RetroGraph":
        return cls(
            target_smiles=d.get("target_smiles", ""),
            nodes={
                k: RetroGraphNode.from_dict(v)
                for k, v in d.get("nodes", {}).items()
            },
            edges={
                k: RetroGraphEdge.from_dict(v)
                for k, v in d.get("edges", {}).items()
            },
            metadata=dict(d.get("metadata", {})),
        )
