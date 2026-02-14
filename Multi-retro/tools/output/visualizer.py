"""
Visualizer — synthesis tree drawing and molecule image generation.

Migrated from ``core/taskflow/visualizer.py`` and ``core/graph/visualizer.py``.

Bug fixes applied:
- ``_fallback_hierarchical_labels()``: Added MAX_DEPTH recursion limit
  (``RetroLimits.MAX_GRAPH_DEPTH``) to prevent stack overflow on cyclic or
  deeply nested graphs.  (Requirement 10.3)
- ``_add_bromine_to_aryl()``: Added atom index validation before accessing
  atoms to prevent ``IndexError`` / ``RuntimeError``.  (Requirement 10.1)

Public API:
    draw_synthesis_tree(retro_graph, output_path)
    generate_molecule_image(smiles, output_path)
"""

from __future__ import annotations

import logging
import os
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

from tools.common.constants import RetroLimits
from tools.models.output_models import RetroGraph

logger = logging.getLogger(__name__)

try:
    from PIL import Image, ImageDraw
    PIL_AVAILABLE = True
except Exception:
    PIL_AVAILABLE = False

try:
    from rdkit import Chem
    from rdkit.Chem import Draw
    RDKIT_AVAILABLE = True
except Exception:
    RDKIT_AVAILABLE = False


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------


def _label_parts(label: str) -> List[int]:
    """Parse a hierarchical label like ``"1.2.3"`` into sortable ints."""
    out: List[int] = []
    for part in str(label).split("."):
        try:
            out.append(int(part))
        except Exception:
            out.append(9999)
    return out


def _build_precursor_children(edges: Dict[str, Any]) -> Dict[str, List[str]]:
    """Build parent → children mapping from edges.

    Supports two edge formats:
      - Standard: ``{"product_id": "X", "precursor_ids": ["A", "B"]}``
      - Simplified: ``{"parent": "X", "children": ["A", "B"]}``
    """
    children: Dict[str, List[str]] = {}
    for _, edge in sorted(edges.items(), key=lambda kv: str(kv[0])):
        product_id = edge.get("product_id") or edge.get("parent")
        precursor_ids = edge.get("precursor_ids") or edge.get("children", [])
        if not isinstance(product_id, str):
            continue
        if not isinstance(precursor_ids, list):
            continue
        children.setdefault(product_id, [])
        for pid in precursor_ids:
            if isinstance(pid, str) and pid not in children[product_id]:
                children[product_id].append(pid)
    return children


def _fallback_hierarchical_labels(
    nodes: Dict[str, Any],
    edges: Dict[str, Any],
) -> Dict[str, str]:
    """Generate hierarchical labels (``"1"``, ``"1.1"``, …) via DFS.

    **Bug fix (Requirement 10.3):** The original implementation had no
    recursion depth limit, causing ``RecursionError`` on cyclic or very
    deep graphs.  We now cap traversal at ``RetroLimits.MAX_GRAPH_DEPTH``
    (default 20).
    """
    max_depth: int = RetroLimits.MAX_GRAPH_DEPTH

    # Find the target / root node.
    target_id: Optional[str] = None
    for nid, ndata in nodes.items():
        if ndata.get("is_target") or ndata.get("type") == "target":
            target_id = nid
            break
    if target_id is None and nodes:
        target_id = sorted(nodes.keys())[0]
    if target_id is None:
        return {}

    children = _build_precursor_children(edges)
    labels: Dict[str, str] = {target_id: "1"}
    visited: set = set()

    def _dfs(nid: str, depth: int = 0) -> None:
        if nid in visited:
            return
        if depth >= max_depth:
            return
        visited.add(nid)
        base = labels[nid]
        ordered = children.get(nid, [])
        for i, child in enumerate(ordered, start=1):
            labels.setdefault(child, f"{base}.{i}")
            _dfs(child, depth + 1)

    _dfs(target_id, 0)
    return labels


# ---------------------------------------------------------------------------
# Bromine helper (migrated from reaction_repairer.py with bug fix)
# ---------------------------------------------------------------------------


def _get_aromatic_carbon_idx(mol: Any) -> Optional[int]:
    """Find first aromatic carbon with an available H for functionalization."""
    for atom in mol.GetAtoms():
        if atom.GetIsAromatic() and atom.GetSymbol() == "C":
            if atom.GetTotalNumHs() > 0:
                return atom.GetIdx()
    return None


def _has_smarts(mol: Any, smarts: str) -> bool:
    """Check whether *mol* contains the given SMARTS pattern."""
    if not RDKIT_AVAILABLE:
        return False
    pat = Chem.MolFromSmarts(smarts)
    if pat is None:
        return False
    return mol.HasSubstructMatch(pat)


def _add_bromine_to_aryl(smiles: str) -> Optional[str]:
    """Add bromine to an aromatic ring.

    **Bug fix (Requirement 10.1):** The original code did not validate
    the atom index returned by ``_get_aromatic_carbon_idx`` against the
    molecule's actual atom count, which could cause ``IndexError`` or
    ``RuntimeError`` when the index was stale or out of range.  We now
    validate the index both before and after creating the editable
    molecule.
    """
    if not RDKIT_AVAILABLE:
        return None

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None

    # Already halogenated — nothing to do.
    if _has_smarts(mol, "[c][Br,I,Cl]"):
        return smiles

    # Find aromatic carbon and validate index.
    idx = _get_aromatic_carbon_idx(mol)
    if idx is None or idx >= mol.GetNumAtoms():
        return None

    # Build editable molecule and re-validate index.
    emol = Chem.RWMol(mol)
    if idx >= emol.GetNumAtoms():
        return None

    br_idx = emol.AddAtom(Chem.Atom(35))  # Br
    emol.AddBond(idx, br_idx, Chem.BondType.SINGLE)

    try:
        Chem.SanitizeMol(emol)
        return Chem.MolToSmiles(emol, canonical=True)
    except Exception:
        return None


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def draw_synthesis_tree(
    retro_graph: RetroGraph,
    output_path: Path,
) -> bool:
    """Draw a synthesis tree image from a *RetroGraph*.

    Renders each node as a panel (with molecule image or SMILES text)
    and connects them with orthogonal lines following the graph edges.

    Args:
        retro_graph: Typed ``RetroGraph`` instance.
        output_path: Destination file path for the PNG image.

    Returns:
        ``True`` on success, ``False`` on failure (logged as warning).
    """
    if not PIL_AVAILABLE:
        logger.warning("PIL not available; cannot draw synthesis tree")
        return False

    graph_dict = retro_graph.to_dict()
    nodes = graph_dict.get("nodes", {})
    edges = graph_dict.get("edges", {})

    if not nodes:
        logger.warning("Empty graph — nothing to draw")
        return False

    # Resolve hierarchical labels.
    hierarchical_labels: Dict[str, str] = (
        graph_dict.get("metadata", {}).get("hierarchical_labels", {})
    )
    if not isinstance(hierarchical_labels, dict) or not hierarchical_labels:
        hierarchical_labels = _fallback_hierarchical_labels(nodes, edges)
    if not hierarchical_labels:
        logger.warning("Could not compute hierarchical labels")
        return False

    # Build parent→children map and BFS levels.
    children_map = _build_precursor_children(edges)

    all_children_set: set = set()
    all_parents_set: set = set()
    for p, cs in children_map.items():
        all_parents_set.add(p)
        all_children_set.update(cs)
    root_nids = [
        nid for nid in hierarchical_labels
        if nid in all_parents_set and nid not in all_children_set
    ]
    if not root_nids:
        root_nids = (
            ["target"] if "target" in hierarchical_labels
            else [next(iter(hierarchical_labels))]
        )

    node_depth: Dict[str, int] = {}
    queue = list(root_nids)
    for r in root_nids:
        node_depth[r] = 0
    while queue:
        nid = queue.pop(0)
        d = node_depth[nid]
        for child in children_map.get(nid, []):
            if child not in node_depth:
                node_depth[child] = d + 1
                queue.append(child)
    for nid in hierarchical_labels:
        if nid not in node_depth:
            node_depth[nid] = len(str(hierarchical_labels[nid]).split("."))

    # Adaptive panel sizing.
    total_nodes = len(hierarchical_labels)
    if total_nodes <= 8:
        panel_w, panel_h = 280, 220
    elif total_nodes <= 16:
        panel_w, panel_h = 240, 200
    else:
        panel_w, panel_h = 200, 170

    # Render each node panel.
    rendered: Dict[str, Image.Image] = {}
    for node_id in hierarchical_labels:
        try:
            panel = Image.new("RGB", (panel_w, panel_h), (255, 255, 255))
            draw_p = ImageDraw.Draw(panel)
            lbl = hierarchical_labels.get(node_id, node_id)
            draw_p.text((8, 6), lbl, fill=(0, 0, 0))

            smiles = str(nodes.get(node_id, {}).get("smiles", ""))
            mol_rendered = False
            if smiles and RDKIT_AVAILABLE:
                try:
                    mol = Chem.MolFromSmiles(smiles)
                    if mol is not None:
                        mol_img = Draw.MolToImage(
                            mol, size=(panel_w - 16, panel_h - 36)
                        )
                        mol_img = mol_img.convert("RGB")
                        x = (panel_w - mol_img.width) // 2
                        y = 24 + max(0, (panel_h - 32 - mol_img.height) // 2)
                        panel.paste(mol_img, (x, y))
                        mol_rendered = True
                except Exception:
                    pass
            if not mol_rendered:
                draw_p.text(
                    (8, panel_h // 2 - 8),
                    smiles[:32] if smiles else node_id,
                    fill=(60, 60, 60),
                )
            draw_p.rectangle(
                [0, 0, panel_w - 1, panel_h - 1],
                outline=(180, 180, 180),
                width=1,
            )
            rendered[node_id] = panel
        except Exception:
            continue

    if not rendered:
        return False

    # Group by BFS level.
    level_nodes: Dict[int, List[str]] = {}
    for nid in hierarchical_labels:
        if nid not in rendered:
            continue
        lv = node_depth.get(nid, 0)
        level_nodes.setdefault(lv, []).append(nid)
    for lv in level_nodes:
        level_nodes[lv] = sorted(
            level_nodes[lv], key=lambda nid: _label_parts(hierarchical_labels[nid])
        )
    levels = sorted(level_nodes.keys())
    if not levels:
        return False

    # Canvas layout.
    x_gap = max(20, 50 - total_nodes)
    y_gap = 60
    margin = 30
    max_cols = max(len(level_nodes[lv]) for lv in levels)
    canvas_w = margin * 2 + max_cols * panel_w + (max_cols - 1) * x_gap
    canvas_h = margin * 2 + len(levels) * panel_h + (len(levels) - 1) * y_gap
    canvas = Image.new("RGB", (canvas_w, canvas_h), (255, 255, 255))
    draw = ImageDraw.Draw(canvas)

    node_pos: Dict[str, Tuple[int, int]] = {}
    for row_idx, lv in enumerate(levels):
        row = level_nodes[lv]
        row_w = len(row) * panel_w + (len(row) - 1) * x_gap
        start_x = (canvas_w - row_w) // 2
        y = margin + row_idx * (panel_h + y_gap)
        for i, nid in enumerate(row):
            x = start_x + i * (panel_w + x_gap)
            node_pos[nid] = (x, y)

    # Draw connector lines.
    line_color = (140, 140, 140)
    line_w = 2
    for parent_nid, child_nids in children_map.items():
        if parent_nid not in node_pos:
            continue
        px, py = node_pos[parent_nid]
        p_bot = (px + panel_w // 2, py + panel_h)
        for child_nid in child_nids:
            if child_nid not in node_pos:
                continue
            cx, cy = node_pos[child_nid]
            c_top = (cx + panel_w // 2, cy)
            route_y = py + panel_h + y_gap // 2
            draw.line([p_bot, (p_bot[0], route_y)], fill=line_color, width=line_w)
            draw.line(
                [(p_bot[0], route_y), (c_top[0], route_y)],
                fill=line_color,
                width=line_w,
            )
            draw.line([(c_top[0], route_y), c_top], fill=line_color, width=line_w)

    # Paste panels on top of lines.
    for nid, (x, y) in node_pos.items():
        if nid in rendered:
            canvas.paste(rendered[nid], (x, y))

    # Save.
    out = Path(output_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    canvas.save(str(out))
    return True


def generate_molecule_image(
    smiles: str,
    output_path: Path,
    size: Tuple[int, int] = (300, 300),
    legend: str = "",
) -> bool:
    """Generate a PNG image of a molecule from its SMILES.

    Args:
        smiles: SMILES string.
        output_path: Destination file path for the PNG.
        size: Image dimensions ``(width, height)``.
        legend: Optional legend text below the molecule.

    Returns:
        ``True`` on success, ``False`` on failure.
    """
    if not RDKIT_AVAILABLE:
        logger.warning("RDKit not available; cannot generate molecule image")
        return False

    if not smiles or not smiles.strip():
        logger.warning("Empty SMILES provided")
        return False

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        logger.warning("Invalid SMILES: %s", smiles)
        return False

    try:
        img = Draw.MolToImage(mol, size=size, legend=legend)
        out = Path(output_path)
        out.parent.mkdir(parents=True, exist_ok=True)
        img.save(str(out))
        return True
    except Exception as exc:
        logger.warning("Failed to generate molecule image for %s: %s", smiles, exc)
        return False
