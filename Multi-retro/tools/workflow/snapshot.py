"""
Snapshot & persistence — save / load / checkpoint for retrosynthesis routes.

Extracted from ``core/taskflow/executor.py._save_snapshot``.  Provides
standalone functions that the orchestrator (or any caller) can use to
persist route state to disk as JSON snapshots and to restore from them
for pause-resume capability.

Snapshot layout on disk::

    {route_dir}/
        snapshots/
            after_task_{task_id}.json   — per-task snapshots
            checkpoint.json             — latest checkpoint for resume

Usage::

    from tools.workflow.snapshot import save_snapshot, load_snapshot, create_checkpoint, load_checkpoint

    save_snapshot(route, route_dir)
    route = load_snapshot(route_dir, task_id="1.2")
    create_checkpoint(route, route_dir)
    route = load_checkpoint(route_dir)
"""

from __future__ import annotations

import json
import logging
from pathlib import Path
from typing import List, Optional

from tools.common.status import TaskStatus
from tools.models.workflow_models import RetroRoute, RetroTask, RetroTaskList

__all__ = [
    "save_snapshot",
    "load_snapshot",
    "list_snapshots",
    "create_checkpoint",
    "load_checkpoint",
    "delete_snapshots",
]

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _snapshots_dir(route_dir: Path) -> Path:
    """Return the snapshots sub-directory for a given route directory."""
    return route_dir / "snapshots"


def _serialize_route(route: RetroRoute) -> str:
    """Serialize a *route* to a JSON string."""
    return json.dumps(route.to_dict(), ensure_ascii=False, indent=2, default=str)


def _deserialize_route(text: str) -> RetroRoute:
    """Deserialize a JSON string back into a ``RetroRoute``."""
    return RetroRoute.from_dict(json.loads(text))


# ---------------------------------------------------------------------------
# Save / load individual snapshots
# ---------------------------------------------------------------------------

def save_snapshot(
    route: RetroRoute,
    route_dir: Path,
    *,
    task_id: Optional[str] = None,
) -> Optional[Path]:
    """Persist *route* state as a snapshot JSON file.

    Parameters
    ----------
    route:
        The route whose state should be saved.
    route_dir:
        Directory for this route (``snapshots/`` will be created inside).
    task_id:
        Explicit task id to use in the filename.  When ``None`` the id of
        the most recently completed / failed / blocked task is used, falling
        back to ``"initial"``.

    Returns
    -------
    Path to the written snapshot file, or ``None`` if saving failed.
    """
    snap_dir = _snapshots_dir(route_dir)
    try:
        snap_dir.mkdir(parents=True, exist_ok=True)

        if task_id is None:
            terminal_statuses = {
                TaskStatus.COMPLETED,
                TaskStatus.FAILED,
                TaskStatus.BLOCKED,
            }
            completed = [
                t for t in route.tasks.tasks
                if t.status in terminal_statuses
            ]
            task_id = completed[-1].task_id if completed else "initial"

        snapshot_path = snap_dir / f"after_task_{task_id}.json"
        snapshot_path.write_text(_serialize_route(route), encoding="utf-8")
        logger.debug("Snapshot saved: %s", snapshot_path)
        return snapshot_path
    except Exception as exc:
        logger.warning("Failed to save snapshot: %s", exc)
        return None


def load_snapshot(
    route_dir: Path,
    task_id: str,
) -> Optional[RetroRoute]:
    """Load a previously saved snapshot for the given *task_id*.

    Returns ``None`` if the snapshot file does not exist or cannot be read.
    """
    snapshot_path = _snapshots_dir(route_dir) / f"after_task_{task_id}.json"
    if not snapshot_path.exists():
        logger.warning("Snapshot not found: %s", snapshot_path)
        return None
    try:
        text = snapshot_path.read_text(encoding="utf-8")
        route = _deserialize_route(text)
        logger.debug("Snapshot loaded: %s", snapshot_path)
        return route
    except Exception as exc:
        logger.warning("Failed to load snapshot %s: %s", snapshot_path, exc)
        return None


def list_snapshots(route_dir: Path) -> List[str]:
    """Return a sorted list of task ids for which snapshots exist.

    E.g. ``["initial", "1", "1.1", "2"]``.
    """
    snap_dir = _snapshots_dir(route_dir)
    if not snap_dir.exists():
        return []
    ids: List[str] = []
    for p in sorted(snap_dir.glob("after_task_*.json")):
        # filename: after_task_{task_id}.json
        stem = p.stem  # "after_task_1.2"
        prefix = "after_task_"
        if stem.startswith(prefix):
            ids.append(stem[len(prefix):])
    return ids


# ---------------------------------------------------------------------------
# Checkpoint (latest resumable state)
# ---------------------------------------------------------------------------

_CHECKPOINT_FILENAME = "checkpoint.json"


def create_checkpoint(
    route: RetroRoute,
    route_dir: Path,
) -> Optional[Path]:
    """Save a checkpoint — the latest resumable state for pause-resume.

    Unlike per-task snapshots, the checkpoint is always written to the
    same file (``checkpoint.json``) so the orchestrator can quickly
    resume without scanning for the latest snapshot.

    Returns the path to the checkpoint file, or ``None`` on failure.
    """
    snap_dir = _snapshots_dir(route_dir)
    try:
        snap_dir.mkdir(parents=True, exist_ok=True)
        checkpoint_path = snap_dir / _CHECKPOINT_FILENAME
        checkpoint_path.write_text(_serialize_route(route), encoding="utf-8")
        logger.debug("Checkpoint created: %s", checkpoint_path)
        return checkpoint_path
    except Exception as exc:
        logger.warning("Failed to create checkpoint: %s", exc)
        return None


def load_checkpoint(route_dir: Path) -> Optional[RetroRoute]:
    """Load the latest checkpoint for pause-resume.

    Returns ``None`` if no checkpoint exists or it cannot be read.
    """
    checkpoint_path = _snapshots_dir(route_dir) / _CHECKPOINT_FILENAME
    if not checkpoint_path.exists():
        logger.debug("No checkpoint found at %s", checkpoint_path)
        return None
    try:
        text = checkpoint_path.read_text(encoding="utf-8")
        route = _deserialize_route(text)
        logger.debug("Checkpoint loaded: %s", checkpoint_path)
        return route
    except Exception as exc:
        logger.warning("Failed to load checkpoint %s: %s", checkpoint_path, exc)
        return None


# ---------------------------------------------------------------------------
# Cleanup
# ---------------------------------------------------------------------------

def delete_snapshots(route_dir: Path) -> int:
    """Remove all snapshot files (including checkpoint) for a route.

    Returns the number of files deleted.
    """
    snap_dir = _snapshots_dir(route_dir)
    if not snap_dir.exists():
        return 0
    count = 0
    for p in snap_dir.iterdir():
        if p.is_file() and p.suffix == ".json":
            try:
                p.unlink()
                count += 1
            except OSError as exc:
                logger.warning("Failed to delete %s: %s", p, exc)
    return count
