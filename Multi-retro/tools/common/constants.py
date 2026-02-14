"""
Centralized Configuration Registry
====================================

Single source of truth for all thresholds, limits, and shared constants.
Replaces hardcoded values previously scattered across executor.py,
availability.py, and internal/common/constants.py.

Usage::

    from tools.common.constants import SAThresholds, RetroLimits

    category = SAThresholds.classify(sa_score)
    if depth >= RetroLimits.MAX_DEPTH:
        ...
"""


class SAThresholds:
    """SA Score classification thresholds (global single definition).

    Unifies the inconsistent thresholds from executor.py (4.0) and
    availability.py (5.0) into one authoritative source.

    Classification:
        - purchasable:           SA < 2.2  (terminal node)
        - easily_synthesizable:  SA < 3.5  (optional further decomposition)
        - complex:               SA >= 3.5 (must decompose further)
    """

    PURCHASABLE_MAX: float = 2.2
    EASILY_SYNTHESIZABLE_MAX: float = 3.5
    COMPLEX_MIN: float = 3.5

    @classmethod
    def classify(cls, sa_score: float) -> str:
        """Classify a molecule by its SA score.

        Returns one of ``"purchasable"``, ``"easily_synthesizable"``,
        or ``"complex"``.
        """
        if sa_score < cls.PURCHASABLE_MAX:
            return "purchasable"
        if sa_score < cls.EASILY_SYNTHESIZABLE_MAX:
            return "easily_synthesizable"
        return "complex"


class RetroLimits:
    """Retrosynthesis recursion and resource limits."""

    MAX_DEPTH: int = 7                    # max recursion depth for decomposition
    MW_TERMINAL: float = 120.0            # molecular weight (Da) below which a precursor is terminal
    MAX_REPAIR_RETRIES: int = 3           # max repair attempts per reaction
    MAX_EXPLORATION_CALLS: int = 5        # max exploration calls per decision point
    MAX_GRAPH_DEPTH: int = 20             # max recursion depth for visualization / graph traversal
    MAX_TOTAL_TASKS: int = 50             # max tasks in a single route


# ---------------------------------------------------------------------------
# Migrated from internal/common/constants.py
# ---------------------------------------------------------------------------

CONFIDENCE_THRESHOLD: float = 0.3
"""Classifier confidence threshold — below this returns 'unknown'."""

DEFAULT_MAX_BYPRODUCTS: int = 3
"""Maximum byproducts to infer during reaction validation."""

DEFAULT_MAX_PROPOSALS: int = 5
"""Default max disconnection proposals."""
