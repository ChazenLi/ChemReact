"""
Simplified Skill Base Class
============================

Abstract base for all MultiRetro skills.  Removes the Mesh-specific
``can_run`` / ``apply`` three-phase interface from the legacy
``core.skills.base.BaseSkill`` and keeps only the pure
``execute(args) -> dict`` contract.

Usage::

    from tools.skills.base import BaseSkill, SkillResult

    class MySkill(BaseSkill):
        name = "my_skill"
        description = "Does something useful."

        def execute(self, args):
            return SkillResult(success=True, data={"key": "value"}).to_dict()
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from typing import Any, Dict


@dataclass
class SkillResult:
    """Standardised result envelope returned by every skill.

    Attributes:
        success: Whether the skill execution succeeded.
        data: Arbitrary JSON-serializable payload.
        error: Human-readable error string (empty on success).
    """

    success: bool
    data: Dict[str, Any] = field(default_factory=dict)
    error: str = ""

    def to_dict(self) -> Dict[str, Any]:
        """Serialize to a plain dictionary."""
        return {
            "success": self.success,
            "data": self.data,
            "error": self.error,
        }

    @classmethod
    def from_dict(cls, d: Dict[str, Any]) -> "SkillResult":
        """Deserialize from a plain dictionary."""
        return cls(
            success=d["success"],
            data=d.get("data", {}),
            error=d.get("error", ""),
        )


class BaseSkill(ABC):
    """Simplified Skill base class — no Mesh-specific interfaces.

    Subclasses must set ``name`` and ``description`` and implement
    :meth:`execute`.
    """

    name: str = ""
    description: str = ""

    @abstractmethod
    def execute(self, args: Any) -> Dict[str, Any]:
        """Execute the skill and return a JSON-serializable dict.

        Args:
            args: Skill-specific arguments (type depends on the skill).

        Returns:
            A JSON-serializable dictionary with the skill's output.
            Prefer returning ``SkillResult.to_dict()`` for consistency.

        Raises:
            NotImplementedError: If the subclass forgets to implement this.
        """
        raise NotImplementedError
