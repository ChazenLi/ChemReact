"""
Unified Error Hierarchy
=======================
Exception-based error system for all MultiRetro tool modules.

Replaces the dataclass-based ToolError/ToolResult pattern from
internal/common/errors.py with a proper exception hierarchy.
"""

from enum import Enum
from typing import Any, Dict


class ErrorSeverity(Enum):
    """Error severity levels for the unified error system."""
    LOW = "low"
    MEDIUM = "medium"
    HIGH = "high"
    CRITICAL = "critical"


class RetroError(Exception):
    """Unified error base class for all MultiRetro errors.

    Attributes:
        code: Machine-readable error code (e.g. "INVALID_SMILES").
        message: Human-readable error description.
        severity: Error severity level.
    """

    def __init__(
        self,
        code: str,
        message: str,
        severity: ErrorSeverity = ErrorSeverity.MEDIUM,
    ):
        self.code = code
        self.message = message
        self.severity = severity
        super().__init__(message)

    def to_dict(self) -> Dict[str, Any]:
        """Serialize to a JSON-compatible dictionary."""
        return {
            "code": self.code,
            "message": self.message,
            "severity": self.severity.value,
        }


class ValidationError(RetroError):
    """Raised when reaction or input validation fails."""
    pass


class ExecutionError(RetroError):
    """Raised when a skill or tool execution fails."""
    pass


class SkillError(RetroError):
    """Raised when a skill encounters a logic error (bad args, unmet preconditions)."""
    pass


class ChemistryError(RetroError):
    """Raised when an RDKit or chemistry operation fails."""
    pass
