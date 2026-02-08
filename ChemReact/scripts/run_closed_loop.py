"""CLI + compatibility entrypoint for ChemReact closed-loop pipeline."""

import sys
from pathlib import Path


_HERE = Path(__file__).resolve().parent
if str(_HERE) not in sys.path:
    sys.path.insert(0, str(_HERE))

# Backward compatibility: tests import helper symbols from this module.
from core_pipeline import *  # noqa: F401,F403
from core_pipeline import _build_repair_prompt, main  # noqa: F401
from planner_audit import _audit_single_step  # noqa: F401


if __name__ == "__main__":
    sys.exit(main())
