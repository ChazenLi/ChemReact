"""
Render Report Skill
===================

Thin skill wrapper around :func:`tools.output.report_generator.generate_synthesis_report`.
Accepts a RetroGraph (or dict) and output directory, delegates to the tool module,
and returns a :class:`SkillResult` envelope.
"""

from __future__ import annotations

import logging
from typing import Any, Dict

from tools.skills.base import BaseSkill, SkillResult

logger = logging.getLogger(__name__)


class RenderReportSkill(BaseSkill):
    """Generate synthesis report from completed RetroGraph.

    Delegates to :func:`tools.output.report_generator.generate_synthesis_report`.
    """

    name = "render_report"
    description = "Generate synthesis report from completed RetroGraph"

    def execute(self, args: Any) -> Dict[str, Any]:
        """Run report generation and return a SkillResult dict.

        Args:
            args: dict with keys:
                - retro_graph (dict or RetroGraph): Required. The graph to render.
                - output_dir (str): Optional. Defaults to ``"outputs/report"``.

        Returns:
            ``SkillResult.to_dict()`` with report path and content under ``data``.
        """
        if not isinstance(args, dict) or "retro_graph" not in args:
            # Support namespace/dataclass args from skill_dispatch
            if hasattr(args, "retro_graph"):
                args = {
                    "retro_graph": args.retro_graph,
                    "output_dir": getattr(args, "output_dir", "outputs/report"),
                }
            else:
                return SkillResult(
                    success=False,
                    error="args must be a dict with 'retro_graph'",
                ).to_dict()

        output_dir = args.get("output_dir", "outputs/report") or "outputs/report"

        try:
            from tools.models.output_models import RetroGraph
            from tools.output.report_generator import generate_synthesis_report

            graph_input = args["retro_graph"]
            if isinstance(graph_input, dict):
                graph = RetroGraph.from_dict(graph_input)
            elif isinstance(graph_input, RetroGraph):
                graph = graph_input
            else:
                return SkillResult(
                    success=False,
                    error="retro_graph must be a dict or RetroGraph instance",
                ).to_dict()

            report_content = generate_synthesis_report(graph, output_dir)
        except Exception as exc:
            logger.exception("render_report failed")
            return SkillResult(success=False, error=str(exc)).to_dict()

        return SkillResult(
            success=True,
            data={
                "output_dir": str(output_dir),
                "report_content": report_content,
            },
        ).to_dict()
