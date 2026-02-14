"""
Guidance tools for MultiRetro.

Provides strategic analysis and decision-making tools:

- strategy_advisor:        Comprehensive strategy decision (linear/convergent/mixed)
- selectivity_analyzer:    Chemoselectivity analysis and protection requirements
- feasibility_assessor:    Forward synthesis feasibility assessment
- disconnection_proposer:  Retrosynthetic disconnection proposals (35+ rules)
- repair_advisor:          Repair guidance based on validation issue codes
- retro_guide:             Reaction-aware retro analysis guide builder
"""

from .strategy_advisor import decide_strategy, analyze_scaffold, plan_ring_strategy
from .selectivity_analyzer import analyze_selectivity
from .feasibility_assessor import assess_forward_feasibility
from .disconnection_proposer import propose_disconnections
from .repair_advisor import generate_repair_guidance, build_repair_prompt
from .retro_guide import build_retro_guide, build_fragment_profiles
