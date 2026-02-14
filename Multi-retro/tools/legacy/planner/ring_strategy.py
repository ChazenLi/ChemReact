"""
MultiRetro - Ring Formation Strategy
=====================================
环形成策略分析器：决定成环顺序、开环断裂、成环反应选择。
"""

from __future__ import annotations
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Tuple
import json


# ─────────────────────────────────────────────────────────────────────────────
# Ring-Forming Reactions Database
# ─────────────────────────────────────────────────────────────────────────────

RING_FORMING_REACTIONS = {
    # Carbocycle formation
    "diels_alder": {
        "name": "Diels-Alder",
        "ring_size": 6,
        "bond_type": "C-C",
        "retro_bonds": 2,
        "conditions": "thermal or Lewis acid",
        "stereo_control": "high",
    },
    "aldol_cyclization": {
        "name": "Aldol Cyclization",
        "ring_size": [5, 6],
        "bond_type": "C-C",
        "retro_bonds": 1,
        "conditions": "base or acid",
        "stereo_control": "moderate",
    },
    "rcm": {
        "name": "Ring-Closing Metathesis",
        "ring_size": [5, 6, 7, 8],
        "bond_type": "C=C",
        "retro_bonds": 1,
        "conditions": "Grubbs catalyst",
        "stereo_control": "moderate",
    },
    "robinson_annulation": {
        "name": "Robinson Annulation",
        "ring_size": 6,
        "bond_type": "C-C",
        "retro_bonds": 2,
        "conditions": "base, heat",
        "stereo_control": "moderate",
    },
    "nazarov": {
        "name": "Nazarov Cyclization",
        "ring_size": 5,
        "bond_type": "C-C",
        "retro_bonds": 1,
        "conditions": "Lewis acid",
        "stereo_control": "good",
    },
    "prins": {
        "name": "Prins Cyclization",
        "ring_size": 6,
        "bond_type": "C-C/C-O",
        "retro_bonds": 2,
        "conditions": "acid",
        "stereo_control": "good",
    },
    "friedel_crafts": {
        "name": "Friedel-Crafts Cyclization",
        "ring_size": [5, 6],
        "bond_type": "C-C",
        "retro_bonds": 1,
        "conditions": "Lewis acid",
        "stereo_control": "low",
    },
    "pictet_spengler": {
        "name": "Pictet-Spengler",
        "ring_size": 6,
        "bond_type": "C-N/C-C",
        "retro_bonds": 2,
        "conditions": "acid",
        "stereo_control": "moderate",
    },
    # Heterocycle formation
    "fischer_indole": {
        "name": "Fischer Indole Synthesis",
        "ring_size": 5,
        "bond_type": "C-N/C-C",
        "retro_bonds": 2,
        "conditions": "acid, heat",
        "stereo_control": "N/A",
    },
    "click_triazole": {
        "name": "Click Chemistry (Triazole)",
        "ring_size": 5,
        "bond_type": "N-N/C-N",
        "retro_bonds": 2,
        "conditions": "Cu(I) catalyst",
        "stereo_control": "regioselective",
    },
    "hantzsch_pyridine": {
        "name": "Hantzsch Pyridine",
        "ring_size": 6,
        "bond_type": "C-N/C-C",
        "retro_bonds": 3,
        "conditions": "heat",
        "stereo_control": "N/A",
    },
}


# ─────────────────────────────────────────────────────────────────────────────
# Data Classes
# ─────────────────────────────────────────────────────────────────────────────

@dataclass
class RingConstruction:
    """单个环的构建计划"""
    ring_id: int
    ring_size: int
    construction_method: str
    step_order: int  # 何时构建
    is_aromatic: bool
    difficulty: str  # easy, moderate, hard
    dependencies: List[int] = field(default_factory=list)  # 依赖的其他环
    notes: str = ""
    
    def to_dict(self) -> Dict[str, Any]:
        return {
            "ring_id": self.ring_id,
            "ring_size": self.ring_size,
            "construction_method": self.construction_method,
            "step_order": self.step_order,
            "is_aromatic": self.is_aromatic,
            "difficulty": self.difficulty,
            "dependencies": self.dependencies,
            "notes": self.notes,
        }


@dataclass
class RingOpeningDisconnection:
    """开环断裂提案"""
    ring_id: int
    bond_atoms: Tuple[int, int]
    opening_reaction: str
    priority: str
    fragments_description: str
    retro_reaction: str
    
    def to_dict(self) -> Dict[str, Any]:
        return {
            "ring_id": self.ring_id,
            "bond_atoms": list(self.bond_atoms),
            "opening_reaction": self.opening_reaction,
            "priority": self.priority,
            "fragments_description": self.fragments_description,
            "retro_reaction": self.retro_reaction,
        }


@dataclass
class RingPlan:
    """完整环策略"""
    total_rings: int
    construction_order: List[RingConstruction]
    ring_opening_disconnections: List[RingOpeningDisconnection]
    ring_forming_reactions_used: List[str]
    early_stage_rings: List[int]
    late_stage_rings: List[int]
    timing_rationale: str
    
    def to_dict(self) -> Dict[str, Any]:
        return {
            "total_rings": self.total_rings,
            "construction_order": [c.to_dict() for c in self.construction_order],
            "ring_opening_disconnections": [d.to_dict() for d in self.ring_opening_disconnections],
            "ring_forming_reactions_used": self.ring_forming_reactions_used,
            "early_stage_rings": self.early_stage_rings,
            "late_stage_rings": self.late_stage_rings,
            "timing_rationale": self.timing_rationale,
        }


# ─────────────────────────────────────────────────────────────────────────────
# Ring Strategy Analyzer
# ─────────────────────────────────────────────────────────────────────────────

class RingStrategyAnalyzer:
    """
    环形成策略分析器
    
    决策逻辑:
    1. 芳香环优先 (更稳定，先构建)
    2. 简单环 → 复杂环
    3. 稠合环从共享边构建
    4. 大环最后 (RCM等)
    """
    
    def analyze(self, scaffold_report: Dict[str, Any]) -> RingPlan:
        """
        分析环形成策略
        
        Args:
            scaffold_report: 来自 ScaffoldAnalyzer 的报告
            
        Returns:
            RingPlan 环策略
        """
        ring_systems = scaffold_report.get("report", {}).get("ring_systems", [])
        
        if not ring_systems:
            return self._empty_plan()
        
        # 分析每个环
        constructions = []
        reactions_used = set()
        early_rings = []
        late_rings = []
        
        # 排序: 芳香 > 非芳香, 小 > 大
        sorted_rings = sorted(
            ring_systems,
            key=lambda r: (not r.get("is_aromatic", False), r.get("size", 6))
        )
        
        for i, ring in enumerate(sorted_rings):
            ring_id = ring.get("ring_id", i)
            size = ring.get("size", 6)
            is_aromatic = ring.get("is_aromatic", False)
            fused_with = ring.get("fused_with", [])
            
            # 选择构建方法
            method, difficulty = self._select_construction_method(ring)
            reactions_used.add(method)
            
            # 决定时机
            if is_aromatic or size == 6:
                step_order = i + 1
                early_rings.append(ring_id)
            else:
                step_order = len(sorted_rings) - i
                late_rings.append(ring_id)
            
            constructions.append(RingConstruction(
                ring_id=ring_id,
                ring_size=size,
                construction_method=method,
                step_order=step_order,
                is_aromatic=is_aromatic,
                difficulty=difficulty,
                dependencies=fused_with,
            ))
        
        # 生成开环断裂提案
        openings = self._generate_ring_openings(ring_systems)
        
        # 时机理由
        rationale = self._generate_rationale(early_rings, late_rings, ring_systems)
        
        return RingPlan(
            total_rings=len(ring_systems),
            construction_order=constructions,
            ring_opening_disconnections=openings,
            ring_forming_reactions_used=list(reactions_used),
            early_stage_rings=early_rings,
            late_stage_rings=late_rings,
            timing_rationale=rationale,
        )
    
    def _select_construction_method(self, ring: Dict) -> Tuple[str, str]:
        """选择合适的成环反应"""
        size = ring.get("size", 6)
        is_aromatic = ring.get("is_aromatic", False)
        is_hetero = ring.get("is_hetero", False)
        heteroatoms = ring.get("heteroatoms", [])
        
        # 芳香杂环
        if is_aromatic and is_hetero:
            if 'N' in heteroatoms:
                if size == 5:
                    return "Fischer Indole / Pyrrole synthesis", "moderate"
                if size == 6:
                    return "Hantzsch / Chichibabin", "moderate"
            return "Heterocycle synthesis", "hard"
        
        # 芳香碳环
        if is_aromatic and not is_hetero:
            return "Diels-Alder / aromatization", "moderate"
        
        # 非芳香环
        if size == 5:
            return "Nazarov / aldol cyclization", "moderate"
        if size == 6:
            return "Diels-Alder / Robinson annulation", "easy"
        if size >= 7:
            return "Ring-closing metathesis (RCM)", "hard"
        
        return "Cyclization", "moderate"
    
    def _generate_ring_openings(
        self, ring_systems: List[Dict]
    ) -> List[RingOpeningDisconnection]:
        """生成开环断裂提案"""
        openings = []
        
        for ring in ring_systems:
            ring_id = ring.get("ring_id", 0)
            size = ring.get("size", 6)
            is_aromatic = ring.get("is_aromatic", False)
            
            # 芳香环不宜直接开环
            if is_aromatic:
                continue
            
            # 根据环大小选择开环策略
            if size == 6:
                openings.append(RingOpeningDisconnection(
                    ring_id=ring_id,
                    bond_atoms=(0, 0),  # placeholder
                    opening_reaction="retro-Diels-Alder",
                    priority="HIGH",
                    fragments_description="diene + dienophile",
                    retro_reaction="Diels-Alder",
                ))
            elif size == 5:
                openings.append(RingOpeningDisconnection(
                    ring_id=ring_id,
                    bond_atoms=(0, 0),
                    opening_reaction="retro-aldol",
                    priority="MEDIUM",
                    fragments_description="aldehyde/ketone + enolate",
                    retro_reaction="Aldol cyclization",
                ))
            elif size >= 7:
                openings.append(RingOpeningDisconnection(
                    ring_id=ring_id,
                    bond_atoms=(0, 0),
                    opening_reaction="retro-RCM",
                    priority="HIGH",
                    fragments_description="diene precursor",
                    retro_reaction="Ring-closing metathesis",
                ))
        
        return openings
    
    def _generate_rationale(
        self,
        early: List[int],
        late: List[int],
        ring_systems: List[Dict],
    ) -> str:
        """生成时机理由"""
        parts = []
        
        if early:
            aromatic_early = sum(
                1 for r in ring_systems 
                if r.get("ring_id") in early and r.get("is_aromatic")
            )
            if aromatic_early:
                parts.append(f"芳香环优先构建 (稳定性)")
        
        if late:
            parts.append(f"复杂环/大环后期构建 (避免副反应)")
        
        fused_count = sum(1 for r in ring_systems if r.get("fused_with"))
        if fused_count:
            parts.append(f"稠环需考虑构建顺序")
        
        return "; ".join(parts) if parts else "Standard ring construction order"
    
    def _empty_plan(self) -> RingPlan:
        """无环时的空计划"""
        return RingPlan(
            total_rings=0,
            construction_order=[],
            ring_opening_disconnections=[],
            ring_forming_reactions_used=[],
            early_stage_rings=[],
            late_stage_rings=[],
            timing_rationale="No rings present - acyclic molecule",
        )


# ─────────────────────────────────────────────────────────────────────────────
# Convenience Function
# ─────────────────────────────────────────────────────────────────────────────

def plan_ring_strategy(scaffold_report: Dict[str, Any]) -> Dict[str, Any]:
    """
    规划环形成策略
    
    Args:
        scaffold_report: 来自 analyze_scaffold 的报告
        
    Returns:
        环策略字典
    """
    analyzer = RingStrategyAnalyzer()
    plan = analyzer.analyze(scaffold_report)
    return {
        "success": True,
        "plan": plan.to_dict(),
    }


# ─────────────────────────────────────────────────────────────────────────────
# CLI
# ─────────────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description="Ring Strategy Analyzer")
    parser.add_argument("scaffold_json", help="Scaffold report JSON or file path")
    args = parser.parse_args()
    
    # Try to parse as JSON or read file
    try:
        scaffold = json.loads(args.scaffold_json)
    except json.JSONDecodeError:
        from pathlib import Path
        scaffold = json.loads(Path(args.scaffold_json).read_text(encoding="utf-8"))
    
    result = plan_ring_strategy(scaffold)
    print(json.dumps(result, indent=2, ensure_ascii=False))
