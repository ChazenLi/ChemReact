"""
MultiRetro - Chemoselectivity Analyzer
======================================
化学选择性分析器：识别竞争官能团、预测选择性、建议保护策略。

**重要局限性声明 (v3.10)**:
本模块的选择性分析基于:
- SMARTS 模式匹配识别官能团 (不考虑具体化学环境)
- 预定义的冲突对列表 (有限覆盖)
- 启发式的电子/位阻因素推测 (非量化计算)

LLM 应将本模块输出视为"选择性初筛"，对于关键决策需要:
1. 考虑具体分子的三维构象和位阻环境
2. 评估反应条件 (溶剂、温度、催化剂) 对选择性的影响
3. 判断是否存在本工具未覆盖的竞争反应路径
4. 对保护基建议进行正交性和后续步骤兼容性验证
"""

from __future__ import annotations
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Tuple
import json

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False


# ─────────────────────────────────────────────────────────────────────────────
# Functional Group Definitions with Reactivity
# ─────────────────────────────────────────────────────────────────────────────

FUNCTIONAL_GROUPS = {
    "alcohol": {
        "smarts": "[OX2H]",
        "reactivity": {
            "nucleophilicity": 7,
            "acidity": 4,
            "oxidation": 8,
            "acylation": 7,
        },
        "protection_options": ["TBS", "TIPS", "Bn", "Ac", "MOM"],
    },
    "phenol": {
        "smarts": "[c][OX2H]",
        "reactivity": {
            "nucleophilicity": 6,
            "acidity": 8,  # more acidic than alcohol
            "oxidation": 5,
            "acylation": 6,
        },
        "protection_options": ["Bn", "MOM", "Ac", "Me"],
    },
    "primary_amine": {
        "smarts": "[NX3H2][CX4]",
        "reactivity": {
            "nucleophilicity": 9,
            "acidity": 1,
            "acylation": 9,
            "reductive_amination": 8,
        },
        "protection_options": ["Boc", "Fmoc", "Cbz", "Ac"],
    },
    "secondary_amine": {
        "smarts": "[NX3H1]([CX4])[CX4]",
        "reactivity": {
            "nucleophilicity": 8,
            "acidity": 1,
            "acylation": 8,
        },
        "protection_options": ["Boc", "Cbz"],
    },
    "aldehyde": {
        "smarts": "[CX3H1](=O)",
        "reactivity": {
            "electrophilicity": 9,
            "reduction": 9,
            "aldol": 9,
            "wittig": 9,
        },
        "protection_options": ["acetal", "1,3-dioxolane"],
    },
    "ketone": {
        "smarts": "[CX3](=O)([#6])[#6]",
        "reactivity": {
            "electrophilicity": 7,
            "reduction": 7,
            "aldol": 7,
        },
        "protection_options": ["ketal", "1,3-dioxolane"],
    },
    "carboxylic_acid": {
        "smarts": "[CX3](=O)[OX2H]",
        "reactivity": {
            "acidity": 10,
            "acylation": 7,
            "reduction": 5,
        },
        "protection_options": ["tBu ester", "Me ester", "Bn ester"],
    },
    "ester": {
        "smarts": "[CX3](=O)[OX2][#6]",
        "reactivity": {
            "hydrolysis": 6,
            "reduction": 6,
            "transesterification": 5,
        },
        "protection_options": [],
    },
    "amide": {
        "smarts": "[CX3](=O)[NX3]",
        "reactivity": {
            "hydrolysis": 4,
            "reduction": 5,
        },
        "protection_options": [],
    },
    "nitrile": {
        "smarts": "[CX2]#N",
        "reactivity": {
            "hydrolysis": 4,
            "reduction": 6,
        },
        "protection_options": [],
    },
    "alkene": {
        "smarts": "[CX3]=[CX3]",
        "reactivity": {
            "electrophilic_addition": 7,
            "hydrogenation": 8,
            "oxidation": 7,
            "metathesis": 8,
        },
        "protection_options": [],
    },
    "alkyne": {
        "smarts": "[CX2]#[CX2]",
        "reactivity": {
            "electrophilic_addition": 6,
            "hydrogenation": 9,
            "click": 9,
        },
        "protection_options": ["TMS"],
    },
    "halide_aryl": {
        "smarts": "[c][F,Cl,Br,I]",
        "reactivity": {
            "cross_coupling": 9,
            "nucleophilic_substitution": 3,
        },
        "protection_options": [],
    },
    "halide_alkyl": {
        "smarts": "[CX4][F,Cl,Br,I]",
        "reactivity": {
            "nucleophilic_substitution": 8,
            "elimination": 6,
        },
        "protection_options": [],
    },
}


# ─────────────────────────────────────────────────────────────────────────────
# Data Classes
# ─────────────────────────────────────────────────────────────────────────────

@dataclass
class FunctionalGroupMatch:
    """官能团匹配结果"""
    fg_type: str
    atom_indices: List[int]
    count: int
    
    def to_dict(self) -> Dict[str, Any]:
        return {
            "fg_type": self.fg_type,
            "atom_indices": self.atom_indices,
            "count": self.count,
        }


@dataclass
class CompetingGroup:
    """竞争官能团"""
    fg_type: str
    count: int
    reaction_context: str
    differentiation_strategy: str
    electronic_factors: List[str]
    steric_factors: List[str]
    recommendation: str
    
    def to_dict(self) -> Dict[str, Any]:
        return {
            "fg_type": self.fg_type,
            "count": self.count,
            "reaction_context": self.reaction_context,
            "differentiation_strategy": self.differentiation_strategy,
            "electronic_factors": self.electronic_factors,
            "steric_factors": self.steric_factors,
            "recommendation": self.recommendation,
        }


@dataclass
class ProtectionRequirement:
    """保护需求"""
    fg_type: str
    atom_indices: List[int]
    reason: str
    recommended_pg: str
    install_timing: str  # early, middle, late
    remove_timing: str
    orthogonal_to: List[str]
    
    def to_dict(self) -> Dict[str, Any]:
        return {
            "fg_type": self.fg_type,
            "atom_indices": self.atom_indices,
            "reason": self.reason,
            "recommended_pg": self.recommended_pg,
            "install_timing": self.install_timing,
            "remove_timing": self.remove_timing,
            "orthogonal_to": self.orthogonal_to,
        }


@dataclass
class ReactivityConflict:
    """反应性冲突"""
    fg1_type: str
    fg2_type: str
    conflict_reaction: str
    severity: str  # low, medium, high
    resolution: str
    
    def to_dict(self) -> Dict[str, Any]:
        return {
            "fg1_type": self.fg1_type,
            "fg2_type": self.fg2_type,
            "conflict_reaction": self.conflict_reaction,
            "severity": self.severity,
            "resolution": self.resolution,
        }


@dataclass
class SelectivityReport:
    """选择性分析报告"""
    smiles: str
    functional_groups: List[FunctionalGroupMatch]
    competing_groups: List[CompetingGroup]
    protection_requirements: List[ProtectionRequirement]
    reactivity_conflicts: List[ReactivityConflict]
    reaction_order_constraints: List[str]
    overall_complexity: str  # low, medium, high
    
    def to_dict(self) -> Dict[str, Any]:
        return {
            "smiles": self.smiles,
            "functional_groups": [fg.to_dict() for fg in self.functional_groups],
            "competing_groups": [cg.to_dict() for cg in self.competing_groups],
            "protection_requirements": [pr.to_dict() for pr in self.protection_requirements],
            "reactivity_conflicts": [rc.to_dict() for rc in self.reactivity_conflicts],
            "reaction_order_constraints": self.reaction_order_constraints,
            "overall_complexity": self.overall_complexity,
        }


# ─────────────────────────────────────────────────────────────────────────────
# Chemoselectivity Analyzer
# ─────────────────────────────────────────────────────────────────────────────

class ChemoselectivityAnalyzer:
    """
    化学选择性分析器
    
    功能:
    1. 识别所有官能团
    2. 发现竞争官能团 (同类型多个)
    3. 分析反应性冲突
    4. 生成保护策略
    """
    
    def analyze(self, smiles: str) -> SelectivityReport:
        """
        分析分子的化学选择性
        
        Args:
            smiles: 输入 SMILES
            
        Returns:
            SelectivityReport
        """
        if not RDKIT_AVAILABLE:
            return self._fallback_report(smiles)
        
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return self._fallback_report(smiles)
        
        # 1. 识别所有官能团
        fg_matches = self._find_functional_groups(mol)
        
        # 2. 发现竞争官能团
        competing = self._find_competing_groups(fg_matches)
        
        # 3. 分析反应性冲突
        conflicts = self._find_reactivity_conflicts(fg_matches)
        
        # 4. 生成保护需求
        protections = self._determine_protection_needs(fg_matches, conflicts)
        
        # 5. 反应顺序约束
        order_constraints = self._derive_reaction_order(fg_matches, conflicts)
        
        # 6. 总体复杂度
        complexity = self._assess_complexity(competing, conflicts, protections)
        
        return SelectivityReport(
            smiles=smiles,
            functional_groups=fg_matches,
            competing_groups=competing,
            protection_requirements=protections,
            reactivity_conflicts=conflicts,
            reaction_order_constraints=order_constraints,
            overall_complexity=complexity,
        )
    
    def _find_functional_groups(self, mol) -> List[FunctionalGroupMatch]:
        """识别官能团"""
        matches = []
        
        for fg_name, fg_info in FUNCTIONAL_GROUPS.items():
            pattern = Chem.MolFromSmarts(fg_info["smarts"])
            if pattern is None:
                continue
            
            hits = mol.GetSubstructMatches(pattern)
            if hits:
                all_atoms = []
                for hit in hits:
                    all_atoms.extend(list(hit))
                
                matches.append(FunctionalGroupMatch(
                    fg_type=fg_name,
                    atom_indices=all_atoms,
                    count=len(hits),
                ))
        
        return matches
    
    def _find_competing_groups(
        self, fg_matches: List[FunctionalGroupMatch]
    ) -> List[CompetingGroup]:
        """发现竞争官能团 (同类型多个)"""
        competing = []
        
        for fg in fg_matches:
            if fg.count >= 2:
                # 同类型官能团 >= 2 个
                fg_info = FUNCTIONAL_GROUPS.get(fg.fg_type, {})
                reactivity = fg_info.get("reactivity", {})
                
                # 确定主要反应场景
                main_reactions = list(reactivity.keys())[:2]
                
                competing.append(CompetingGroup(
                    fg_type=fg.fg_type,
                    count=fg.count,
                    reaction_context=", ".join(main_reactions),
                    differentiation_strategy="steric/electronic",
                    electronic_factors=self._guess_electronic_factors(fg.fg_type),
                    steric_factors=["位阻差异", "环境取代基"],
                    recommendation=self._recommend_selectivity_approach(fg.fg_type, fg.count),
                ))
        
        return competing
    
    def _guess_electronic_factors(self, fg_type: str) -> List[str]:
        """推测电子因素"""
        if fg_type in ["phenol", "alcohol"]:
            return ["邻位取代基影响酸性", "共轭效应"]
        if fg_type in ["aldehyde", "ketone"]:
            return ["α-取代影响亲电性", "共轭稳定性"]
        if fg_type in ["primary_amine", "secondary_amine"]:
            return ["取代基供/吸电子效应"]
        return ["电子效应待分析"]
    
    def _recommend_selectivity_approach(self, fg_type: str, count: int) -> str:
        """推荐选择性策略"""
        if count == 2:
            return f"选择性保护其一，或利用位阻差异"
        if count >= 3:
            return f"需要正交保护策略，或分步官能团化"
        return "标准反应条件"
    
    def _find_reactivity_conflicts(
        self, fg_matches: List[FunctionalGroupMatch]
    ) -> List[ReactivityConflict]:
        """发现反应性冲突"""
        conflicts = []
        fg_types = [fg.fg_type for fg in fg_matches]
        
        # 常见冲突对 (v3.10 扩展: 4 → 12)
        conflict_pairs = [
            # 原有 4 种
            ("aldehyde", "primary_amine", "reductive_amination", "high",
             "醛与胺会直接反应形成亚胺，需控制条件或保护"),
            ("carboxylic_acid", "alcohol", "esterification", "medium",
             "酸与醇可能自酯化，需控制活化条件"),
            ("alkene", "aldehyde", "conjugate_addition", "medium",
             "共轭体系可能发生Michael加成"),
            ("halide_aryl", "primary_amine", "buchwald_side_reaction", "medium",
             "Buchwald条件下胺可能发生副反应"),
            # v3.10 新增 8 种
            ("aldehyde", "secondary_amine", "enamine_formation", "medium",
             "醛/酮与仲胺可能形成烯胺，需控制条件"),
            ("alkene", "alkyne", "selectivity_competition", "medium",
             "烯烃和炔烃共存时加氢/氧化选择性需要控制"),
            ("alcohol", "primary_amine", "nucleophile_competition", "low",
             "醇和胺同为亲核试剂，酰化反应中胺通常更快但需确认"),
            ("phenol", "primary_amine", "acylation_competition", "medium",
             "酚和胺竞争酰化，胺通常优先但酚的O-酰化不可忽略"),
            ("halide_alkyl", "alcohol", "sn2_competition", "medium",
             "卤代烷与醇的SN2反应可能与目标反应竞争"),
            ("nitrile", "primary_amine", "amidine_formation", "low",
             "腈与胺在强条件下可能形成脒"),
            ("alkyne", "primary_amine", "hydroamination", "low",
             "炔烃与胺在催化条件下可能发生氢胺化"),
            ("ester", "primary_amine", "aminolysis", "medium",
             "酯与胺可能发生氨解生成酰胺，尤其在加热条件下"),
        ]
        
        for fg1, fg2, reaction, severity, resolution in conflict_pairs:
            if fg1 in fg_types and fg2 in fg_types:
                conflicts.append(ReactivityConflict(
                    fg1_type=fg1,
                    fg2_type=fg2,
                    conflict_reaction=reaction,
                    severity=severity,
                    resolution=resolution,
                ))
        
        return conflicts
    
    def _determine_protection_needs(
        self,
        fg_matches: List[FunctionalGroupMatch],
        conflicts: List[ReactivityConflict],
    ) -> List[ProtectionRequirement]:
        """确定保护需求"""
        protections = []
        
        # 基于冲突确定保护
        conflict_fgs = set()
        for c in conflicts:
            if c.severity in ["high", "medium"]:
                conflict_fgs.add(c.fg1_type)
                conflict_fgs.add(c.fg2_type)
        
        for fg in fg_matches:
            fg_info = FUNCTIONAL_GROUPS.get(fg.fg_type, {})
            options = fg_info.get("protection_options", [])
            
            if not options:
                continue
            
            # 需要保护的情况
            needs_protection = (
                fg.fg_type in conflict_fgs or
                fg.count >= 2  # 多个同类型
            )
            
            if needs_protection:
                protections.append(ProtectionRequirement(
                    fg_type=fg.fg_type,
                    atom_indices=fg.atom_indices[:2],  # 代表性原子
                    reason="reactivity_conflict" if fg.fg_type in conflict_fgs else "selectivity",
                    recommended_pg=options[0] if options else "",
                    install_timing="early",
                    remove_timing="late",
                    orthogonal_to=self._get_orthogonal_pgs(options[0] if options else ""),
                ))
        
        return protections
    
    def _get_orthogonal_pgs(self, pg: str) -> List[str]:
        """获取正交保护基"""
        orthogonality = {
            "Boc": ["Fmoc", "Cbz", "TBS"],
            "Fmoc": ["Boc", "Cbz", "TBS"],
            "TBS": ["Boc", "Fmoc", "Ac"],
            "Bn": ["TBS", "TIPS", "Boc"],
        }
        return orthogonality.get(pg, [])
    
    def _derive_reaction_order(
        self,
        fg_matches: List[FunctionalGroupMatch],
        conflicts: List[ReactivityConflict],
    ) -> List[str]:
        """推导反应顺序约束"""
        constraints = []
        
        for c in conflicts:
            if c.severity == "high":
                constraints.append(f"保护 {c.fg1_type} 后再处理 {c.fg2_type}")
        
        # 通用规则
        fg_types = [fg.fg_type for fg in fg_matches]
        if "aldehyde" in fg_types:
            constraints.append("醛基应早期引入或保护")
        if "primary_amine" in fg_types:
            constraints.append("胺基需保护以避免副反应")
        
        return constraints
    
    def _assess_complexity(
        self,
        competing: List[CompetingGroup],
        conflicts: List[ReactivityConflict],
        protections: List[ProtectionRequirement],
    ) -> str:
        """评估总体复杂度"""
        score = 0
        score += len(competing) * 2
        score += sum(1 for c in conflicts if c.severity == "high") * 3
        score += sum(1 for c in conflicts if c.severity == "medium") * 1
        score += len(protections)
        
        if score <= 2:
            return "low"
        if score <= 6:
            return "medium"
        return "high"
    
    def _fallback_report(self, smiles: str) -> SelectivityReport:
        """RDKit 不可用时的回退"""
        return SelectivityReport(
            smiles=smiles,
            functional_groups=[],
            competing_groups=[],
            protection_requirements=[],
            reactivity_conflicts=[],
            reaction_order_constraints=[],
            overall_complexity="unknown",
        )


# ─────────────────────────────────────────────────────────────────────────────
# Convenience Function
# ─────────────────────────────────────────────────────────────────────────────

def analyze_chemoselectivity(smiles: str) -> Dict[str, Any]:
    """
    分析化学选择性
    
    Args:
        smiles: 输入 SMILES
        
    Returns:
        选择性分析报告字典
    """
    analyzer = ChemoselectivityAnalyzer()
    report = analyzer.analyze(smiles)
    return {
        "success": True,
        "report": report.to_dict(),
        "analysis_limitations": {
            "methodology": "smarts_pattern_matching",
            "description": (
                "选择性分析基于 SMARTS 模式匹配和预定义冲突规则，"
                "不考虑具体分子构象、溶剂效应和反应条件。"
            ),
            "coverage_gaps": [
                "冲突对列表有限 — 仅覆盖 12 种常见冲突模式",
                "电子/位阻因素为定性推测 — 不进行量化计算",
                "保护基建议未验证后续步骤兼容性",
                "未考虑溶剂和温度对选择性的影响",
            ],
            "llm_action_required": (
                "LLM 应基于具体反应条件和分子环境独立评估选择性，"
                "特别关注: (1) 位阻差异是否足以区分竞争基团; "
                "(2) 电子效应的定量影响; "
                "(3) 保护基方案的正交性和全局兼容性"
            ),
        },
    }


# ─────────────────────────────────────────────────────────────────────────────
# CLI
# ─────────────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description="Chemoselectivity Analyzer")
    parser.add_argument("smiles", help="SMILES to analyze")
    args = parser.parse_args()
    
    result = analyze_chemoselectivity(args.smiles)
    print(json.dumps(result, indent=2, ensure_ascii=False))
