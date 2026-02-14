"""
MultiRetro - Scaffold Analyzer
==============================
通用骨架分析器：识别核心骨架、环系统、战略断裂点。
适用于所有分子类型，为逆合成策略提供骨架层面的分析。
"""

from __future__ import annotations
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Tuple, Set
import json

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors
    from rdkit.Chem import BRICS
    from rdkit.Chem.Scaffolds import MurckoScaffold
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False


# ─────────────────────────────────────────────────────────────────────────────
# Data Classes
# ─────────────────────────────────────────────────────────────────────────────

@dataclass
class RingSystem:
    """环系统信息"""
    ring_id: int
    size: int
    atoms: List[int]
    is_aromatic: bool
    is_hetero: bool
    heteroatoms: List[str]
    fused_with: List[int] = field(default_factory=list)
    
    def to_dict(self) -> Dict[str, Any]:
        return {
            "ring_id": self.ring_id,
            "size": self.size,
            "atoms": self.atoms,
            "is_aromatic": self.is_aromatic,
            "is_hetero": self.is_hetero,
            "heteroatoms": self.heteroatoms,
            "fused_with": self.fused_with,
        }


@dataclass
class StrategicBond:
    """战略断裂键"""
    atom1_idx: int
    atom2_idx: int
    bond_type: str
    in_ring: bool
    ring_id: Optional[int]
    priority: str  # HIGH, MEDIUM, LOW
    reason: str
    fragments_complexity: Tuple[float, float] = (0.0, 0.0)
    
    def to_dict(self) -> Dict[str, Any]:
        return {
            "atom1_idx": self.atom1_idx,
            "atom2_idx": self.atom2_idx,
            "bond_type": self.bond_type,
            "in_ring": self.in_ring,
            "ring_id": self.ring_id,
            "priority": self.priority,
            "reason": self.reason,
            "fragments_complexity": list(self.fragments_complexity),
        }


@dataclass
class ScaffoldReport:
    """骨架分析报告"""
    smiles: str
    canonical_smiles: str
    scaffold_type: str  # fused_polycyclic, linear, branched, spiro, macrocyclic, simple
    murcko_scaffold: str
    ring_systems: List[RingSystem]
    ring_count: int
    fused_ring_count: int
    spiro_centers: List[int]
    bridgehead_atoms: List[int]
    strategic_bonds: List[StrategicBond]
    convergent_feasibility: float
    linear_feasibility: float
    complexity_score: float
    
    def to_dict(self) -> Dict[str, Any]:
        return {
            "smiles": self.smiles,
            "canonical_smiles": self.canonical_smiles,
            "scaffold_type": self.scaffold_type,
            "murcko_scaffold": self.murcko_scaffold,
            "ring_systems": [r.to_dict() for r in self.ring_systems],
            "ring_count": self.ring_count,
            "fused_ring_count": self.fused_ring_count,
            "spiro_centers": self.spiro_centers,
            "bridgehead_atoms": self.bridgehead_atoms,
            "strategic_bonds": [b.to_dict() for b in self.strategic_bonds],
            "convergent_feasibility": self.convergent_feasibility,
            "linear_feasibility": self.linear_feasibility,
            "complexity_score": self.complexity_score,
        }


# ─────────────────────────────────────────────────────────────────────────────
# Scaffold Analyzer
# ─────────────────────────────────────────────────────────────────────────────

class ScaffoldAnalyzer:
    """
    通用骨架分析器
    
    功能:
    1. 识别核心骨架类型 (稠环、线性、分支等)
    2. 分析环系统和融合点
    3. 识别战略断裂键
    4. 评估汇聚/线性合成可行性
    """
    
    def __init__(self):
        self._mol = None
        self._ring_info = None
    
    def analyze(self, smiles: str) -> ScaffoldReport:
        """
        分析分子骨架
        
        Args:
            smiles: 输入 SMILES
            
        Returns:
            ScaffoldReport 包含完整骨架分析
        """
        if not RDKIT_AVAILABLE:
            return self._fallback_report(smiles)
        
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return self._fallback_report(smiles)
        
        canonical = Chem.MolToSmiles(mol, canonical=True)
        # Re-parse canonical SMILES so atom indices are deterministic
        # and consistent with analyze_molecule / break_bond.
        mol_c = Chem.MolFromSmiles(canonical)
        if mol_c is None:
            mol_c = mol  # fallback
        self._mol = mol_c
        
        # 提取 Murcko 骨架
        murcko = self._get_murcko_scaffold(mol_c)
        
        # 分析环系统
        ring_systems = self._analyze_ring_systems(mol_c)
        fused_count = self._count_fused_rings(ring_systems)
        
        # 识别特殊中心
        spiro_centers = self._find_spiro_centers(mol_c)
        bridgehead = self._find_bridgehead_atoms(mol_c)
        
        # 确定骨架类型
        scaffold_type = self._determine_scaffold_type(
            ring_systems, fused_count, spiro_centers, bridgehead
        )
        
        # 识别战略断裂键
        strategic_bonds = self._find_strategic_bonds(mol_c, ring_systems)
        
        # 计算汇聚/线性可行性
        conv_score, lin_score = self._evaluate_strategy_feasibility(
            mol_c, ring_systems, strategic_bonds
        )
        
        # 复杂度评分
        complexity = self._calculate_complexity(mol_c, ring_systems)
        
        return ScaffoldReport(
            smiles=smiles,
            canonical_smiles=canonical,
            scaffold_type=scaffold_type,
            murcko_scaffold=murcko,
            ring_systems=ring_systems,
            ring_count=len(ring_systems),
            fused_ring_count=fused_count,
            spiro_centers=spiro_centers,
            bridgehead_atoms=bridgehead,
            strategic_bonds=strategic_bonds,
            convergent_feasibility=conv_score,
            linear_feasibility=lin_score,
            complexity_score=complexity,
        )
    
    def _get_murcko_scaffold(self, mol) -> str:
        """提取 Murcko 核心骨架"""
        try:
            core = MurckoScaffold.GetScaffoldForMol(mol)
            return Chem.MolToSmiles(core, canonical=True)
        except:
            return ""
    
    def _analyze_ring_systems(self, mol) -> List[RingSystem]:
        """分析所有环系统"""
        ring_info = mol.GetRingInfo()
        atom_rings = ring_info.AtomRings()
        
        ring_systems = []
        for i, ring_atoms in enumerate(atom_rings):
            atoms_list = list(ring_atoms)
            
            # 检查芳香性
            is_aromatic = all(mol.GetAtomWithIdx(a).GetIsAromatic() for a in atoms_list)
            
            # 检查杂原子
            heteroatoms = []
            is_hetero = False
            for atom_idx in atoms_list:
                atom = mol.GetAtomWithIdx(atom_idx)
                symbol = atom.GetSymbol()
                if symbol not in ['C', 'H']:
                    is_hetero = True
                    heteroatoms.append(symbol)
            
            ring_systems.append(RingSystem(
                ring_id=i,
                size=len(atoms_list),
                atoms=atoms_list,
                is_aromatic=is_aromatic,
                is_hetero=is_hetero,
                heteroatoms=heteroatoms,
            ))
        
        # 找融合关系
        for i, ring1 in enumerate(ring_systems):
            for j, ring2 in enumerate(ring_systems):
                if i >= j:
                    continue
                shared = set(ring1.atoms) & set(ring2.atoms)
                if len(shared) >= 2:  # 共享2个或更多原子 = 融合
                    ring1.fused_with.append(j)
                    ring2.fused_with.append(i)
        
        return ring_systems
    
    def _count_fused_rings(self, ring_systems: List[RingSystem]) -> int:
        """计算融合环数量"""
        fused_count = 0
        for ring in ring_systems:
            if ring.fused_with:
                fused_count += 1
        return fused_count
    
    def _find_spiro_centers(self, mol) -> List[int]:
        """找 spiro 中心 (共享单个原子的环)"""
        ring_info = mol.GetRingInfo()
        atom_rings = ring_info.AtomRings()
        
        spiro_centers = []
        atom_ring_count = {}
        
        for ring_atoms in atom_rings:
            for atom_idx in ring_atoms:
                atom_ring_count[atom_idx] = atom_ring_count.get(atom_idx, 0) + 1
        
        # 如果一个原子在多个环中但不是融合点，可能是 spiro
        for atom_idx, count in atom_ring_count.items():
            if count >= 2:
                # 检查是否真正是 spiro (不是融合)
                atom = mol.GetAtomWithIdx(atom_idx)
                if atom.GetDegree() == 4:  # sp3 碳
                    spiro_centers.append(atom_idx)
        
        return spiro_centers
    
    def _find_bridgehead_atoms(self, mol) -> List[int]:
        """找桥头原子"""
        bridgehead = []
        ring_info = mol.GetRingInfo()
        atom_rings = ring_info.AtomRings()
        
        atom_ring_count = {}
        for ring_atoms in atom_rings:
            for atom_idx in ring_atoms:
                atom_ring_count[atom_idx] = atom_ring_count.get(atom_idx, 0) + 1
        
        for atom_idx, count in atom_ring_count.items():
            if count >= 2:
                atom = mol.GetAtomWithIdx(atom_idx)
                # 桥头原子通常在多个环中且连接度高
                if atom.GetDegree() >= 3:
                    bridgehead.append(atom_idx)
        
        return bridgehead
    
    def _determine_scaffold_type(
        self,
        ring_systems: List[RingSystem],
        fused_count: int,
        spiro_centers: List[int],
        bridgehead: List[int],
    ) -> str:
        """确定骨架类型"""
        n_rings = len(ring_systems)
        
        if n_rings == 0:
            return "acyclic"
        if n_rings == 1:
            if ring_systems[0].size >= 12:
                return "macrocyclic"
            return "simple_monocyclic"
        
        if spiro_centers:
            return "spiro"
        if bridgehead and fused_count > 0:
            return "bridged_polycyclic"
        if fused_count >= 2:
            return "fused_polycyclic"
        if fused_count == 1:
            return "bicyclic_fused"
        
        return "polycyclic_separated"
    
    def _find_strategic_bonds(
        self,
        mol,
        ring_systems: List[RingSystem],
    ) -> List[StrategicBond]:
        """识别战略断裂键"""
        strategic = []
        
        # 1. BRICS 断裂点 (已知合成反应)
        try:
            brics_bonds = list(BRICS.FindBRICSBonds(mol))
            for (a1, a2), (label1, label2) in brics_bonds:
                bond = mol.GetBondBetweenAtoms(a1, a2)
                if bond:
                    strategic.append(StrategicBond(
                        atom1_idx=a1,
                        atom2_idx=a2,
                        bond_type=str(bond.GetBondType()),
                        in_ring=bond.IsInRing(),
                        ring_id=None,
                        priority="HIGH",
                        reason=f"BRICS:{label1}-{label2}",
                    ))
        except:
            pass
        
        # 2. 环融合键 (开环断裂)
        ring_atoms_all = set()
        for ring in ring_systems:
            ring_atoms_all.update(ring.atoms)
        
        for ring in ring_systems:
            if ring.fused_with:
                for other_id in ring.fused_with:
                    other_ring = ring_systems[other_id]
                    shared = set(ring.atoms) & set(other_ring.atoms)
                    if len(shared) == 2:
                        atoms = list(shared)
                        strategic.append(StrategicBond(
                            atom1_idx=atoms[0],
                            atom2_idx=atoms[1],
                            bond_type="ring_fusion",
                            in_ring=True,
                            ring_id=ring.ring_id,
                            priority="HIGH",
                            reason="ring_fusion_bond",
                        ))
        
        # 3. 连接环与侧链的键
        for bond in mol.GetBonds():
            a1, a2 = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
            a1_in_ring = a1 in ring_atoms_all
            a2_in_ring = a2 in ring_atoms_all
            
            if a1_in_ring != a2_in_ring:  # 一个在环内一个在环外
                strategic.append(StrategicBond(
                    atom1_idx=a1,
                    atom2_idx=a2,
                    bond_type="ring_linker",
                    in_ring=False,
                    ring_id=None,
                    priority="MEDIUM",
                    reason="ring_sidechain_linker",
                ))
        
        # 去重
        seen = set()
        unique = []
        for bond in strategic:
            key = tuple(sorted([bond.atom1_idx, bond.atom2_idx]))
            if key not in seen:
                seen.add(key)
                unique.append(bond)
        
        return unique
    
    def _evaluate_strategy_feasibility(
        self,
        mol,
        ring_systems: List[RingSystem],
        strategic_bonds: List[StrategicBond],
    ) -> Tuple[float, float]:
        """评估汇聚/线性可行性"""
        n_atoms = mol.GetNumHeavyAtoms()
        n_rings = len(ring_systems)
        n_strategic = len([b for b in strategic_bonds if b.priority == "HIGH"])
        
        # 汇聚可行性: 多个高优先级断裂点 → 适合汇聚
        conv_score = min(1.0, n_strategic / 3.0)
        
        # 线性可行性: 少环、链状 → 适合线性
        lin_score = max(0.0, 1.0 - n_rings / 5.0)
        
        # 大分子更适合汇聚
        if n_atoms > 30:
            conv_score = min(1.0, conv_score + 0.2)
            lin_score = max(0.0, lin_score - 0.2)
        
        return round(conv_score, 2), round(lin_score, 2)
    
    def _calculate_complexity(
        self,
        mol,
        ring_systems: List[RingSystem],
    ) -> float:
        """计算复杂度评分 (1-10)"""
        score = 1.0
        
        # 原子数
        n_atoms = mol.GetNumHeavyAtoms()
        score += n_atoms / 20.0
        
        # 环数
        score += len(ring_systems) * 0.5
        
        # 融合环
        fused = sum(1 for r in ring_systems if r.fused_with)
        score += fused * 0.3
        
        # 立体中心
        try:
            stereo = len(Chem.FindMolChiralCenters(mol))
            score += stereo * 0.4
        except:
            pass
        
        # 杂环
        hetero = sum(1 for r in ring_systems if r.is_hetero)
        score += hetero * 0.2
        
        return min(10.0, round(score, 2))
    
    def _fallback_report(self, smiles: str) -> ScaffoldReport:
        """RDKit 不可用时的回退"""
        return ScaffoldReport(
            smiles=smiles,
            canonical_smiles=smiles,
            scaffold_type="unknown",
            murcko_scaffold="",
            ring_systems=[],
            ring_count=0,
            fused_ring_count=0,
            spiro_centers=[],
            bridgehead_atoms=[],
            strategic_bonds=[],
            convergent_feasibility=0.5,
            linear_feasibility=0.5,
            complexity_score=5.0,
        )


# ─────────────────────────────────────────────────────────────────────────────
# Convenience Function
# ─────────────────────────────────────────────────────────────────────────────

def analyze_scaffold(smiles: str) -> Dict[str, Any]:
    """
    分析分子骨架
    
    Args:
        smiles: 输入 SMILES
        
    Returns:
        骨架分析报告字典
    """
    analyzer = ScaffoldAnalyzer()
    report = analyzer.analyze(smiles)
    return {
        "success": True,
        "report": report.to_dict(),
    }


# ─────────────────────────────────────────────────────────────────────────────
# CLI
# ─────────────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description="Scaffold Analyzer")
    parser.add_argument("smiles", help="SMILES to analyze")
    args = parser.parse_args()
    
    result = analyze_scaffold(args.smiles)
    print(json.dumps(result, indent=2, ensure_ascii=False))
