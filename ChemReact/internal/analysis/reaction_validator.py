"""
Reaction Validation Module
==========================
Provides comprehensive reaction validation with structured output including:
- Atom balance checking
- Bond changes (formed/broken)
- Selectivity notes
- Reaction type classification

核心功能：验证反应的化学可行性并输出结构化信息
"""

import sys
from pathlib import Path
from typing import Dict, Any, List, Optional
from collections import Counter

# 确保可以导入项目模块
_PROJECT_ROOT = Path(__file__).parent.parent.parent
if str(_PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(_PROJECT_ROOT))

from internal.normalize.smiles_normalizer import normalize_reaction_smiles

try:
    from rdkit import Chem
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False


def _load_atom_mapper_utils():
    """延迟加载 atom_mappers/utils 模块"""
    from internal.atom_mappers import utils
    return utils


def _load_reaction_classifier():
    """延迟加载 reaction_classifier 模块"""
    from internal.reaction_classifier import classifier
    return classifier


def _element_counter_from_smiles(smiles: str) -> Counter:
    """从 SMILES 计算元素计数（含隐式氢）"""
    counter = Counter()
    if not RDKIT_AVAILABLE:
        return counter
    
    for frag in smiles.split("."):
        frag = frag.strip()
        if not frag:
            continue
        mol = Chem.MolFromSmiles(frag)
        if mol is None:
            continue
        mol = Chem.AddHs(mol)
        for atom in mol.GetAtoms():
            counter[atom.GetSymbol()] += 1
    
    return counter


def validate_reaction(
    reaction_smiles: str,
    reaction_class: str = "",
    include_bond_changes: bool = True,
    include_selectivity: bool = True
) -> Dict[str, Any]:
    """
    对反应进行完整验证，输出结构化 ReactionValidation 结果
    
    Args:
        reaction_smiles: 反应 SMILES (rs >> ps)
        reaction_class: 可选的反应类型提示
        include_bond_changes: 是否包含成键断键分析
        include_selectivity: 是否包含选择性简评
        
    Returns:
        ReactionValidation 结构：
        {
            "reaction_smiles": {original, canonical, display},
            "atom_balance": {pass, delta, inferred_byproducts},
            "bond_changes": {formed, broken},
            "reaction_type": {predicted, confidence},
            "selectivity_notes": [...],
            "verdict": "PASS"|"FAIL"|"CONDITIONAL",
            "issues": [...]
        }
    """
    result = {
        "reaction_smiles": normalize_reaction_smiles(reaction_smiles),
        "atom_balance": {
            "pass": True,
            "delta": {},
            "inferred_byproducts": []
        },
        "bond_changes": {
            "formed": [],
            "broken": []
        },
        "reaction_type": {
            "predicted": reaction_class or "unknown",
            "confidence": 0.0,
            "mechanism_hint": ""
        },
        "selectivity_notes": [],
        "verdict": "CONDITIONAL",
        "issues": [],
        "success": False
    }
    
    if ">>" not in reaction_smiles:
        result["issues"].append({
            "code": "REACTION_FORMAT_INVALID",
            "message": "Missing '>>' separator"
        })
        result["verdict"] = "FAIL"
        return result
    
    try:
        parts = reaction_smiles.split(">>")
        if len(parts) != 2:
            result["issues"].append({
                "code": "REACTION_FORMAT_INVALID",
                "message": "Invalid reaction format"
            })
            result["verdict"] = "FAIL"
            return result
        
        reactants, products = parts
        
        # 1. 原子守恒检查
        reactant_counter = _element_counter_from_smiles(reactants)
        product_counter = _element_counter_from_smiles(products)
        
        # 计算原子差分 (正数=反应物过剩，负数=产物过剩)
        # 注意: Counter 的减法会丢弃负值，需要手动计算完整差分
        all_elements = set(reactant_counter.keys()) | set(product_counter.keys())
        delta = {elem: reactant_counter.get(elem, 0) - product_counter.get(elem, 0) 
                 for elem in all_elements}
        
        # 过滤零值
        delta = {k: v for k, v in delta.items() if v != 0}
        result["atom_balance"]["delta"] = dict(delta)
        
        if delta:
            result["atom_balance"]["pass"] = False
            result["issues"].append({
                "code": "ATOM_CONSERVATION_FAIL",
                "message": f"Atom imbalance: {dict(delta)}"
            })
            
            # 推断可能的副产物
            result["atom_balance"]["inferred_byproducts"] = _infer_byproducts(delta)
        
        # 2. 成键断键分析
        if include_bond_changes:
            try:
                utils = _load_atom_mapper_utils()
                bond_changes = utils.get_bond_changes(reaction_smiles)
                result["bond_changes"] = bond_changes
            except Exception as e:
                result["issues"].append({
                    "code": "BOND_ANALYSIS_FAILED",
                    "message": str(e)
                })
        
        # 3. 反应类型分类
        try:
            classifier = _load_reaction_classifier()
            type_info = classifier.classify_reaction(
                reactants, products, reaction_class
            )
            result["reaction_type"]["predicted"] = type_info.get("reaction_type", reaction_class or "unknown")
            result["reaction_type"]["confidence"] = type_info.get("confidence", 0.0)
            
            # 提取机理提示
            candidates = type_info.get("byproduct_candidates", [])
            if candidates:
                result["reaction_type"]["mechanism_hint"] = candidates[0].get("mechanism", "")
        except Exception:
            pass
        
        # 4. 选择性简评
        if include_selectivity:
            result["selectivity_notes"] = _generate_selectivity_notes(
                result["reaction_type"]["predicted"],
                result["bond_changes"]
            )
        
        # 5. 判定结果
        if not result["atom_balance"]["pass"]:
            if result["atom_balance"]["inferred_byproducts"]:
                result["verdict"] = "CONDITIONAL"
            else:
                result["verdict"] = "FAIL"
        elif result["issues"]:
            result["verdict"] = "CONDITIONAL"
        else:
            result["verdict"] = "PASS"
        
        result["success"] = True
        
    except Exception as e:
        result["issues"].append({
            "code": "VALIDATION_ERROR",
            "message": str(e)
        })
        result["verdict"] = "FAIL"
    
    return result


def _infer_byproducts(delta: Dict[str, int]) -> List[Dict[str, Any]]:
    """
    根据原子差分推断可能的副产物
    
    正数 delta = 反应物过剩 → 推断为副产物
    负数 delta = 产物过剩 → 推断为缺失反应物
    
    Returns:
        包含 byproducts 和 missing_reactants 的推断结果
    """
    byproducts = []
    
    # 分离正负 delta
    excess_reactant = {k: v for k, v in delta.items() if v > 0}  # 需要副产物消耗
    excess_product = {k: -v for k, v in delta.items() if v < 0}  # 缺失反应物
    
    # 常见副产物模式 (按分子量排序，从小到大匹配)
    patterns = [
        {"formula": "H2", "elements": {"H": 2}, "name": "Hydrogen", "smiles": "[H][H]"},
        {"formula": "H2O", "elements": {"H": 2, "O": 1}, "name": "Water", "smiles": "O"},
        {"formula": "HCl", "elements": {"H": 1, "Cl": 1}, "name": "Hydrochloric acid", "smiles": "Cl"},
        {"formula": "HBr", "elements": {"H": 1, "Br": 1}, "name": "Hydrobromic acid", "smiles": "Br"},
        {"formula": "HI", "elements": {"H": 1, "I": 1}, "name": "Hydroiodic acid", "smiles": "I"},
        {"formula": "NH3", "elements": {"N": 1, "H": 3}, "name": "Ammonia", "smiles": "N"},
        {"formula": "CO2", "elements": {"C": 1, "O": 2}, "name": "Carbon dioxide", "smiles": "O=C=O"},
        {"formula": "N2", "elements": {"N": 2}, "name": "Nitrogen", "smiles": "N#N"},
        {"formula": "CH3OH", "elements": {"C": 1, "H": 4, "O": 1}, "name": "Methanol", "smiles": "CO"},
        {"formula": "C2H5OH", "elements": {"C": 2, "H": 6, "O": 1}, "name": "Ethanol", "smiles": "CCO"},
        {"formula": "AcOH", "elements": {"C": 2, "H": 4, "O": 2}, "name": "Acetic acid", "smiles": "CC(=O)O"},
    ]
    
    remaining = dict(excess_reactant)
    
    # 贪心匹配副产物（可匹配多次）
    for pattern in patterns:
        while True:
            # 检查是否可以用这个副产物解释部分差分
            can_use = all(
                remaining.get(elem, 0) >= count 
                for elem, count in pattern["elements"].items()
            )
            
            if not can_use:
                break
                
            # 找到已有副产物或新建
            existing = next((bp for bp in byproducts if bp["formula"] == pattern["formula"]), None)
            if existing:
                existing["count"] = existing.get("count", 1) + 1
            else:
                byproducts.append({
                    "formula": pattern["formula"],
                    "name": pattern["name"],
                    "smiles": pattern["smiles"],
                    "count": 1,
                    "type": "byproduct",
                    "confidence": 0.8
                })
            
            # 更新剩余差分
            for elem, count in pattern["elements"].items():
                remaining[elem] = remaining.get(elem, 0) - count
                if remaining[elem] == 0:
                    del remaining[elem]
    
    # 如果有产物过剩，推断缺失反应物
    if excess_product:
        for pattern in patterns:
            while True:
                can_use = all(
                    excess_product.get(elem, 0) >= count 
                    for elem, count in pattern["elements"].items()
                )
                
                if not can_use:
                    break
                    
                byproducts.append({
                    "formula": pattern["formula"],
                    "name": pattern["name"],
                    "smiles": pattern["smiles"],
                    "count": 1,
                    "type": "missing_reactant",
                    "confidence": 0.6  # Lower confidence for missing reactants
                })
                
                for elem, count in pattern["elements"].items():
                    excess_product[elem] = excess_product.get(elem, 0) - count
                    if excess_product[elem] == 0:
                        del excess_product[elem]
    
    # 调整置信度：如果有未解释的原子，降低置信度
    unexplained = sum(abs(v) for v in remaining.values()) + sum(excess_product.values())
    if unexplained > 0 and byproducts:
        penalty = min(0.4, unexplained * 0.1)  # 每个未解释原子降低 10%，最多降 40%
        for bp in byproducts:
            bp["confidence"] = round(bp["confidence"] - penalty, 2)
    
    return byproducts


def _generate_selectivity_notes(
    reaction_type: str,
    bond_changes: Dict[str, List]
) -> List[str]:
    """生成选择性简评"""
    notes = []
    
    # 基于反应类型的选择性提示
    selectivity_hints = {
        "esterification": [
            "Carboxylic acid is more reactive than alcohol",
            "Consider acid catalyst for faster reaction"
        ],
        "amide_coupling": [
            "Primary amines more reactive than secondary",
            "Watch for racemization at alpha-carbon"
        ],
        "nucleophilic_substitution": [
            "SN2: steric hindrance affects rate",
            "SN1: carbocation stability determines pathway"
        ],
        "oxidation": [
            "Primary alcohols → aldehydes → carboxylic acids",
            "Control conditions to stop at aldehyde stage"
        ],
        "reduction": [
            "Nitriles reduce to amines",
            "Aldehydes reduce to alcohols"
        ],
        "cross_coupling": [
            "Palladium-catalyzed: halide leaving group order Br > I > Cl",
            "Aryl halides more reactive than alkyl"
        ]
    }
    
    # 查找匹配的反应类型
    for key, hints in selectivity_hints.items():
        if key.lower() in reaction_type.lower():
            notes.extend(hints)
            break
    
    # 基于键变化的通用提示
    if bond_changes.get("formed"):
        if len(bond_changes["formed"]) > 1:
            notes.append("Multiple bonds forming - check regioselectivity")
    
    if bond_changes.get("broken"):
        if len(bond_changes["broken"]) > 1:
            notes.append("Multiple bonds breaking - consider reaction order")
    
    return notes


# CLI 入口
if __name__ == "__main__":
    import argparse
    import json
    
    parser = argparse.ArgumentParser(description="Validate a reaction")
    parser.add_argument("reaction", help="Reaction SMILES (rs>>ps)")
    parser.add_argument("--class", "-c", dest="reaction_class", default="",
                       help="Reaction class hint")
    args = parser.parse_args()
    
    result = validate_reaction(args.reaction, args.reaction_class)
    print(json.dumps(result, indent=2, ensure_ascii=False, default=str))
