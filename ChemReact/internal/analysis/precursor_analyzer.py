"""
Precursor Analysis Module
=========================
Provides comprehensive chemical analysis for precursor molecules in synthesis routes.

核心功能：
- 对 rs（反应物）集合进行完整化学分析
- 复用现有的 rdkit_utils 模块
- 输出规范化的 MoleculeInfo 结构
"""

import sys
from pathlib import Path
from typing import Dict, Any, List, Optional

# 确保可以导入项目模块
_PROJECT_ROOT = Path(__file__).parent.parent.parent
if str(_PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(_PROJECT_ROOT))

from internal.normalize.smiles_normalizer import normalize_smiles


def _load_rdkit_modules():
    """延迟加载 rdkit_utils 模块"""
    from internal.rdkit_utils import analyze_structure
    from internal.rdkit_utils import calculate_features
    from internal.rdkit_utils import functional_groups
    return analyze_structure, calculate_features, functional_groups


def analyze_molecule(smiles: str) -> Dict[str, Any]:
    """
    对单个分子进行完整的化学分析，输出规范化的 MoleculeInfo 结构
    
    Args:
        smiles: 分子 SMILES（可含映射）
        
    Returns:
        MoleculeInfo 结构：
        {
            "smiles": {original, canonical, display},
            "analysis": {structure, features, functional_groups},
            "alerts": [...],
            "success": bool
        }
    """
    result = {
        "smiles": normalize_smiles(smiles),
        "analysis": {
            "structure": None,
            "features": None,
            "functional_groups": None
        },
        "alerts": [],
        "success": False,
        "error": None
    }
    
    try:
        analyze_mod, features_mod, fg_mod = _load_rdkit_modules()
        
        # 使用 display 版本进行分析（无映射）
        display_smiles = result["smiles"]["display"]
        
        # 结构分析
        structure_result = analyze_mod.analyze_smiles(display_smiles)
        if structure_result.get("success"):
            result["analysis"]["structure"] = structure_result
        
        # 特征计算
        features_result = features_mod.calculate_features(display_smiles)
        if features_result.get("success"):
            result["analysis"]["features"] = features_result
        
        # 官能团分析
        fg_result = fg_mod.analyze_functional_groups(display_smiles)
        if fg_result.get("success"):
            result["analysis"]["functional_groups"] = fg_result.get("functional_groups", {})
            result["alerts"].extend(fg_result.get("alerts", []))
        
        result["success"] = True
        
    except Exception as e:
        result["error"] = str(e)
    
    return result


def analyze_precursor_set(
    precursor_smiles_list: List[str],
    include_compatibility_check: bool = True
) -> Dict[str, Any]:
    """
    对前体集合进行完整的化学分析
    
    Args:
        precursor_smiles_list: 前体 SMILES 列表
        include_compatibility_check: 是否检查前体间兼容性
        
    Returns:
        {
            "precursors": [MoleculeInfo, ...],
            "summary": {
                "total": int,
                "analyzed": int,
                "failed": int,
                "all_functional_groups": {fg_name: count},
                "total_alerts": int
            },
            "compatibility_issues": [...],  # 前体间可能的问题
            "success": bool
        }
    """
    result = {
        "precursors": [],
        "summary": {
            "total": len(precursor_smiles_list),
            "analyzed": 0,
            "failed": 0,
            "all_functional_groups": {},
            "total_alerts": 0
        },
        "compatibility_issues": [],
        "success": False
    }
    
    all_fg = {}
    all_alerts = []
    
    for smi in precursor_smiles_list:
        if not smi or not smi.strip():
            continue
        
        mol_info = analyze_molecule(smi.strip())
        result["precursors"].append(mol_info)
        
        if mol_info["success"]:
            result["summary"]["analyzed"] += 1
            
            # 汇总官能团
            fgs = mol_info["analysis"].get("functional_groups", {})
            if isinstance(fgs, dict):
                for fg_name, count in fgs.items():
                    if count > 0:
                        all_fg[fg_name] = all_fg.get(fg_name, 0) + count
            
            # 汇总告警
            alerts = mol_info.get("alerts", [])
            all_alerts.extend(alerts)
        else:
            result["summary"]["failed"] += 1
    
    result["summary"]["all_functional_groups"] = all_fg
    result["summary"]["total_alerts"] = len(all_alerts)
    
    # 兼容性检查
    if include_compatibility_check:
        result["compatibility_issues"] = _check_compatibility(result["precursors"], all_fg)
    
    result["success"] = result["summary"]["analyzed"] > 0
    
    return result


def _check_compatibility(precursors: List[Dict[str, Any]], all_fg: Dict[str, int]) -> List[Dict[str, Any]]:
    """
    检查前体间的兼容性问题
    
    Args:
        precursors: 前体分析结果列表
        all_fg: 汇总的官能团计数
        
    Returns:
        兼容性问题列表
    """
    issues = []
    
    # 规则1：酸碱共存
    if all_fg.get("carboxylic_acid", 0) > 0 and all_fg.get("amine", 0) > 0:
        issues.append({
            "code": "ACID_BASE_COEXIST",
            "severity": "warning",
            "message": "Carboxylic acid and amine coexist - may form salt or require protection"
        })
    
    # 规则2：活性酰基与亲核试剂
    acyl_groups = all_fg.get("ester", 0) + all_fg.get("amide", 0)
    nucleophiles = all_fg.get("amine", 0) + all_fg.get("alcohol", 0) + all_fg.get("thiol", 0)
    if acyl_groups > 0 and nucleophiles > 1:
        issues.append({
            "code": "COMPETING_NUCLEOPHILES",
            "severity": "info",
            "message": "Multiple nucleophiles detected - consider chemoselectivity"
        })
    
    # 规则3：多卤素
    if all_fg.get("halide", 0) > 2:
        issues.append({
            "code": "MULTIPLE_HALIDES",
            "severity": "info",
            "message": "Multiple halides detected - competing substitution reactions possible"
        })
    
    # 规则4：氧化敏感基团
    if all_fg.get("aldehyde", 0) > 0:
        issues.append({
            "code": "OXIDATION_SENSITIVE",
            "severity": "info",
            "message": "Aldehyde present - may oxidize to carboxylic acid"
        })
    
    return issues


# CLI 入口
if __name__ == "__main__":
    import argparse
    import json
    
    parser = argparse.ArgumentParser(description="Analyze precursor molecules")
    parser.add_argument("smiles", nargs="+", help="SMILES strings to analyze")
    parser.add_argument("--single", "-s", action="store_true",
                       help="Analyze as single molecule instead of set")
    args = parser.parse_args()
    
    if args.single and len(args.smiles) == 1:
        result = analyze_molecule(args.smiles[0])
    else:
        result = analyze_precursor_set(args.smiles)
    
    print(json.dumps(result, indent=2, ensure_ascii=False, default=str))
