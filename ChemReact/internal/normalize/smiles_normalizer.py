"""
SMILES Normalization Module
===========================
Provides utilities to normalize SMILES strings into multiple formats:
- original: 原始输入
- canonical: RDKit 规范化后的 SMILES
- display: 去除原子映射的版本，用于可视化

核心目的：解决带映射 SMILES (如 [CH3:1]) 导致的可视化失败问题
"""

import re
from typing import Dict, Any, List, Optional
from rdkit import Chem


def remove_atom_numbers(smiles: str) -> str:
    """
    Remove atom mapping numbers from SMILES.
    
    复用自 atom_mappers/utils.py，保持独立避免循环依赖
    
    Args:
        smiles: SMILES with or without atom numbers
        
    Returns:
        SMILES without atom numbers
        
    Example:
        >>> remove_atom_numbers("[CH3:1][OH:2]")
        '[CH3][OH]'
    """
    return re.sub(r'\[([A-Za-z][A-Za-z0-9+@-]*):\d+\]', r'[\1]', smiles)


def normalize_smiles(smiles: str) -> Dict[str, str]:
    """
    Normalize a SMILES string into three formats.
    
    Args:
        smiles: Input SMILES string (may contain atom mapping)
        
    Returns:
        {
            "original": 原始输入,
            "canonical": RDKit 规范化 (保留映射),
            "display": 去除映射，用于可视化
        }
        
    Example:
        >>> result = normalize_smiles("[CH3:1]O")
        >>> result["display"]  # 可直接用于 RDKit 可视化
        'CO'
    """
    result = {
        "original": smiles,
        "canonical": smiles,
        "display": smiles
    }
    
    if not smiles:
        return result
    
    try:
        # 尝试解析原始 SMILES（可能含映射）
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            # 无法解析，返回原始值
            return result
        
        # canonical: 保留原子映射的规范化 SMILES
        result["canonical"] = Chem.MolToSmiles(mol)
        
        # display: 去除映射后重新规范化，确保可视化兼容
        display_smiles = remove_atom_numbers(smiles)
        display_mol = Chem.MolFromSmiles(display_smiles)
        if display_mol:
            result["display"] = Chem.MolToSmiles(display_mol)
        else:
            result["display"] = display_smiles
            
    except Exception:
        # 保留原始值
        pass
    
    return result


def normalize_reaction_smiles(reaction_smiles: str) -> Dict[str, str]:
    """
    Normalize a reaction SMILES (Reactants>>Products) into multiple formats.
    
    Args:
        reaction_smiles: Reaction SMILES with >> separator
        
    Returns:
        {
            "original": 原始输入,
            "canonical": 规范化后的反应 SMILES (保留映射),
            "display": 去除映射的版本，用于可视化
        }
        
    Example:
        >>> result = normalize_reaction_smiles("[CH3:1]Br.[OH:2]>>[CH3:1][OH:2]")
        >>> result["display"]  # 可直接用于 RDKit 反应图渲染
        'CBr.O>>CO'
    """
    result = {
        "original": reaction_smiles,
        "canonical": reaction_smiles,
        "display": reaction_smiles
    }
    
    if not reaction_smiles or ">>" not in reaction_smiles:
        return result
    
    try:
        parts = reaction_smiles.split(">>")
        if len(parts) != 2:
            return result
        
        reactants_raw, products_raw = parts
        
        # 分别规范化反应物和产物
        reactants_normalized = []
        products_normalized = []
        reactants_display = []
        products_display = []
        
        for reactant in reactants_raw.split("."):
            reactant = reactant.strip()
            if reactant:
                norm = normalize_smiles(reactant)
                reactants_normalized.append(norm["canonical"])
                reactants_display.append(norm["display"])
        
        for product in products_raw.split("."):
            product = product.strip()
            if product:
                norm = normalize_smiles(product)
                products_normalized.append(norm["canonical"])
                products_display.append(norm["display"])
        
        result["canonical"] = ".".join(reactants_normalized) + ">>" + ".".join(products_normalized)
        result["display"] = ".".join(reactants_display) + ">>" + ".".join(products_display)
        
    except Exception:
        pass
    
    return result


def normalize_smiles_list(smiles_list: List[str]) -> List[Dict[str, str]]:
    """
    Batch normalize a list of SMILES strings.
    
    Args:
        smiles_list: List of SMILES strings
        
    Returns:
        List of normalized SMILES dictionaries
    """
    return [normalize_smiles(smi) for smi in smiles_list]


# CLI 入口
if __name__ == "__main__":
    import argparse
    import json
    
    parser = argparse.ArgumentParser(description="Normalize SMILES strings")
    parser.add_argument("smiles", help="Input SMILES string")
    parser.add_argument("--reaction", "-r", action="store_true", 
                       help="Treat as reaction SMILES")
    args = parser.parse_args()
    
    if args.reaction:
        result = normalize_reaction_smiles(args.smiles)
    else:
        result = normalize_smiles(args.smiles)
    
    print(json.dumps(result, indent=2, ensure_ascii=False))
