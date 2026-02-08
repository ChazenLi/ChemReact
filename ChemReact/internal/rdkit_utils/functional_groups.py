import argparse
import json
from typing import Dict, Any, List

from rdkit import Chem


_FG_SMARTS = {
    "alcohol": "[OX2H][CX4]",
    "phenol": "c[OX2H]",
    "amine": "[NX3;H2,H1;!$(NC=O)]",
    "amide": "[NX3][CX3](=[OX1])[#6]",
    "carboxylic_acid": "[CX3](=O)[OX2H1]",
    "ester": "[CX3](=O)[OX2][#6]",
    "ether": "[OD2]([#6])[#6]",
    "aldehyde": "[CX3H1](=O)[#6]",
    "ketone": "[#6][CX3](=O)[#6]",
    "halide": "[F,Cl,Br,I]",
    "nitrile": "[CX2]#N",
    "nitro": "[NX3+](=O)[O-]",
    "thiol": "[SX2H]",
    "sulfone": "[SX4](=O)(=O)([#6])[#6]",
    "boronic_acid": "[BX3]([OX2H])[OX2H]",
}

# 保护基识别 SMARTS 模式
_PROTECTING_GROUP_SMARTS = {
    # 氨基保护基
    "Boc": "[NX3][CX3](=O)OC(C)(C)C",  # tert-butyloxycarbonyl
    "Cbz": "[NX3][CX3](=O)OCc1ccccc1",  # benzyloxycarbonyl
    "Fmoc": "[NX3][CX3](=O)OCC1c2ccccc2-c2ccccc12",  # fluorenylmethyloxycarbonyl
    "Acetyl_N": "[NX3][CX3](=O)C",  # N-acetyl
    "Tosyl_N": "[NX3]S(=O)(=O)c1ccc(C)cc1",  # N-tosyl
    
    # 羟基保护基
    "TBS": "[OX2][Si](C)(C)C(C)(C)C",  # tert-butyldimethylsilyl
    "TIPS": "[OX2][Si](C(C)C)(C(C)C)C(C)C",  # triisopropylsilyl
    "TMS": "[OX2][Si](C)(C)C",  # trimethylsilyl
    "Bn": "[OX2]Cc1ccccc1",  # benzyl ether
    "PMB": "[OX2]Cc1ccc(OC)cc1",  # para-methoxybenzyl
    "THP": "[OX2]C1CCCCO1",  # tetrahydropyranyl
    "Acetyl_O": "[OX2][CX3](=O)C",  # O-acetyl
    
    # 羧酸保护基
    "MethylEster": "[CX3](=O)OC",  # methyl ester
    "EthylEster": "[CX3](=O)OCC",  # ethyl ester
    "tBuEster": "[CX3](=O)OC(C)(C)C",  # tert-butyl ester
    "BnEster": "[CX3](=O)OCc1ccccc1",  # benzyl ester
}


def _detect_protecting_groups(mol) -> Dict[str, Any]:
    """检测分子中的保护基团"""
    detected: List[Dict[str, Any]] = []
    pg_counts: Dict[str, int] = {}
    
    for name, smarts in _PROTECTING_GROUP_SMARTS.items():
        patt = Chem.MolFromSmarts(smarts)
        if patt is None:
            continue
        matches = mol.GetSubstructMatches(patt, uniquify=True)
        count = len(matches)
        if count > 0:
            pg_counts[name] = count
            # 分类保护基类型
            if name in ["Boc", "Cbz", "Fmoc", "Acetyl_N", "Tosyl_N"]:
                pg_type = "amine_protection"
            elif name in ["TBS", "TIPS", "TMS", "Bn", "PMB", "THP", "Acetyl_O"]:
                pg_type = "hydroxyl_protection"
            else:
                pg_type = "carboxyl_protection"
                
            detected.append({
                "name": name,
                "count": count,
                "type": pg_type,
                "atom_indices": [list(m) for m in matches],
            })
    
    return {
        "detected": detected,
        "counts": pg_counts,
        "has_protected_amine": any(d["type"] == "amine_protection" for d in detected),
        "has_protected_hydroxyl": any(d["type"] == "hydroxyl_protection" for d in detected),
        "has_protected_carboxyl": any(d["type"] == "carboxyl_protection" for d in detected),
    }


def analyze_functional_groups(smiles: str) -> Dict[str, Any]:
    result: Dict[str, Any] = {
        "input_smiles": smiles,
        "success": False,
        "error": None,
        "functional_groups": {},
        "protecting_groups": {},
        "alerts": [],
    }
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            result["error"] = "Invalid SMILES string"
            return result

        counts: Dict[str, int] = {}
        for name, smarts in _FG_SMARTS.items():
            patt = Chem.MolFromSmarts(smarts)
            if patt is None:
                continue
            matches = mol.GetSubstructMatches(patt, uniquify=True)
            counts[name] = len(matches)
        result["functional_groups"] = counts
        
        # 保护基检测
        result["protecting_groups"] = _detect_protecting_groups(mol)

        # 警告生成
        if counts.get("nitro", 0) > 0:
            result["alerts"].append("Nitro functionality may require chemoselective conditions.")
        if counts.get("halide", 0) > 2:
            result["alerts"].append("Multiple halides detected; competing substitutions may occur.")
        if counts.get("carboxylic_acid", 0) > 0 and counts.get("amine", 0) > 0:
            result["alerts"].append("Acid and amine coexistence may need protection/salt control.")
        
        # 保护基相关警告
        pg = result["protecting_groups"]
        if pg.get("has_protected_amine") and counts.get("amine", 0) > 0:
            result["alerts"].append("Both protected and unprotected amines present - selective deprotection may be needed.")
        if len(pg.get("detected", [])) > 2:
            result["alerts"].append("Multiple protecting groups detected - consider orthogonal deprotection strategy.")

        result["success"] = True
        return result
    except Exception as exc:
        result["error"] = str(exc)
        return result


def main() -> None:
    parser = argparse.ArgumentParser(description="Analyze functional groups for a SMILES string.")
    parser.add_argument("smiles", help="Input SMILES")
    args = parser.parse_args()
    print(json.dumps(analyze_functional_groups(args.smiles), indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
