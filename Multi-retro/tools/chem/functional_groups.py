"""
Functional Group Analysis Module
================================
Identifies functional groups, protecting groups, and generates compatibility alerts.

Migrated from internal/rdkit_utils/functional_groups.py.
v2: 27 FG SMARTS + protecting group detection + dangerous combination alerts.

Public API:
    analyze_functional_groups(smiles) → Dict[str, Any]
"""

import argparse
import json
from typing import Any, Dict, List

from rdkit import Chem

from tools.chem.mol_parser import parse_smiles


# ─── 原始 15 FG ───
_FG_SMARTS_CORE = {
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

# ─── 新增 12 FG (v2) ───
_FG_SMARTS_V2 = {
    "epoxide": "[CX4]1[OX2][CX4]1",
    "acyl_halide": "[CX3](=O)[F,Cl,Br,I]",
    "anhydride": "[CX3](=O)[OX2][CX3](=O)",
    "azide": "[N-]=[N+]=[N]",
    "isocyanate": "[NX2]=C=O",
    "imine": "[CX3]=[NX2]",
    "phosphonate": "[PX4](=O)([OX2])([OX2])",
    "sulfoxide": "[SX3](=O)([#6])[#6]",
    "vinyl": "[CX3]=[CX3]",
    "enol": "[OX2H][CX3]=[CX3]",
    "diazo": "[#6]=[N+]=[N-]",
    "hydrazine": "[NX3][NX3]",
}

# 合并完整 FG 字典
_FG_SMARTS: Dict[str, str] = {**_FG_SMARTS_CORE, **_FG_SMARTS_V2}

# 保护基识别 SMARTS 模式
_PROTECTING_GROUP_SMARTS = {
    # 氨基保护基
    "Boc": "[NX3][CX3](=O)OC(C)(C)C",
    "Cbz": "[NX3][CX3](=O)OCc1ccccc1",
    "Fmoc": "[NX3][CX3](=O)OCC1c2ccccc2-c2ccccc12",
    "Acetyl_N": "[NX3][CX3](=O)C",
    "Tosyl_N": "[NX3]S(=O)(=O)c1ccc(C)cc1",
    # 羟基保护基
    "TBS": "[OX2][Si](C)(C)C(C)(C)C",
    "TIPS": "[OX2][Si](C(C)C)(C(C)C)C(C)C",
    "TMS": "[OX2][Si](C)(C)C",
    "Bn": "[OX2]Cc1ccccc1",
    "PMB": "[OX2]Cc1ccc(OC)cc1",
    "THP": "[OX2]C1CCCCO1",
    "Acetyl_O": "[OX2][CX3](=O)C",
    # 羧酸保护基
    "MethylEster": "[CX3](=O)OC",
    "EthylEster": "[CX3](=O)OCC",
    "tBuEster": "[CX3](=O)OC(C)(C)C",
    "BnEster": "[CX3](=O)OCc1ccccc1",
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
            if name in ("Boc", "Cbz", "Fmoc", "Acetyl_N", "Tosyl_N"):
                pg_type = "amine_protection"
            elif name in ("TBS", "TIPS", "TMS", "Bn", "PMB", "THP", "Acetyl_O"):
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


def _generate_alerts(counts: Dict[str, int], pg: Dict[str, Any]) -> List[str]:
    """生成告警列表（原有 + v2 + 危险组合检测）"""
    alerts: List[str] = []

    # 原有告警
    if counts.get("nitro", 0) > 0:
        alerts.append("Nitro functionality may require chemoselective conditions.")
    if counts.get("halide", 0) > 2:
        alerts.append("Multiple halides detected; competing substitutions may occur.")
    if counts.get("carboxylic_acid", 0) > 0 and counts.get("amine", 0) > 0:
        alerts.append("Acid and amine coexistence may need protection/salt control.")
    if pg.get("has_protected_amine") and counts.get("amine", 0) > 0:
        alerts.append("Both protected and unprotected amines present - selective deprotection may be needed.")
    if len(pg.get("detected", [])) > 2:
        alerts.append("Multiple protecting groups detected - consider orthogonal deprotection strategy.")

    # ─── v2 告警 ───
    if counts.get("azide", 0) > 0:
        alerts.append("Azide detected — potential safety hazard, handle with care.")
    if counts.get("acyl_halide", 0) > 0:
        alerts.append("Acyl halide present — highly reactive, moisture-sensitive.")
    if counts.get("epoxide", 0) > 0 and counts.get("amine", 0) > 0:
        alerts.append("Epoxide + amine coexist — uncontrolled ring opening risk.")
    if counts.get("isocyanate", 0) > 0:
        alerts.append("Isocyanate detected — moisture-sensitive, toxic.")
    if counts.get("diazo", 0) > 0:
        alerts.append("Diazo compound detected — potentially explosive, handle with extreme care.")
    if counts.get("anhydride", 0) > 0 and (counts.get("amine", 0) > 0 or counts.get("alcohol", 0) > 0):
        alerts.append("Anhydride + nucleophile coexist — competing acylation risk.")
    if counts.get("epoxide", 0) > 0 and counts.get("thiol", 0) > 0:
        alerts.append("Epoxide + thiol coexist — rapid ring opening expected.")

    # ─── 危险官能团组合检测 (Req 10.4) ───
    if counts.get("nitro", 0) > 0 and counts.get("amine", 0) > 0:
        alerts.append("硝基+胺组合检测 — 潜在爆炸风险，需要严格控制反应条件")
    if counts.get("azide", 0) > 0 and (counts.get("ketone", 0) > 0 or counts.get("aldehyde", 0) > 0):
        alerts.append("叠氮+羰基组合检测 — Curtius重排风险，注意安全")

    return alerts


def analyze_functional_groups(smiles: str) -> Dict[str, Any]:
    """
    分析分子的官能团、保护基和告警。

    Args:
        smiles: 分子 SMILES

    Returns:
        {
            "input_smiles": str,
            "success": bool,
            "error": str | None,
            "functional_groups": {fg_name: count, ...},
            "protecting_groups": {...},
            "alerts": [str, ...]
        }
    """
    result: Dict[str, Any] = {
        "input_smiles": smiles,
        "success": False,
        "error": None,
        "functional_groups": {},
        "protecting_groups": {},
        "alerts": [],
    }
    try:
        mol = parse_smiles(smiles)
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

        pg = _detect_protecting_groups(mol)
        result["protecting_groups"] = pg

        result["alerts"] = _generate_alerts(counts, pg)
        result["success"] = True
        return result
    except Exception as exc:
        result["error"] = str(exc)
        return result


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Analyze functional groups for a SMILES string."
    )
    parser.add_argument("smiles", help="Input SMILES")
    args = parser.parse_args()
    print(json.dumps(analyze_functional_groups(args.smiles), indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
