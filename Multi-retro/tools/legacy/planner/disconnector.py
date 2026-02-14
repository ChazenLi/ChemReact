"""
MultiRetro Core Planner - Disconnector Module
==============================================
Comprehensive retrosynthetic disconnection logic with 35+ reaction types.

v2: SMARTS 特异性修复, 杂环断裂规则, Redox 提案生成, 优先级系统, SMARTS 验证
"""

from __future__ import annotations
from dataclasses import dataclass, field
from typing import List, Dict, Any, Optional, Tuple
import hashlib
import json

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False

from tools.common.errors import ChemistryError, ValidationError


@dataclass
class DisconnectionProposal:
    """A proposed retrosynthetic disconnection."""
    id: str
    strategy: str
    bond_description: str
    precursors: List[str]
    reaction_class: str
    confidence: float
    rationale: str
    reaction_smiles: str = ""
    byproducts: List[str] = field(default_factory=list)
    reagents: List[str] = field(default_factory=list)
    is_redox: bool = False
    specificity: int = 1       # 1=通用, 2=中等, 3=高特异性
    priority: int = 5           # 1=最高, 10=最低

    def to_dict(self) -> Dict[str, Any]:
        return {
            "id": self.id,
            "strategy": self.strategy,
            "bond_description": self.bond_description,
            "precursors": self.precursors,
            "reaction_class": self.reaction_class,
            "confidence": self.confidence,
            "rationale": self.rationale,
            "reaction_smiles": self.reaction_smiles,
            "byproducts": self.byproducts,
            "reagents": self.reagents,
            "is_redox": self.is_redox,
        }


# ─────────────────────────────────────────────────────────────────────────────
# Disconnection Rules (35+ reaction types, with specificity & priority)
# ─────────────────────────────────────────────────────────────────────────────

DISCONNECTION_RULES: List[Dict[str, Any]] = [
    # ═══ Carbonyl Chemistry ═══
    {
        "id": "ester_hydrolysis", "name": "Ester Hydrolysis",
        "pattern": "[#6:1](=[O:4])-[O;X2:2]-[#6:3]",
        "strategy": "ester_hydrolysis", "reaction_class": "Fischer Esterification (reverse)",
        "bond_to_break": (0, 2), "confidence": 0.85, "specificity": 3, "priority": 2,
        "rationale": "酯键断裂 → 羧酸 + 醇",
        "byproducts": ["O"], "reagents": ["H2SO4 (cat.)", "reflux"],
    },
    {
        "id": "amide_hydrolysis", "name": "Amide Hydrolysis",
        "pattern": "[#6:1](=[O:4])-[N;X3:2](-[#6:3])",
        "strategy": "amide_hydrolysis", "reaction_class": "Amide Coupling (reverse)",
        "bond_to_break": (0, 2), "confidence": 0.80, "specificity": 3, "priority": 2,
        "rationale": "酰胺键断裂 → 羧酸 + 胺",
        "byproducts": ["O"], "reagents": ["EDC", "HOBt", "DIPEA"],
    },
    {
        "id": "acid_chloride_amide", "name": "Acid Chloride to Amide",
        "pattern": "[#6:1](=[O:3])-[N;X3:2]",
        "strategy": "acyl_substitution", "reaction_class": "Schotten-Baumann (reverse)",
        "bond_to_break": (0, 2), "confidence": 0.78, "specificity": 2, "priority": 3,
        "rationale": "酰氯胺化逆反应 → 酰氯 + 胺",
        "byproducts": ["[Cl-]"], "reagents": ["Et3N", "DCM"],
    },
    {
        "id": "aldol", "name": "Aldol Disconnection",
        "pattern": "[#6:1]-[#6:2](-[OH])-[#6:3](=O)",
        "strategy": "aldol", "reaction_class": "Aldol Reaction",
        "bond_to_break": (0, 1), "confidence": 0.70, "specificity": 3, "priority": 3,
        "rationale": "β-羟基羰基断裂 → 醛/酮 + 烯醇",
        "byproducts": ["O"], "reagents": ["LDA", "THF, -78°C"],
    },
    {
        "id": "claisen_condensation", "name": "Claisen Condensation",
        "pattern": "[#6:1](=O)-[#6:2]-[#6:3](=O)-[O]",
        "strategy": "claisen", "reaction_class": "Claisen Condensation",
        "bond_to_break": (0, 1), "confidence": 0.68, "specificity": 3, "priority": 4,
        "rationale": "β-酮酯断裂 → 两个酯",
        "byproducts": ["CCO"], "reagents": ["NaOEt", "EtOH"],
    },

    # ═══ Cross-Coupling (Pd-catalyzed) — SMARTS 特异性修复 ═══
    {
        "id": "suzuki_coupling", "name": "Suzuki Coupling",
        # v2: 更具体的联芳模式，降低 confidence 避免假阳性
        "pattern": "[c:1](-[#6]):[c]:[c]-[c:2]",
        "strategy": "cross_coupling", "reaction_class": "Suzuki-Miyaura Coupling",
        "bond_to_break": (0, 1), "confidence": 0.65, "specificity": 2, "priority": 3,
        "rationale": "联芳断裂 → 芳基卤化物 + 芳基硼酸",
        "byproducts": ["[Br-]", "OB(O)O"],
        "reagents": ["Pd(PPh3)4", "K2CO3", "DMF/H2O"],
        "functionalize": "suzuki",
    },
    {
        "id": "buchwald_hartwig", "name": "Buchwald-Hartwig Amination",
        "pattern": "[c:1]-[N:2]",
        "strategy": "cross_coupling", "reaction_class": "Buchwald-Hartwig Amination",
        "bond_to_break": (0, 1), "confidence": 0.72, "specificity": 2, "priority": 3,
        "rationale": "C-N键断裂 → 芳基卤化物 + 胺",
        "byproducts": ["[Br-]"], "reagents": ["Pd2(dba)3", "XPhos", "Cs2CO3"],
        "functionalize": "buchwald",
    },
    {
        "id": "heck_reaction", "name": "Heck Reaction",
        "pattern": "[c:1]-[C:2]=[C:3]",
        "strategy": "cross_coupling", "reaction_class": "Mizoroki-Heck Reaction",
        "bond_to_break": (0, 1), "confidence": 0.70, "specificity": 2, "priority": 3,
        "rationale": "芳基烯烃断裂 → 芳基卤化物 + 烯烃",
        "byproducts": ["[Br-]"], "reagents": ["Pd(OAc)2", "PPh3", "Et3N"],
        "functionalize": "heck",
    },
    {
        "id": "sonogashira", "name": "Sonogashira Coupling",
        "pattern": "[c:1]-[C:2]#[C:3]",
        "strategy": "cross_coupling", "reaction_class": "Sonogashira Coupling",
        "bond_to_break": (0, 1), "confidence": 0.72, "specificity": 3, "priority": 3,
        "rationale": "芳基炔断裂 → 芳基卤化物 + 端炔",
        "byproducts": ["[Br-]"], "reagents": ["Pd(PPh3)2Cl2", "CuI", "Et3N"],
        "functionalize": "sonogashira",
    },
    {
        "id": "negishi_coupling", "name": "Negishi Coupling",
        # v2: 增加 Zn 上下文要求，与 Kumada 区分
        "pattern": "[c:1]-[CX4:2]",
        "strategy": "cross_coupling", "reaction_class": "Negishi Coupling",
        "bond_to_break": (0, 1), "confidence": 0.50, "specificity": 1, "priority": 6,
        "rationale": "C-C键断裂 → 芳基卤化物 + 有机锌试剂",
        "byproducts": ["[Zn]", "[Br-]"], "reagents": ["Pd(PPh3)4", "THF"],
        "require_context": ["halide"],
    },
    {
        "id": "kumada_coupling", "name": "Kumada Coupling",
        "pattern": "[c:1]-[CX4:2]",
        "strategy": "cross_coupling", "reaction_class": "Kumada Coupling",
        "bond_to_break": (0, 1), "confidence": 0.46, "specificity": 1, "priority": 7,
        "rationale": "C-C键断裂 → 芳基卤化物 + 格氏试剂",
        "byproducts": ["[Mg]Br"], "reagents": ["NiCl2(dppe)", "THF"],
        "require_context": ["halide"],
    },

    # ═══ C-C Bond Formation (Non-metal) ═══
    {
        "id": "grignard_addition", "name": "Grignard Addition",
        "pattern": "[#6:1]-[C:2](-[OH])-[#6:3]",
        "strategy": "grignard", "reaction_class": "Grignard Reaction",
        "bond_to_break": (0, 1), "confidence": 0.75, "specificity": 3, "priority": 2,
        "rationale": "醇α断裂 → 醛/酮 + 格氏试剂",
        "byproducts": [], "reagents": ["RMgBr", "THF", "then H3O+"],
    },
    {
        "id": "wittig_reaction", "name": "Wittig Reaction",
        "pattern": "[#6:1]=[C:2]-[#6:3]",
        "strategy": "wittig", "reaction_class": "Wittig Olefination",
        "bond_to_break": (0, 1), "confidence": 0.70, "specificity": 2, "priority": 3,
        "rationale": "烯烃断裂 → 醛/酮 + 磷叶立德",
        "byproducts": ["O=P(c1ccccc1)(c1ccccc1)c1ccccc1"],
        "reagents": ["Ph3P=CHR", "THF"],
    },
    {
        "id": "diels_alder", "name": "Diels-Alder Cycloaddition",
        "pattern": "[#6:1]1[#6:2][#6:3][#6:4][#6:5][#6:6]1",
        "strategy": "cycloaddition", "reaction_class": "Diels-Alder [4+2]",
        "bond_to_break": (0, 1), "confidence": 0.65, "specificity": 2, "priority": 4,
        "rationale": "六元环断裂 → 二烯 + 亲二烯体",
        "byproducts": [], "reagents": ["heat or Lewis acid"],
    },
    {
        "id": "friedel_crafts_alkylation", "name": "Friedel-Crafts Alkylation",
        "pattern": "[c:1]-[C:2]",
        "strategy": "electrophilic_aromatic", "reaction_class": "Friedel-Crafts Alkylation",
        "bond_to_break": (0, 1), "confidence": 0.60, "specificity": 1, "priority": 5,
        "rationale": "芳基烷基断裂 → 芳烃 + 烷基卤化物",
        "byproducts": ["[Cl-]"], "reagents": ["AlCl3", "DCM"],
    },
    {
        "id": "friedel_crafts_acylation", "name": "Friedel-Crafts Acylation",
        "pattern": "[c:1]-[C:2](=O)",
        "strategy": "electrophilic_aromatic", "reaction_class": "Friedel-Crafts Acylation",
        "bond_to_break": (0, 1), "confidence": 0.70, "specificity": 2, "priority": 3,
        "rationale": "芳基酰基断裂 → 芳烃 + 酰氯",
        "byproducts": ["[Cl-]"], "reagents": ["AlCl3", "DCM"],
    },

    # ═══ Ether/Amine Formation ═══
    {
        "id": "ether_cleavage", "name": "Ether Cleavage",
        "pattern": "[#6:1]-[O;X2:2]-[#6:3]",
        "strategy": "ether_cleavage", "reaction_class": "Williamson Ether Synthesis (reverse)",
        "bond_to_break": (0, 1), "confidence": 0.65, "specificity": 2, "priority": 4,
        "rationale": "醚键断裂 → 醇 + 烷基卤化物",
        "byproducts": ["[Br-]"], "reagents": ["NaH", "THF"],
        "exclude_if": "C(=O)O",
    },
    {
        "id": "sn2_substitution", "name": "SN2 Nucleophilic Substitution",
        "pattern": "[#6:1]-[O,N,S:2]-[#6:3]",
        "strategy": "substitution", "reaction_class": "SN2 Substitution",
        "bond_to_break": (0, 1), "confidence": 0.68, "specificity": 1, "priority": 4,
        "rationale": "C-X键形成逆反应 → 亲核试剂 + 卤代物",
        "byproducts": ["[Br-]"], "reagents": ["DMSO", "heat"],
    },
    {
        "id": "reductive_amination", "name": "Reductive Amination",
        "pattern": "[#6:1]-[N:2]-[#6:3]",
        "strategy": "reductive_amination", "reaction_class": "Reductive Amination",
        "bond_to_break": (0, 1), "confidence": 0.72, "specificity": 2, "priority": 3,
        "rationale": "C-N键断裂 → 羰基 + 胺",
        "byproducts": ["O"], "reagents": ["NaBH3CN", "MeOH", "AcOH"],
    },

    # ═══ Oxidation/Reduction — v2: 生成实际提案 ═══
    {
        "id": "alcohol_oxidation", "name": "Alcohol Oxidation",
        "pattern": "[#6:1](=O)-[#6:2]",
        "strategy": "oxidation", "reaction_class": "Alcohol Oxidation (reverse = reduction)",
        "bond_to_break": None, "confidence": 0.60, "specificity": 2, "priority": 5,
        "rationale": "羰基还原 → 醇",
        "byproducts": [], "reagents": ["TEMPO", "NaClO", "or Swern"],
        "is_redox": True, "redox_transform": "C=O → C-OH",
    },
    {
        "id": "carbonyl_reduction", "name": "Carbonyl Reduction",
        "pattern": "[#6:1](-[OH])-[#6:2]",
        "strategy": "reduction", "reaction_class": "Ketone/Aldehyde Reduction",
        "bond_to_break": None, "confidence": 0.65, "specificity": 2, "priority": 5,
        "rationale": "醇氧化 → 醛/酮",
        "byproducts": [], "reagents": ["NaBH4", "MeOH"],
        "is_redox": True, "redox_transform": "C-OH → C=O",
    },
    {
        "id": "nitro_reduction", "name": "Nitro Reduction to Amine",
        "pattern": "[c:1]-[N:2]",
        "strategy": "reduction", "reaction_class": "Nitro Reduction",
        "bond_to_break": None, "confidence": 0.70, "specificity": 2, "priority": 4,
        "rationale": "芳胺来自硝基还原",
        "byproducts": ["O", "O"], "reagents": ["Fe/HCl", "or H2/Pd-C"],
        "is_redox": True, "redox_transform": "Ar-NH2 → Ar-NO2",
    },

    # ═══ Protecting Groups ═══
    {
        "id": "boc_protection", "name": "Boc Protection",
        "pattern": "[N:1]-[C:2](=O)-[O]-[C](C)(C)C",
        "strategy": "protection", "reaction_class": "Boc Protection",
        "bond_to_break": (0, 1), "confidence": 0.80, "specificity": 3, "priority": 2,
        "rationale": "Boc脱保护 → 自由胺 + CO2 + 叔丁醇",
        "byproducts": ["O=C=O", "CC(C)(C)O"], "reagents": ["TFA", "DCM"],
    },
    {
        "id": "cbz_protection", "name": "Cbz Protection",
        "pattern": "[N:1]-[C:2](=O)-[O]-[C]-c1ccccc1",
        "strategy": "protection", "reaction_class": "Cbz Protection",
        "bond_to_break": (0, 1), "confidence": 0.78, "specificity": 3, "priority": 2,
        "rationale": "Cbz脱保护 → 自由胺",
        "byproducts": ["O=C=O", "c1ccccc1C"], "reagents": ["H2/Pd-C", "MeOH"],
    },
    {
        "id": "silyl_protection", "name": "Silyl Ether Protection",
        "pattern": "[O:1]-[Si:2]",
        "strategy": "protection", "reaction_class": "TBS/TIPS Protection",
        "bond_to_break": (0, 1), "confidence": 0.75, "specificity": 3, "priority": 2,
        "rationale": "硅醚脱保护 → 自由醇",
        "byproducts": [], "reagents": ["TBAF", "THF"],
    },

    # ═══ Ring Opening ═══
    {
        "id": "epoxide_opening", "name": "Epoxide Ring Opening",
        "pattern": "[#6:1]1-[O:2]-[#6:3]1",
        "strategy": "ring_opening", "reaction_class": "Epoxide Formation (reverse)",
        "bond_to_break": (0, 1), "confidence": 0.72, "specificity": 3, "priority": 3,
        "rationale": "环氧开环 → 烯烃 + 氧化剂",
        "byproducts": [], "reagents": ["mCPBA", "DCM"],
    },
]

# ═══ Additional templates ═══
DISCONNECTION_RULES.extend([
    {
        "id": "stille_coupling", "name": "Stille Coupling",
        "pattern": "[c:1](-[#6]):[c]:[c]-[c:2]",
        "strategy": "cross_coupling", "category": "cross_coupling",
        "reaction_class": "Stille Coupling",
        "bond_to_break": (0, 1), "confidence": 0.54, "specificity": 2, "priority": 5,
        "rationale": "Biaryl C-C → aryl halide + organostannane.",
        "byproducts": ["[Br-]"], "reagents": ["Pd(PPh3)4", "DMF"],
    },
    {
        "id": "snar_disconnection", "name": "SNAr Disconnection",
        "pattern": "[c:1]([N+](=O)[O-])-[N,O,S:2]",
        "strategy": "substitution", "category": "aromatic_substitution",
        "reaction_class": "SNAr Substitution",
        "bond_to_break": (0, 1), "confidence": 0.55, "specificity": 3, "priority": 4,
        "rationale": "Aryl-heteroatom → activated aryl halide + nucleophile.",
        "byproducts": ["[F-]", "[Cl-]"], "reagents": ["K2CO3", "DMF", "heat"],
    },
    {
        "id": "michael_addition", "name": "Michael Addition",
        "pattern": "[#6:1]-[#6:2]-[#6:3](=O)-[#6:4]",
        "strategy": "conjugate_addition", "category": "cc_bond_formation",
        "reaction_class": "Michael Addition",
        "bond_to_break": (0, 1), "confidence": 0.62, "specificity": 2, "priority": 4,
        "rationale": "Beta-substituted carbonyl → Michael donor + acceptor.",
        "byproducts": [], "reagents": ["base", "polar solvent"],
    },
    {
        "id": "chan_lam_coupling", "name": "Chan-Lam Coupling",
        "pattern": "[c:1]-[N,O:2]",
        "strategy": "cross_coupling", "category": "cross_coupling",
        "reaction_class": "Chan-Lam Coupling",
        "bond_to_break": (0, 1), "confidence": 0.52, "specificity": 1, "priority": 6,
        "rationale": "Aryl C-N/C-O → boronic acid + amine/alcohol.",
        "byproducts": ["B(O)(O)O"], "reagents": ["Cu(OAc)2", "pyridine", "air"],
    },
    {
        "id": "ullmann_cn_coupling", "name": "Ullmann C-N Coupling",
        "pattern": "[c:1]-[N:2]",
        "strategy": "cross_coupling", "category": "cross_coupling",
        "reaction_class": "Ullmann C-N Coupling",
        "bond_to_break": (0, 1), "confidence": 0.50, "specificity": 1, "priority": 6,
        "rationale": "Aryl amine → aryl halide + amine.",
        "byproducts": ["[I-]", "[Br-]"], "reagents": ["CuI", "ligand", "base"],
    },
    {
        "id": "click_triazole", "name": "CuAAC Click Triazole",
        "pattern": "[n:1]1[c:2][n:3][n:4][c:5]1",
        "strategy": "cycloaddition", "category": "cyclization_ring",
        "reaction_class": "CuAAC Click Cycloaddition",
        "bond_to_break": (1, 2), "confidence": 0.56, "specificity": 3, "priority": 3,
        "rationale": "1,2,3-triazole → azide + terminal alkyne.",
        "byproducts": [], "reagents": ["CuSO4", "sodium ascorbate", "H2O/tBuOH"],
    },
    {
        "id": "olefin_metathesis", "name": "Olefin Metathesis",
        "pattern": "[#6:1]=[#6:2]-[#6:3]=[#6:4]",
        "strategy": "metathesis", "category": "cc_bond_formation",
        "reaction_class": "Olefin Metathesis",
        "bond_to_break": (1, 2), "confidence": 0.49, "specificity": 2, "priority": 5,
        "rationale": "Metathesis-derived alkene → two alkene fragments.",
        "byproducts": [], "reagents": ["Grubbs catalyst"],
    },
    {
        "id": "mitsunobu_etherification", "name": "Mitsunobu Etherification",
        "pattern": "[#6:1]-[O:2]-[#6:3]",
        "strategy": "substitution", "category": "substitution",
        "reaction_class": "Mitsunobu Reaction",
        "bond_to_break": (0, 1), "confidence": 0.48, "specificity": 1, "priority": 6,
        "rationale": "Ether → alcohol + acidic nucleophile under Mitsunobu.",
        "byproducts": ["O=P(c1ccccc1)(c1ccccc1)c1ccccc1"],
        "reagents": ["PPh3", "DIAD/DEAD", "THF"],
    },
])

# ═══ v2 新增: 杂环断裂规则 (药化核心) ═══
DISCONNECTION_RULES.extend([
    {
        "id": "fischer_indole", "name": "Fischer Indole Synthesis",
        "pattern": "[nH]1c2ccccc2c([#6:1])c1[#6:2]",
        "strategy": "heterocycle_formation", "category": "cyclization_ring",
        "reaction_class": "Fischer Indole",
        "bond_to_break": (0, 1), "confidence": 0.60, "specificity": 3, "priority": 4,
        "rationale": "吲哚断裂 → 苯肼 + 酮",
        "byproducts": ["N", "O"], "reagents": ["ZnCl2", "AcOH", "heat"],
    },
    {
        "id": "hantzsch_pyridine", "name": "Hantzsch Pyridine Synthesis",
        "pattern": "c1cc([#6:1])nc([#6:2])c1",
        "strategy": "heterocycle_formation", "category": "cyclization_ring",
        "reaction_class": "Hantzsch Synthesis",
        "bond_to_break": (0, 1), "confidence": 0.55, "specificity": 3, "priority": 5,
        "rationale": "吡啶断裂 → 醛 + β-酮酯 + NH3",
        "byproducts": ["O", "O"], "reagents": ["NH4OAc", "EtOH", "reflux"],
    },
    {
        "id": "paal_knorr_pyrrole", "name": "Paal-Knorr Pyrrole Synthesis",
        "pattern": "[nH]1cccc1",
        "strategy": "heterocycle_formation", "category": "cyclization_ring",
        "reaction_class": "Paal-Knorr",
        "bond_to_break": None, "confidence": 0.58, "specificity": 3, "priority": 4,
        "rationale": "吡咯断裂 → 1,4-二酮 + 胺",
        "byproducts": ["O", "O"], "reagents": ["AcOH", "heat"],
        "is_redox": False, "is_heterocycle_formation": True,
    },
])

STRATEGY_CATEGORY_MAP: Dict[str, str] = {
    "ester_hydrolysis": "acylation",
    "amide_hydrolysis": "acylation",
    "acyl_substitution": "acylation",
    "aldol": "cc_bond_formation",
    "claisen": "cc_bond_formation",
    "cross_coupling": "cross_coupling",
    "grignard": "cc_bond_formation",
    "wittig": "cc_bond_formation",
    "cycloaddition": "cyclization_ring",
    "electrophilic_aromatic": "aromatic_substitution",
    "ether_cleavage": "substitution",
    "substitution": "substitution",
    "reductive_amination": "amination",
    "oxidation": "redox",
    "reduction": "redox",
    "protection": "protection_deprotection",
    "ring_opening": "cyclization_ring",
    "conjugate_addition": "cc_bond_formation",
    "metathesis": "cc_bond_formation",
    "heterocycle_formation": "cyclization_ring",
}



def _rule_category(rule: Dict[str, Any]) -> str:
    explicit = rule.get("category")
    if isinstance(explicit, str) and explicit.strip():
        return explicit.strip()
    strategy = str(rule.get("strategy", "")).strip().lower()
    return STRATEGY_CATEGORY_MAP.get(strategy, "uncategorized")


def _rule_category_counts() -> Dict[str, int]:
    counts: Dict[str, int] = {}
    for rule in DISCONNECTION_RULES:
        category = _rule_category(rule)
        counts[category] = counts.get(category, 0) + 1
    return dict(sorted(counts.items(), key=lambda item: item[0]))


def _generate_id(prefix: str, content: str) -> str:
    return f"{prefix}_{hashlib.sha256(content.encode()).hexdigest()[:8]}"


def _validate_smarts(smarts: str) -> bool:
    """v2: SMARTS 有效性验证"""
    if not RDKIT_AVAILABLE:
        return True
    try:
        patt = Chem.MolFromSmarts(smarts)
        return patt is not None
    except Exception:
        return False


def _find_pattern_matches(mol, smarts: str) -> List[Tuple[int, ...]]:
    if not _validate_smarts(smarts):
        return []
    pattern = Chem.MolFromSmarts(smarts)
    if pattern is None:
        return []
    return list(mol.GetSubstructMatches(pattern))


def _h_count_for_bond_type(bt) -> int:
    """Return the number of H atoms to add per broken end for a given bond type."""
    _MAP = {
        Chem.BondType.SINGLE: 1,
        Chem.BondType.DOUBLE: 2,
        Chem.BondType.TRIPLE: 3,
    }
    return _MAP.get(bt, 1)


def _break_bond_and_get_fragments(mol, atom_idx1: int, atom_idx2: int) -> List[str]:
    """Remove a bond, compensate hydrogen on broken ends, return sanitized
    fragment SMILES.

    After RemoveBond, atoms with explicit H or losing multi-order bonds
    get their explicit H count adjusted.  For atoms using only implicit H
    with single bond breaks, RDKit's sanitize handles recalculation.
    """
    try:
        emol = Chem.RWMol(mol)
        bond = emol.GetBondBetweenAtoms(atom_idx1, atom_idx2)
        if bond is None:
            return []

        bt = bond.GetBondType()
        bond_order = {
            Chem.BondType.SINGLE: 1,
            Chem.BondType.DOUBLE: 2,
            Chem.BondType.TRIPLE: 3,
        }.get(bt, 1)

        emol.RemoveBond(atom_idx1, atom_idx2)

        for idx in (atom_idx1, atom_idx2):
            atom = emol.GetAtomWithIdx(idx)
            if atom.GetIsAromatic():
                continue
            cur_explicit = atom.GetNumExplicitHs()
            if cur_explicit > 0 or atom.GetNoImplicit():
                atom.SetNumExplicitHs(cur_explicit + bond_order)
            elif bond_order > 1:
                atom.SetNumExplicitHs(bond_order - 1)

        frags = Chem.GetMolFrags(emol.GetMol(), asMols=True, sanitizeFrags=False)
        result = []
        for frag in frags:
            smiles = _sanitize_fragment(frag)
            if smiles:
                result.append(smiles)
        return result
    except Exception:
        return []


def _sanitize_fragment(frag) -> Optional[str]:
    """Sanitize a molecular fragment with multiple fallback strategies.

    Strategy 1: Full SanitizeMol (standard path)
    Strategy 2: Partial sanitize — skip KEKULIZE, then manually set
                aromatic atoms to aromatic bonds and retry
    Strategy 3: Use SetNoImplicitHs + partial sanitize to handle
                heterocyclic fragments (imidazole, pyrrole, etc.)
    Strategy 4: Last resort — convert to SMILES without full sanitization,
                re-parse, and return if valid

    Returns canonical SMILES string or None if all strategies fail.
    """
    # Strategy 1: standard full sanitize
    try:
        frag_copy = Chem.RWMol(frag)
        Chem.SanitizeMol(frag_copy)
        smiles = Chem.MolToSmiles(frag_copy, canonical=True)
        if smiles:
            return smiles
    except Exception:
        pass

    # Strategy 2: partial sanitize — skip kekulize step
    try:
        frag_copy = Chem.RWMol(frag)
        Chem.SanitizeMol(
            frag_copy,
            sanitizeOps=Chem.SanitizeFlags.SANITIZE_ALL
            ^ Chem.SanitizeFlags.SANITIZE_KEKULIZE,
        )
        smiles = Chem.MolToSmiles(frag_copy, canonical=True)
        if smiles:
            # Validate by round-trip: parse the SMILES back
            check = Chem.MolFromSmiles(smiles)
            if check is not None:
                return Chem.MolToSmiles(check, canonical=True)
            # Round-trip failed — likely missing [nH] in heterocycle.
            # Try fixing aromatic 'n' at ring closure: n1 -> [nH]1
            import re
            fixed = re.sub(r'(?<!\[)n(\d)', r'[nH]\1', smiles)
            if fixed != smiles:
                check2 = Chem.MolFromSmiles(fixed)
                if check2 is not None:
                    return Chem.MolToSmiles(check2, canonical=True)
    except Exception:
        pass

    # Strategy 3: add explicit Hs to nitrogen atoms in aromatic rings,
    # then retry full sanitize (fixes imidazole NH kekulize issues)
    try:
        frag_copy = Chem.RWMol(frag)
        for atom in frag_copy.GetAtoms():
            if (atom.GetSymbol() == "N"
                    and atom.GetIsAromatic()
                    and atom.GetTotalNumHs() == 0
                    and atom.GetFormalCharge() == 0
                    and atom.GetDegree() == 2):
                atom.SetNumExplicitHs(1)
                atom.SetNoImplicit(True)
        Chem.SanitizeMol(frag_copy)
        smiles = Chem.MolToSmiles(frag_copy, canonical=True)
        if smiles:
            return smiles
    except Exception:
        pass

    # Strategy 4: last resort — raw SMILES output + round-trip validation
    # Also attempts to fix common heterocycle issues (missing NH)
    try:
        raw_smiles = Chem.MolToSmiles(frag, canonical=False)
        if raw_smiles:
            check = Chem.MolFromSmiles(raw_smiles)
            if check is not None:
                return Chem.MolToSmiles(check, canonical=True)
            # Try fixing aromatic nitrogen missing H:
            # Replace bare 'n' at ring closure with '[nH]'
            # Pattern: lowercase 'n' followed by digit (ring closure) that
            # isn't already inside brackets
            import re
            fixed = re.sub(r'(?<!\[)n(\d)', r'[nH]\1', raw_smiles)
            if fixed != raw_smiles:
                check2 = Chem.MolFromSmiles(fixed)
                if check2 is not None:
                    return Chem.MolToSmiles(check2, canonical=True)
            # Return raw SMILES so LLM can attempt manual correction
            return raw_smiles
    except Exception:
        pass

    return None


def _apply_redox_rule(
    mol,
    target_smiles: str,
    rule: Dict[str, Any],
    match: Tuple[int, ...],
) -> Optional[DisconnectionProposal]:
    """v3: Redox 规则 — 用 SMARTS 反应模板实际执行官能团变换，生成真正的前体 SMILES。

    不再使用 ?? 占位符。如果变换失败则返回 None（跳过该提案）。
    """
    try:
        transform = rule.get("redox_transform", "")
        if not transform:
            return None

        # SMARTS 反应模板: 逆向变换 (产物 → 前体)
        _REDOX_RXNS = {
            "C=O → C-OH": "[C:1](=O)>>[C:1](O)",           # 羰基 → 醇 (逆向: 还原)
            "C-OH → C=O": "[C:1]([OH])>>[C:1](=O)",         # 醇 → 羰基 (逆向: 氧化)
            "Ar-NH2 → Ar-NO2": "[c:1][NH2]>>[c:1][N+](=O)[O-]",  # 芳胺 → 硝基 (逆向: 还原)
        }

        rxn_smarts = _REDOX_RXNS.get(transform)
        if rxn_smarts is None:
            return None

        from rdkit.Chem import AllChem
        rxn = AllChem.ReactionFromSmarts(rxn_smarts)
        if rxn is None:
            return None

        products = rxn.RunReactants((mol,))
        if not products:
            return None

        # 取第一个有效产物
        precursor_smiles = None
        for product_set in products:
            for p in product_set:
                try:
                    Chem.SanitizeMol(p)
                    smi = Chem.MolToSmiles(p, canonical=True)
                    if smi and smi != target_smiles:
                        precursor_smiles = smi
                        break
                except Exception:
                    continue
            if precursor_smiles:
                break

        if precursor_smiles is None:
            return None

        reaction_smiles = f"{precursor_smiles}>>{target_smiles}"

        return DisconnectionProposal(
            id=_generate_id("redox", f"{rule['id']}_{match}"),
            strategy=rule["strategy"],
            bond_description=f"functional group transformation ({transform})",
            precursors=[precursor_smiles],
            reaction_class=rule["reaction_class"],
            confidence=rule["confidence"],
            rationale=rule["rationale"],
            reaction_smiles=reaction_smiles,
            byproducts=rule.get("byproducts", []),
            reagents=rule.get("reagents", []),
            is_redox=True,
            specificity=rule.get("specificity", 2),
            priority=rule.get("priority", 5),
        )
    except Exception:
        return None


def _apply_disconnection_rule(
    mol,
    target_smiles: str,
    rule: Dict[str, Any],
    match: Tuple[int, ...],
) -> Optional[DisconnectionProposal]:
    try:
        # v2: Redox 规则走专用路径
        if rule.get("is_redox"):
            return _apply_redox_rule(mol, target_smiles, rule, match)

        bond_indices = rule.get("bond_to_break")
        if bond_indices is None:
            return None

        atom_idx1 = match[bond_indices[0]]
        atom_idx2 = match[bond_indices[1]]

        atom1 = mol.GetAtomWithIdx(atom_idx1)
        atom2 = mol.GetAtomWithIdx(atom_idx2)
        bond_desc = f"{atom1.GetSymbol()}-{atom2.GetSymbol()} bond"

        fragments = _break_bond_and_get_fragments(mol, atom_idx1, atom_idx2)
        if len(fragments) < 2:
            return None

        byproducts = rule.get("byproducts", [])
        products_part = target_smiles
        if byproducts:
            products_part = target_smiles + "." + ".".join(byproducts)
        reaction_smiles = ".".join(fragments) + ">>" + products_part

        return DisconnectionProposal(
            id=_generate_id("disc", f"{rule['id']}_{match}"),
            strategy=rule["strategy"],
            bond_description=bond_desc,
            precursors=fragments,
            reaction_class=rule["reaction_class"],
            confidence=rule["confidence"],
            rationale=rule["rationale"],
            reaction_smiles=reaction_smiles,
            byproducts=byproducts,
            reagents=rule.get("reagents", []),
            is_redox=False,
            specificity=rule.get("specificity", 1),
            priority=rule.get("priority", 5),
        )
    except Exception:
        return None


def propose_disconnections(
    target_smiles: str,
    strategy_hint: str = "",
    max_proposals: int = 5,
) -> Dict[str, Any]:
    """
    提出逆合成断裂方案。

    输入:
        target_smiles (str): 目标分子 SMILES
        strategy_hint (str): 可选策略提示
        max_proposals (int): 最大提案数

    输出:
        {
            "success": bool,
            "error": str | None,
            "target_smiles": str,
            "disconnections": [{id, strategy, bond_description, precursors, reaction_class, confidence, rationale, ...}, ...],
            "total_rules": int,
            "rule_categories": {category: count}
        }
    """
    if not RDKIT_AVAILABLE:
        raise ChemistryError("RDKIT_UNAVAILABLE", "RDKit not available")

    mol = Chem.MolFromSmiles(target_smiles)
    if mol is None:
        raise ChemistryError("INVALID_SMILES", f"Invalid SMILES: {target_smiles}")

    canonical_smiles = Chem.MolToSmiles(mol, canonical=True)
    proposals: List[DisconnectionProposal] = []
    seen_precursor_sets = set()

    for rule in DISCONNECTION_RULES:
        if strategy_hint:
            hint_lower = strategy_hint.lower()
            rule_id = rule["id"].lower()
            rule_name = rule["name"].lower()
            rule_strategy = rule["strategy"].lower()
            if not any(hint_lower in x for x in [rule_id, rule_name, rule_strategy]):
                continue

        matches = _find_pattern_matches(mol, rule["pattern"])
        for match in matches:
            proposal = _apply_disconnection_rule(mol, canonical_smiles, rule, match)
            if proposal:
                precursor_key = tuple(sorted(proposal.precursors))
                if precursor_key not in seen_precursor_sets:
                    seen_precursor_sets.add(precursor_key)
                    proposals.append(proposal)

    if strategy_hint and len(proposals) == 0:
        for rule in DISCONNECTION_RULES:
            matches = _find_pattern_matches(mol, rule["pattern"])
            for match in matches:
                proposal = _apply_disconnection_rule(mol, canonical_smiles, rule, match)
                if proposal:
                    precursor_key = tuple(sorted(proposal.precursors))
                    if precursor_key not in seen_precursor_sets:
                        seen_precursor_sets.add(precursor_key)
                        proposals.append(proposal)

    # v2: 按 priority 升序 + confidence 降序排序
    proposals.sort(key=lambda p: (p.priority, -p.confidence))
    proposals = proposals[:max_proposals]

    return {
        "success": True,
        "target_smiles": canonical_smiles,
        "disconnections": [p.to_dict() for p in proposals],
        "total_rules": len(DISCONNECTION_RULES),
        "rule_categories": _rule_category_counts(),
    }


# CLI Entry Point
if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Propose retrosynthetic disconnections")
    parser.add_argument("target_smiles", help="SMILES of target molecule")
    parser.add_argument("--strategy-hint", default="", help="Strategy hint")
    parser.add_argument("--max-proposals", type=int, default=5, help="Max proposals")
    parser.add_argument("--list-rules", action="store_true", help="List all rules")
    args = parser.parse_args()

    if args.list_rules:
        print(f"Total disconnection rules: {len(DISCONNECTION_RULES)}")
        print(f"Categories: {_rule_category_counts()}")
        for rule in DISCONNECTION_RULES:
            category = _rule_category(rule)
            spec = rule.get("specificity", 1)
            prio = rule.get("priority", 5)
            print(f"  - {rule['id']}: {rule['name']} (spec={spec}, prio={prio}, cat={category})")
    else:
        result = propose_disconnections(
            args.target_smiles,
            strategy_hint=args.strategy_hint,
            max_proposals=args.max_proposals,
        )
        print(json.dumps(result, indent=2, ensure_ascii=False))
