"""
Repair Advisor Module
=====================
Intelligent repair suggestions based on validation issue codes.
Migrated from internal/analysis/repair_guidance.py.

Bug fix (Req 10.2): priority field is always converted to int before sorting,
preventing mixed str/int comparison errors.

Public API:
    generate_repair_guidance(issues) -> Dict[str, Any]
    build_repair_prompt(rejected_routes, strategy, target_smiles) -> str
"""

from __future__ import annotations

from typing import Any, Dict, List


# 26 issue code -> repair template mapping
REPAIR_TEMPLATES: Dict[str, Dict[str, Any]] = {
    "ATOM_CONSERVATION_FAIL": {
        "priority": 1, "category": "atom_balance",
        "suggestions": [
            "Check for missing byproducts (H2O, HCl, HBr, NH3, CO2, NaBr, LiCl)",
            "Verify all reactants are included",
            "Consider implicit hydrogen atoms",
            "For condensation reactions, add H2O as byproduct",
        ],
        "examples": [
            "Esterification: add H2O to products",
            "Amide coupling: add H2O to products",
            "SN2: add leaving group (Br-, Cl-) to products",
        ],
    },
    "ATOM_BALANCE_INFERRED_BYPRODUCT": {
        "priority": 2, "category": "atom_balance",
        "suggestions": [
            "Verify the inferred byproduct is chemically reasonable",
            "Consider adding the byproduct explicitly",
        ],
    },
    "REACTION_FORMAT_INVALID": {
        "priority": 0, "category": "format",
        "suggestions": [
            "Ensure format is: Reactants>>Products",
            "Check '>>' separator is present",
            "Separate multiple reactants/products with '.'",
        ],
    },
    "REACTION_SMILES_INVALID": {
        "priority": 0, "category": "format",
        "suggestions": [
            "Verify all SMILES strings are valid",
            "Check for unclosed brackets or parentheses",
        ],
    },
    "MAPPING_INVALID": {
        "priority": 2, "category": "mapping",
        "suggestions": [
            "Check atom mapping numbers for consistency",
            "Ensure mapped atoms exist on both sides",
        ],
    },
    "MAPPING_INCOMPLETE": {
        "priority": 3, "category": "mapping",
        "suggestions": [
            "Add atom mapping to track bond changes",
            "Consider using RXNMapper for automatic mapping",
        ],
    },
    "NO_BOND_CHANGE_DETECTED": {
        "priority": 2, "category": "bond",
        "suggestions": ["Verify reaction involves bond changes", "Check atom mapping"],
    },
    "BOND_CHANGE_INCONSISTENT": {
        "priority": 2, "category": "bond",
        "suggestions": [
            "Verify bond changes match reaction type",
            "For SN2: one C-X breaks, one C-Nu forms",
        ],
    },
    "ROUTE_COHERENCE_FAIL": {
        "priority": 1, "category": "route",
        "suggestions": [
            "Verify each step's product is used in subsequent steps",
            "Check step order is logical",
        ],
    },
    "TARGET_MISMATCH": {
        "priority": 0, "category": "route",
        "suggestions": [
            "Verify final step product matches target molecule",
            "Check for stereochemistry differences",
        ],
    },
    "RING_TOPOLOGY_MISMATCH": {
        "priority": 2, "category": "ring",
        "suggestions": ["Check ring opening/closing steps", "Verify ring sizes"],
    },
    "HETEROCYCLE_TOPOLOGY_MISMATCH": {
        "priority": 2, "category": "ring",
        "suggestions": ["Verify heterocycle formation is feasible", "Check heteroatom positions"],
    },
    "SCAFFOLD_CONTINUITY_FAIL": {
        "priority": 1, "category": "scaffold",
        "suggestions": ["Verify key scaffold is preserved", "Identify steps that alter core structure"],
    },
    "STRATEGY_MISMATCH": {
        "priority": 3, "category": "strategy",
        "suggestions": ["Align route steps with specified strategy"],
    },
    "RISK_BUDGET_EXCEEDED": {
        "priority": 2, "category": "strategy",
        "suggestions": ["Reduce high-risk steps", "Consider alternative reactions"],
    },
    "PRECURSOR_ANALYSIS_FAIL": {
        "priority": 2, "category": "precursor",
        "suggestions": ["Verify all precursor SMILES are valid"],
    },
    "COMPETING_NUCLEOPHILES": {
        "priority": 3, "category": "selectivity",
        "suggestions": ["Consider protecting one nucleophile", "Use selective conditions"],
    },
    "ACID_BASE_COEXIST": {
        "priority": 3, "category": "selectivity",
        "suggestions": ["Consider salt form control", "Use protecting groups"],
    },
    "REDOX_INCOMPATIBLE": {
        "priority": 1, "category": "compatibility",
        "suggestions": ["Separate oxidizing and reducing steps", "Use protecting groups"],
    },
    "AZIDE_SAFETY": {
        "priority": 0, "category": "safety",
        "suggestions": ["Handle azides with extreme care - explosion hazard"],
    },
    "ACYL_HALIDE_HYDROLYSIS": {
        "priority": 2, "category": "compatibility",
        "suggestions": ["Use anhydrous conditions"],
    },
    "EPOXIDE_NUCLEOPHILE_RISK": {
        "priority": 2, "category": "compatibility",
        "suggestions": ["Control temperature", "Use protecting groups on nucleophiles"],
    },
    "CATALYST_POISONING_RISK": {
        "priority": 2, "category": "compatibility",
        "suggestions": ["Remove or protect thiols before Pd-catalyzed reactions"],
    },
    "BASE_SENSITIVE_ESTER": {
        "priority": 3, "category": "compatibility",
        "suggestions": ["Use mild bases or buffered conditions"],
    },
    "ALDEHYDE_AMINE_IMINE": {
        "priority": 2, "category": "compatibility",
        "suggestions": ["Control stoichiometry", "Protect aldehyde or amine if imine undesired"],
    },
}


def generate_repair_guidance(issues: List[Dict[str, Any]]) -> Dict[str, Any]:
    """Generate repair guidance based on validation issues.

    Bug fix: priority is always cast to int() before sorting to prevent
    mixed-type comparison errors.

    Args:
        issues: List of {code, message} dicts from validation.

    Returns:
        {
            "guidance": [...],
            "categories_affected": [...],
            "repair_strategy": str,
            "total_issues": int,
            "critical_issues": int,
        }
    """
    guidance: List[Dict[str, Any]] = []
    categories: set = set()

    for issue in issues:
        code = issue.get("code", "UNKNOWN")
        message = issue.get("message", "")

        template = REPAIR_TEMPLATES.get(code, {
            "priority": 5, "category": "unknown",
            "suggestions": [f"Review and fix: {message}"],
        })

        # Bug fix (Req 10.2): ensure priority is always int
        priority = template.get("priority", 5)
        try:
            priority = int(priority)
        except (TypeError, ValueError):
            priority = 5

        guidance.append({
            "issue_code": code,
            "issue_message": message,
            "priority": priority,
            "category": template.get("category", "unknown"),
            "suggestions": template.get("suggestions", []),
            "examples": template.get("examples", []),
        })
        categories.add(template.get("category", "unknown"))

    # Sort by priority (int, lower = more critical)
    guidance.sort(key=lambda x: x["priority"])

    strategy = _generate_repair_strategy(guidance, categories)

    return {
        "guidance": guidance,
        "categories_affected": list(categories),
        "repair_strategy": strategy,
        "total_issues": len(issues),
        "critical_issues": sum(1 for g in guidance if g["priority"] <= 1),
    }


def _generate_repair_strategy(
    guidance: List[Dict[str, Any]], categories: set
) -> str:
    if not guidance:
        return "No issues detected. Route is ready for execution."

    critical = [g for g in guidance if g["priority"] <= 1]
    if critical:
        codes = [g["issue_code"] for g in critical]
        if "ATOM_CONSERVATION_FAIL" in codes:
            return "Critical: Fix atom balance by adding missing byproducts or reactants."
        if "TARGET_MISMATCH" in codes:
            return "Critical: Final product does not match target. Verify last step."
        if "REACTION_FORMAT_INVALID" in codes:
            return "Critical: Fix reaction SMILES format before further analysis."

    if "selectivity" in categories:
        return "Moderate: Address selectivity concerns with protecting groups or condition optimization."
    if "mapping" in categories:
        return "Optional: Add atom mapping for better analysis quality."

    return "Minor issues detected. Review suggestions and apply fixes."


def build_repair_prompt(
    rejected_routes: List[Dict[str, Any]],
    strategy: Dict[str, Any],
    target_smiles: str,
) -> str:
    """Build a repair prompt for the LLM with specific guidance."""
    lines = [
        "# Synthesis Route Repair Request",
        "",
        f"**Target Molecule**: `{target_smiles}`",
        "",
        "## Rejected Routes and Issues",
        "",
    ]

    for route in rejected_routes:
        route_id = route.get("route_id", "?")
        issues = route.get("audit_details", {}).get("issues", [])
        lines.append(f"### Route {route_id}")

        result = generate_repair_guidance(issues)
        lines.append(f"**Strategy**: {result['repair_strategy']}")
        lines.append("")

        for g in result["guidance"]:
            lines.append(f"- **{g['issue_code']}**: {g['issue_message']}")
            for s in g["suggestions"][:2]:
                lines.append(f"  - {s}")
        lines.append("")

    lines.extend([
        "## Repair Instructions",
        "",
        "1. Address all critical issues (priority 0-1) first",
        "2. Fix atom balance by adding appropriate byproducts",
        "3. Ensure all steps connect logically",
        "4. Verify final product matches target exactly",
    ])

    return "\n".join(lines)
