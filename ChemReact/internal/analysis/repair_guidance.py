"""
Repair Guidance Module
======================
Provides intelligent repair suggestions based on issue codes from validation.

核心功能：根据审计失败的 issue code 生成具体的修复建议
"""

from typing import Dict, Any, List


# Issue code -> repair suggestions mapping
REPAIR_TEMPLATES: Dict[str, Dict[str, Any]] = {
    "ATOM_CONSERVATION_FAIL": {
        "priority": 1,
        "category": "atom_balance",
        "suggestions": [
            "Check for missing byproducts (common: H2O, HCl, HBr, NH3, CO2)",
            "Verify all reactants are included in the reaction",
            "Consider implicit hydrogen atoms",
            "For condensation reactions, add H2O as byproduct"
        ],
        "examples": [
            "Esterification: add H2O to products",
            "Amide coupling: add H2O to products",
            "SN2: add leaving group (Br-, Cl-) to products"
        ]
    },
    "ATOM_BALANCE_INFERRED_BYPRODUCT": {
        "priority": 2,
        "category": "atom_balance",
        "suggestions": [
            "Verify the inferred byproduct is chemically reasonable",
            "Consider adding the byproduct explicitly to the reaction"
        ]
    },
    "REACTION_FORMAT_INVALID": {
        "priority": 0,
        "category": "format",
        "suggestions": [
            "Ensure reaction format is: Reactants>>Products",
            "Check that '>>' separator is present",
            "Separate multiple reactants/products with '.' (dot)"
        ]
    },
    "REACTION_SMILES_INVALID": {
        "priority": 0,
        "category": "format",
        "suggestions": [
            "Verify all SMILES strings are valid",
            "Check for unclosed brackets or parentheses",
            "Verify atom symbols and bond notations"
        ]
    },
    "MAPPING_INVALID": {
        "priority": 2,
        "category": "mapping",
        "suggestions": [
            "Check atom mapping numbers for consistency",
            "Ensure mapped atoms exist on both sides",
            "Remove duplicate mapping numbers"
        ]
    },
    "MAPPING_INCOMPLETE": {
        "priority": 3,
        "category": "mapping",
        "suggestions": [
            "Add atom mapping to track bond changes",
            "Map reaction center atoms at minimum",
            "Consider using RXNMapper for automatic mapping"
        ]
    },
    "NO_BOND_CHANGE_DETECTED": {
        "priority": 2,
        "category": "bond",
        "suggestions": [
            "Verify reaction actually involves bond changes",
            "Check atom mapping is correct",
            "This may indicate identity reaction or error"
        ]
    },
    "BOND_CHANGE_INCONSISTENT": {
        "priority": 2,
        "category": "bond",
        "suggestions": [
            "Verify bond changes match reaction type",
            "For SN2: one C-X bond breaks, one C-Nu bond forms",
            "For amide coupling: one N-H breaks, one C-N forms"
        ]
    },
    "ROUTE_COHERENCE_FAIL": {
        "priority": 1,
        "category": "route",
        "suggestions": [
            "Verify each step's product is used in subsequent steps",
            "Check that step order is logical",
            "Ensure no orphan intermediates"
        ]
    },
    "TARGET_MISMATCH": {
        "priority": 0,
        "category": "route",
        "suggestions": [
            "Verify final step product matches target molecule",
            "Check for stereochemistry differences",
            "Verify canonical SMILES representation"
        ]
    },
    "RING_TOPOLOGY_MISMATCH": {
        "priority": 2,
        "category": "ring",
        "suggestions": [
            "Check ring opening/closing steps are chemically valid",
            "Verify ring sizes are consistent with target",
            "For cyclization: ensure appropriate ring size (5, 6 preferred)"
        ]
    },
    "HETEROCYCLE_TOPOLOGY_MISMATCH": {
        "priority": 2,
        "category": "ring",
        "suggestions": [
            "Verify heterocycle formation is chemically feasible",
            "Check heteroatom positions in ring",
            "Consider alternative cyclization approaches"
        ]
    },
    "SCAFFOLD_CONTINUITY_FAIL": {
        "priority": 1,
        "category": "scaffold",
        "suggestions": [
            "Verify key scaffold is preserved through synthesis",
            "Check for unintended scaffold modifications",
            "Identify steps that alter core structure"
        ]
    },
    "STRATEGY_MISMATCH": {
        "priority": 3,
        "category": "strategy",
        "suggestions": [
            "Align route steps with specified strategy",
            "Check if convergent/linear pattern is maintained",
            "Verify key disconnections match strategy"
        ]
    },
    "RISK_BUDGET_EXCEEDED": {
        "priority": 2,
        "category": "strategy",
        "suggestions": [
            "Reduce number of high-risk steps",
            "Consider alternative reactions with lower risk",
            "Add protective measures for risky transformations"
        ]
    },
    "PRECURSOR_ANALYSIS_FAIL": {
        "priority": 2,
        "category": "precursor",
        "suggestions": [
            "Verify all precursor SMILES are valid",
            "Check precursor structure compatibility",
            "Consider alternative precursor choices"
        ]
    },
    "COMPETING_NUCLEOPHILES": {
        "priority": 3,
        "category": "selectivity",
        "suggestions": [
            "Consider protecting one nucleophile temporarily",
            "Use more selective reaction conditions",
            "Adjust pH or temperature for chemoselectivity"
        ]
    },
    "ACID_BASE_COEXIST": {
        "priority": 3,
        "category": "selectivity",
        "suggestions": [
            "Consider salt form control",
            "Use protecting groups for one functionality",
            "Adjust reaction order to avoid premature reaction"
        ]
    }
}


def generate_repair_guidance(issues: List[Dict[str, Any]]) -> Dict[str, Any]:
    """
    Generate repair guidance based on validation issues.
    
    Args:
        issues: List of issues from validation, each with 'code' and 'message'
        
    Returns:
        {
            "guidance": [
                {
                    "issue_code": str,
                    "issue_message": str,
                    "priority": int,
                    "suggestions": [...],
                    "examples": [...]
                }
            ],
            "categories_affected": [...],
            "repair_strategy": str
        }
    """
    guidance = []
    categories = set()
    
    for issue in issues:
        code = issue.get("code", "UNKNOWN")
        message = issue.get("message", "")
        
        template = REPAIR_TEMPLATES.get(code, {
            "priority": 5,
            "category": "unknown",
            "suggestions": [f"Review and fix: {message}"]
        })
        
        guidance.append({
            "issue_code": code,
            "issue_message": message,
            "priority": template["priority"],
            "category": template["category"],
            "suggestions": template.get("suggestions", []),
            "examples": template.get("examples", [])
        })
        
        categories.add(template.get("category", "unknown"))
    
    # Sort by priority (lower = more critical)
    guidance.sort(key=lambda x: x["priority"])
    
    # Generate overall repair strategy
    strategy = _generate_repair_strategy(guidance, categories)
    
    return {
        "guidance": guidance,
        "categories_affected": list(categories),
        "repair_strategy": strategy,
        "total_issues": len(issues),
        "critical_issues": sum(1 for g in guidance if g["priority"] <= 1)
    }


def _generate_repair_strategy(
    guidance: List[Dict[str, Any]],
    categories: set
) -> str:
    """Generate an overall repair strategy description."""
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
    target_smiles: str
) -> str:
    """
    Build a repair prompt for the LLM with specific guidance.
    
    Args:
        rejected_routes: List of routes that failed validation
        strategy: Original strategy document
        target_smiles: Target molecule SMILES
        
    Returns:
        Formatted repair prompt string
    """
    lines = [
        "# Synthesis Route Repair Request",
        "",
        f"**Target Molecule**: `{target_smiles}`",
        "",
        "## Rejected Routes and Issues",
        ""
    ]
    
    for route in rejected_routes:
        route_id = route.get("route_id", "?")
        audit_details = route.get("audit_details", {})
        issues = audit_details.get("issues", [])
        
        lines.append(f"### Route {route_id}")
        
        # Generate guidance for this route's issues
        guidance_result = generate_repair_guidance(issues)
        
        lines.append(f"**Strategy**: {guidance_result['repair_strategy']}")
        lines.append("")
        
        for g in guidance_result["guidance"]:
            lines.append(f"- **{g['issue_code']}**: {g['issue_message']}")
            for suggestion in g["suggestions"][:2]:  # Limit to 2 suggestions
                lines.append(f"  - {suggestion}")
        
        lines.append("")
    
    lines.extend([
        "## Repair Instructions",
        "",
        "1. Address all critical issues (priority 0-1) first",
        "2. Fix atom balance by adding appropriate byproducts",
        "3. Ensure all steps connect logically",
        "4. Verify final product matches target exactly",
        "",
        "## Output Format",
        "",
        "Return repaired routes in JSON format with corrected reaction SMILES.",
        "Each step should have validated atom balance."
    ])
    
    return "\n".join(lines)


# CLI 入口
if __name__ == "__main__":
    import json
    import argparse
    
    parser = argparse.ArgumentParser(description="Generate repair guidance")
    parser.add_argument("--demo", action="store_true", help="Run demo")
    args = parser.parse_args()
    
    if args.demo:
        # Demo with sample issues
        sample_issues = [
            {"code": "ATOM_CONSERVATION_FAIL", "message": "Delta: {H: 2, O: 1}"},
            {"code": "COMPETING_NUCLEOPHILES", "message": "Multiple nucleophiles detected"}
        ]
        
        result = generate_repair_guidance(sample_issues)
        print(json.dumps(result, indent=2, ensure_ascii=False))
