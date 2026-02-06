import os
from typing import List, Dict, Any

def generate_report(target_smiles: str, 
                    strategy_analysis: Dict[str, Any],
                    top_routes: List[Dict[str, Any]],
                    image_paths: Dict[str, str],
                    output_file: str = "RETRO_REPORT.md"):
    """
    Generates a Markdown report from the retrosynthesis analysis data.
    """
    
    md = []
    
    # 1. Header
    md.append(f"# Retrosynthesis Analysis Report")
    md.append(f"**Target Molecule**: `{target_smiles}`")
    if "target_image" in image_paths:
        md.append(f"![Target]({image_paths['target_image']})\n")
    
    # 2. Executive Summary (Strategist)
    analysis = strategy_analysis.get("analysis", {})
    complexity_features = analysis.get("complexity_features", [])
    complexity_text = ", ".join(complexity_features) if complexity_features else "N/A"
    md.append("## 1. Executive Summary (Deep Chemical Audit)")
    md.append(f"- **Core Skeleton**: {analysis.get('core_skeleton', 'N/A')}")
    md.append(f"- **Complexity**: {complexity_text}")
    md.append(f"- **Strategy**: **{analysis.get('strategy_type', 'Linear')}**")
    
    # 3. Top Routes
    md.append("\n## 2. Recommended Routes")
    
    for i, route in enumerate(top_routes):
        route_id = route.get("route_id", i+1)
        score = route.get("score", 0)
        verdict = route.get("audit_verdict", "UNKNOWN")
        
        md.append(f"### Route {route_id} (Score: {score}/10 - {verdict})")
        
        # Route Overview Image
        route_img_key = f"route_{route_id}_overview"
        if route_img_key in image_paths:
            md.append(f"![Route {route_id}]({image_paths[route_img_key]})\n")
        
        # Reaction Tree Image
        tree_img_key = f"route_{route_id}_tree"
        if tree_img_key in image_paths:
             md.append(f"![Route {route_id} Tree]({image_paths[tree_img_key]})\n")
        
        # Auditor Notes
        md.append("#### Auditor's Verdict")
        critical_issues = route.get("critical_issues", [])
        critical_issues_text = ", ".join(critical_issues) if critical_issues else "None"
        md.append(f"- **Critical Issues**: {critical_issues_text}")
        if "pg_audit" in route:
            pg = route["pg_audit"]
            md.append(f"- **Protection Groups**: {'Needed' if pg.get('needs_pg') else 'Not Needed'}")
            if pg.get("pg_suggestions"):
                md.append(f"  - Suggestions: {', '.join(pg['pg_suggestions'])}")
        
        # Steps Detail
        md.append("\n#### Detailed Steps")
        steps = route.get("steps", [])
        if not steps:
            md.append("*No detailed steps provided.*")
            
        for j, step in enumerate(steps):
            md.append(f"**Step {j+1}: {step.get('reaction_class', 'Reaction')}**")
            reagents = step.get("reagents", [])
            reagents_text = ", ".join(reagents) if reagents else "N/A"
            md.append(f"- **Reagents**: {reagents_text}")
            md.append(f"- **Conditions**: {step.get('conditions', 'N/A')}")
            
            # Step Image
            step_img_key = f"route_{route_id}_step_{j+1}"
            if step_img_key in image_paths:
                md.append(f"![Step {j+1}]({image_paths[step_img_key]})")
            md.append("\n")

    # 4. Next Steps
    md.append("## 3. Recommended Next Steps")
    md.append("1. Verify availability of Key Starting Materials (KSMs) for Route 1.")
    md.append("2. Run *Conformer Generation* (Module 4) on late-stage intermediates to check steric hindrance.")
    md.append("3. Review safety flags for Scale-up.")

    # Write to file
    out_dir = os.path.dirname(output_file)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    with open(output_file, "w", encoding="utf-8") as f:
        f.write("\n".join(md))
    
    print(f"Report generated: {output_file}")
    return output_file

if __name__ == "__main__":
    # Test Data
    mock_strategy = {
        "analysis": {
            "core_skeleton": "Biaryl-imidazole",
            "complexity_features": ["Tetrazole ring", "Biphenyl system"],
            "strategy_type": "Convergent"
        }
    }
    mock_routes = [
        {
            "route_id": 1,
            "score": 8.5,
            "audit_verdict": "PASS",
            "critical_issues": ["Tetrazole formation requires Azide"],
            "pg_audit": {"needs_pg": True, "pg_suggestions": ["Trityl for Tetrazole"]},
            "steps": [
                {"reaction_class": "Suzuki Coupling", "reagents": ["Pd(dppf)Cl2", "K2CO3"], "conditions": "DME/H2O, 80C"},
                {"reaction_class": "Deprotection", "reagents": ["HCl/MeOH"], "conditions": "RT, 2h"}
            ]
        }
    ]
    mock_images = {
        "target_image": "test_mol.png" 
    }
    
    generate_report("CCCCC1=NC(Cl)=...", mock_strategy, mock_routes, mock_images, "TEST_REPORT.md")
