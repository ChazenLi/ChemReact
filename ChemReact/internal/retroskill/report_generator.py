import os
from typing import List, Dict, Any


def _append_unique(items: List[str], value: Any) -> None:
    if not isinstance(value, str):
        return
    cleaned = value.strip()
    if cleaned and cleaned not in items:
        items.append(cleaned)


def _collect_precursors(route: Dict[str, Any]) -> List[str]:
    precursors: List[str] = []
    for value in route.get("precursors", []):
        _append_unique(precursors, value)
    for step in route.get("steps", []):
        if not isinstance(step, dict):
            continue
        for value in step.get("reactants", []):
            _append_unique(precursors, value)
        reaction_smiles = step.get("reaction_smiles", "")
        if isinstance(reaction_smiles, str) and ">>" in reaction_smiles:
            left = reaction_smiles.split(">>", 1)[0]
            for value in left.split("."):
                _append_unique(precursors, value)
    return precursors

def generate_report(target_smiles: str, 
                    strategy_analysis: Dict[str, Any],
                    top_routes: List[Dict[str, Any]],
                    image_paths: Dict[str, str],
                    output_file: str = "RETRO_REPORT.md",
                    rejected_routes: List[Dict[str, Any]] = None,
                    audit_summary: Dict[str, Any] = None,
                    all_routes_visualized: bool = False,
                    planner_loop_summary: Dict[str, Any] = None,
                    planner_loop_enabled: bool = False,
                    route_count_policy: Dict[str, Any] = None):
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
    md.append("\n## 2. Candidate Routes (Visualized)" if all_routes_visualized else "\n## 2. Recommended Routes")

    if audit_summary:
        md.append(
            f"- **Audit Gate**: strict={audit_summary.get('strict_mode', False)}, "
            f"recommended={audit_summary.get('recommended_routes', 0)}, "
            f"rejected={audit_summary.get('rejected_routes', 0)}"
        )
        counts = audit_summary.get("verdict_counts", {})
        if isinstance(counts, dict) and counts:
            md.append(f"- **Verdict Counts**: {counts}")
    loop = planner_loop_summary if isinstance(planner_loop_summary, dict) else {}
    rounds = loop.get("rounds", []) if isinstance(loop.get("rounds", []), list) else []
    final_round = loop.get("final_round")
    if final_round is None and rounds:
        last = rounds[-1]
        if isinstance(last, dict) and isinstance(last.get("round"), int):
            final_round = last.get("round")
    md.append(f"- **Repair Loop Enabled**: {bool(planner_loop_enabled)}")
    md.append(f"- **Final Result Round**: {final_round if final_round is not None else 'N/A'}")
    if isinstance(loop.get("max_iterations"), int):
        md.append(f"- **Max Repair Attempts**: {loop.get('max_iterations')}")
    if "reached_repair_limit" in loop:
        md.append(f"- **Reached Repair Limit**: {bool(loop.get('reached_repair_limit'))}")
    if isinstance(loop.get("termination_reason"), str) and loop.get("termination_reason"):
        md.append(f"- **Repair Loop Termination**: {loop.get('termination_reason')}")
    if isinstance(route_count_policy, dict):
        md.append(
            f"- **Route Count Policy**: required>={route_count_policy.get('minimum_required')}, "
            f"provided={route_count_policy.get('provided')} "
            f"(reason={route_count_policy.get('partial_allowed_due_to')})"
        )

    if not top_routes:
        md.append("*No routes passed recommendation gate.*")

    for i, route in enumerate(top_routes):
        route_id = route.get("route_id", i+1)
        score = route.get("score", 0)
        verdict = route.get("audit_verdict", "UNKNOWN")
        
        md.append(f"### Route {route_id} (Score: {score}/10 - {verdict})")
        precursors = _collect_precursors(route)
        md.append(f"- **Precursors (RS side)**: {', '.join(precursors) if precursors else 'N/A'}")
        
        # Route Overview Image
        route_img_key = f"route_{route_id}_overview"
        if route_img_key in image_paths:
            md.append(f"![Route {route_id}]({image_paths[route_img_key]})\n")
        p_idx = 1
        found_precursor_image = False
        while True:
            p_key = f"route_{route_id}_precursor_{p_idx}"
            if p_key not in image_paths:
                break
            found_precursor_image = True
            md.append(f"![Route {route_id} RS {p_idx}]({image_paths[p_key]})")
            p_idx += 1
        if found_precursor_image:
            md.append("")
        
        # Reaction Tree Image
        tree_img_key = f"route_{route_id}_tree"
        if tree_img_key in image_paths:
             md.append(f"![Route {route_id} Tree]({image_paths[tree_img_key]})\n")
        
        # Auditor Notes
        md.append("#### Auditor's Verdict")
        critical_issues = route.get("critical_issues", [])
        critical_issues_text = ", ".join(critical_issues) if critical_issues else "None"
        md.append(f"- **Critical Issues**: {critical_issues_text}")
        if "quality_score" in route:
            md.append(f"- **Quality Score**: {route.get('quality_score')}/10")
        issue_codes = route.get("audit_details", {}).get("issue_codes", [])
        if issue_codes:
            md.append(f"- **Issue Codes**: {', '.join(issue_codes)}")
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
            reaction_smiles = step.get("reaction_smiles")
            if isinstance(reaction_smiles, str) and reaction_smiles.strip():
                md.append(f"- **Reaction (RS>>PS)**: `{reaction_smiles}`")
            
            # Step Image
            step_img_key = f"route_{route_id}_step_{j+1}"
            if step_img_key in image_paths:
                md.append(f"![Step {j+1}]({image_paths[step_img_key]})")
            md.append("\n")

    rejected_routes = rejected_routes or []
    shown_ids = {route.get("route_id") for route in top_routes if isinstance(route, dict)}
    tail_rejected = [r for r in rejected_routes if r.get("route_id") not in shown_ids]
    if rejected_routes:
        md.append("## 3. Strict Audit Reminders" if all_routes_visualized else "## 3. Rejected Routes (Gate Blocked)")
        source = tail_rejected if all_routes_visualized else rejected_routes
        if all_routes_visualized and not source:
            source = rejected_routes
        for i, route in enumerate(source):
            route_id = route.get("route_id", i + 1)
            score = route.get("score", 0)
            verdict = route.get("audit_verdict", "FAIL")
            md.append(f"### Route {route_id} (Score: {score}/10 - {verdict})")
            critical_issues = route.get("critical_issues", [])
            critical_issues_text = ", ".join(critical_issues) if critical_issues else "None"
            md.append(f"- **Audit Reminders**: {critical_issues_text}" if all_routes_visualized else f"- **Block Reasons**: {critical_issues_text}")
            if "quality_score" in route:
                md.append(f"- **Quality Score**: {route.get('quality_score')}/10")

    # 4. Next Steps
    md.append("## 4. Recommended Next Steps")
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
