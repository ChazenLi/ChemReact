[PROJECT_REPORT.md](https://github.com/user-attachments/files/25139398/PROJECT_REPORT.md)
# ChemReact: Integrated Retrosynthesis Planning & Visualization System

## 1. Project Overview

**ChemReact** is a sophisticated chemoinformatics and retrosynthesis automation framework designed to bridge the gap between high-level chemical reasoning and executable laboratory procedures. By leveraging the advanced molecular manipulation capabilities of **RDKit** and a **Multi-Persona LLM Architecture**, ChemReact provides a closed-loop system for molecular design, auditing, and high-fidelity visualization.

## 2. Core Architecture: Competitive Multi-Persona Reasoning

Unlike monolithic planning systems, ChemReact utilizes a "checks and balances" approach through specialized agent personas, ensuring that every synthetic route is both strategically sound and tactically executable.

| Persona | Responsibility | Focus |
| :--- | :--- | :--- |
| **Top-Level Designer** | Strategic Planning | Skeletal disconnection, convergent vs. linear strategy, core ring system construction. |
| **Reaction Designer** | Tactical Execution | Detailed reagent selection, solvent/catalyst optimization, selectivity control (Regio/Stereo). |
| **Auditor** | Quality Assurance | Mass balance verification, protecting group (PG) loop detection, safety/toxicity screening. |
| **Visualization Specialist** | Visual Communication | Creative direction for molecular rendering, identifying key intermediates and reaction centers for highlighting. |

## 3. Visual Results and Demonstrations

The system's strength lies in its ability to translate abstract JSON data into intuitive visual artifacts. Below are representative outputs from a standard run (`run_001`) performing the retrosynthesis of **Losartan**.

### 3.1 Target Molecule Analysis
The system identifies the core biphenyl-tetrazole-imidazole scaffold, providing a clean 2D representation for strategic mapping.

<img width="400" height="400" alt="target" src="https://github.com/user-attachments/assets/9c853df2-d31c-41f4-bb1b-1ee1241fd8f7" />

*Fig 1. Target Molecule (Losartan) 2D Rendering.*

### 3.2 Retrosynthesis Tree (The "Tree View")
A cornerstone of the ChemReact system is the orthogonal tree visualization, which clearly maps the target to its primary precursors, illustrating the convergent nature of the synthesis.

<img width="1000" height="800" alt="route_1_tree" src="https://github.com/user-attachments/assets/a55bc067-6f43-461a-96f3-2c6a37f1544f" />
<img width="1100" height="800" alt="route_2_tree" src="https://github.com/user-attachments/assets/072ae0c2-3b4d-4234-998f-5bac2be572c3" />
*Fig 2. Route 1&2: Strategic Disconnection Tree.*

### 3.3 Reaction Step Detailing
For each step, the system generates high-fidelity reaction mappings, highlighting the transformation of functional groups and atomic changes.

<img width="800" height="300" alt="route_1_step_1" src="https://github.com/user-attachments/assets/9f857d75-b280-40fd-8086-3b1ea88207ae" />
<img width="800" height="300" alt="route_1_step_2" src="https://github.com/user-attachments/assets/d2672520-694f-49e0-8df5-3d77fe06dbcc" />

*Fig 3. Detailed Reaction Mapping: Step 1 (Coupling).*

<img width="800" height="300" alt="route_2_step_1" src="https://github.com/user-attachments/assets/09bc8cd9-4406-4a9d-bc9d-a74ad0dd258a" />
<img width="800" height="300" alt="route_2_step_2" src="https://github.com/user-attachments/assets/57f76fdd-954f-4b23-865f-1b0826283e50" />


*Fig 4. Detailed Reaction Mapping: Step 2 (Coupling).*
## 4. Key Technical Innovations

*   **Closed-Loop Verification**: Integration with `verify_skill.py` ensures that all RDKit-derived properties (LogP, MW, Fingerprints) are consistent throughout the planning process.
*   **Persona-Driven Prompting**: Specialized prompts in `prompts_personas.py` reduce hallucination by forcing agents to focus on their specific domain (e.g., the Auditor cannot ignore PG loops).
*   **Automated Reporting**: The `report_generator.py` compiles JSON audit trails and PNG assets into a single, cohesive Markdown document for peer review.

## 5. Conclusion

ChemReact transforms retrosynthesis from a solo "guessing" game into an audited, visual, and documented engineering process. By combining the precision of RDKit with the flexibility of multi-persona LLMs, it offers a scalable solution for early-stage drug discovery and process chemistry.

---
**Released Version**: v0.1.0      the multi-step skill is coming
[PROJECT_REPORT.pdf](https://github.com/user-attachments/files/25139404/PROJECT_REPORT.pdf)

