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
## 4. Another Visual Results and Demonstrations 
**Target Molecule**: `Ic1ccc(c(c1)F)Nc1c(ccc(c1F)F)C(=O)N1CC(O)(C1)[C@@H]1CCCCN1`
<img width="400" height="400" alt="target" src="https://github.com/user-attachments/assets/79943cc5-a2d5-411d-abed-a312b8246b88" />
### 4.1 Executive Summary (Deep Chemical Audit)
- **Core Skeleton**: Biaryl-amide with fused azabicyclic system
- **Complexity**: Halogenated aryl rings (I, F), Chiral azabicyclo[2.2.1]heptane, Pyrrolidine moiety, Multiple amide linkages
- **Strategy**: **Convergent**

### 4.2 Recommended Routes
#### Route 1 (Score: 8.5/10 - PASS)
<img width="900" height="300" alt="route_1_overview" src="https://github.com/user-attachments/assets/f630825f-6bec-4ae8-90e6-51cf83def5d3" />

<img width="1000" height="800" alt="route_1_tree" src="https://github.com/user-attachments/assets/fa341cbe-3f39-444a-b529-e2d7d40bce8d" />


##### Auditor's Verdict
- **Critical Issues**: None

##### Detailed Steps
**Step 1: Nucleophilic Aromatic Substitution**
- **Reagents**: K2CO3, DMF
- **Conditions**: 80°C, 12h
<img width="800" height="300" alt="route_1_step_1" src="https://github.com/user-attachments/assets/90eb17fa-5af5-4f6d-948e-107f1a9761bb" />


**Step 2: Amide Coupling**
- **Reagents**: HATU, DIPEA, DMF
- **Conditions**: RT, 4h
<img width="800" height="300" alt="route_1_step_2" src="https://github.com/user-attachments/assets/231282a0-efd7-463f-9e6f-f8b2f87d1f50" />



#### Route 2 (Score: 7.8/10 - PASS)
<img width="900" height="300" alt="route_2_overview" src="https://github.com/user-attachments/assets/9a352809-a064-4a63-96d6-e977d797f940" />


<img width="1000" height="800" alt="route_2_tree" src="https://github.com/user-attachments/assets/c1b02b47-aa5c-4618-98da-c3bb371c21fc" />


##### Auditor's Verdict                  
- **Critical Issues**: Requires protection/deprotection                

#### Detailed Steps
**Step 1: Buchwald-Hartwig Coupling**                       
- **Reagents**: Pd2(dba)3, XPhos, Cs2CO3, Toluene                       
- **Conditions**: 100°C, 16h                      
<img width="800" height="300" alt="route_2_step_1" src="https://github.com/user-attachments/assets/16a31cac-823a-4bd3-8752-e9224c3754a6" />



### 4.3 Recommended Next Steps
1. Verify availability of Key Starting Materials (KSMs) for Route 1.                          
2. Run *Conformer Generation* (Module 4) on late-stage intermediates to check steric hindrance.                     
3. Review safety flags for Scale-up.              

## 5. Key Technical Innovations

*   **Closed-Loop Verification**: Integration with `verify_skill.py` ensures that all RDKit-derived properties (LogP, MW, Fingerprints) are consistent throughout the planning process.
*   **Persona-Driven Prompting**: Specialized prompts in `prompts_personas.py` reduce hallucination by forcing agents to focus on their specific domain (e.g., the Auditor cannot ignore PG loops).
*   **Automated Reporting**: The `report_generator.py` compiles JSON audit trails and PNG assets into a single, cohesive Markdown document for peer review.

## 6. Conclusion

ChemReact transforms retrosynthesis from a solo "guessing" game into an audited, visual, and documented engineering process. By combining the precision of RDKit with the flexibility of multi-persona LLMs, it offers a scalable solution for early-stage drug discovery and process chemistry.

---
**Released Version**: v0.1.0      the multi-step skill is coming
[PROJECT_REPORT.pdf](https://github.com/user-attachments/files/25139404/PROJECT_REPORT.pdf)

