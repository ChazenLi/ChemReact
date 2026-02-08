"""Simple host planner that generates pre-defined routes."""
import argparse
import json
from pathlib import Path

# Hard-coded routes for the target molecule
ROUTES_DATA = [
  {
    "route_id": 1,
    "synthesis_style": "convergent",
    "score": 4.8,
    "audit_verdict": "CONDITIONAL",
    "critical_issues": ["Planner output requires validation."],
    "precursors": ["CCCCCc1nc(Cl)c(CO)n1", "Cc1ccc(-c2ccccc2-c2n[nH]nn2)cc1"],
    "steps": [
      {
        "reaction_class": "N-Alkylation",
        "reagents": ["Cs2CO3", "CH3CN"],
        "conditions": "reflux, 12h",
        "reaction_smiles": "CCCCCc1nc(Cl)c(CO)n1.Cc1ccc(-c2ccccc2-c2n[nH]nn2)cc1>>CCCCCc1nc(Cl)c(CO)n1Cc1ccc(-c2ccccc2-c2n[nH]nn2)cc1"
      },
      {
        "reaction_class": "Reduction",
        "reagents": ["NaBH4", "MeOH"],
        "conditions": "0°C to rt, 2h",
        "reaction_smiles": "CCCCCc1nc(Cl)c(C=O)n1>>CCCCCc1nc(Cl)c(CO)n1"
      },
      {
        "reaction_class": "[3+2] Cycloaddition",
        "reagents": ["NaN3", "NH4Cl", "DMF"],
        "conditions": "100°C, 18h",
        "reaction_smiles": "Cc1ccc(-c2ccccc2C#N)cc1.NaN3>>Cc1ccc(-c2ccccc2-c2n[nH]nn2)cc1.Na"
      },
      {
        "reaction_class": "Suzuki Coupling",
        "reagents": ["Pd(PPh3)4", "K2CO3"],
        "conditions": "90°C, toluene/H2O, 12h",
        "reaction_smiles": "Brc1ccc(C)cc1.Brc2ccccc2>>Cc1ccc(-c2ccccc2)cc1.PdBr2"
      }
    ]
  },
  {
    "route_id": 2,
    "synthesis_style": "linear",
    "score": 4.5,
    "audit_verdict": "CONDITIONAL",
    "critical_issues": ["Planner output requires validation."],
    "precursors": ["CCCCCc1nc(Cl)c(CO)n1", "Brc1ccc(C#N)cc1", "Brc2ccccc2B(OH)2"],
    "steps": [
      {
        "reaction_class": "Reduction",
        "reagents": ["NaBH4", "MeOH"],
        "conditions": "0°C to rt, 2h",
        "reaction_smiles": "CCCCCc1nc(Cl)c(C=O)n1>>CCCCCc1nc(Cl)c(CO)n1"
      },
      {
        "reaction_class": "N-Alkylation",
        "reagents": ["K2CO3", "CH3CN"],
        "conditions": "reflux, 10h",
        "reaction_smiles": "CCCCCc1nc(Cl)c(CO)n1.Brc1ccc(C#N)cc1>>CCCCCc1nc(Cl)c(CO)n1Cc1ccc(C#N)cc1.KBr"
      },
      {
        "reaction_class": "Suzuki Coupling",
        "reagents": ["Pd(dppf)Cl2", "K3PO4"],
        "conditions": "95°C, dioxane/H2O, 16h",
        "reaction_smiles": "CCCCCc1nc(Cl)c(CO)n1Cc1ccc(C#N)cc1.Brc2ccccc2B(OH)2>>CCCCCc1nc(Cl)c(CO)n1Cc1ccc(-c2ccccc2C#N)cc1.PdBr2"
      },
      {
        "reaction_class": "[3+2] Cycloaddition",
        "reagents": ["TMSN3", "Bu2SnO", "MeOH"],
        "conditions": "reflux, 24h",
        "reaction_smiles": "CCCCCc1nc(Cl)c(CO)n1Cc1ccc(-c2ccccc2C#N)cc1.TMSN3>>CCCCCc1nc(Cl)c(CO)n1Cc1ccc(-c2ccccc2-c2n[nH]nn2)cc1.TMSOH"
      }
    ]
  },
  {
    "route_id": 3,
    "synthesis_style": "convergent",
    "score": 4.6,
    "audit_verdict": "CONDITIONAL",
    "critical_issues": ["Planner output requires validation."],
    "precursors": ["CCCCCc1nc(Cl)c(C=O)n1", "CCc1ccc(-c2ccccc2-c2n[nH]nn2)cc1"],
    "steps": [
      {
        "reaction_class": "Reduction",
        "reagents": ["LiAlH4", "THF"],
        "conditions": "0°C to rt, 3h",
        "reaction_smiles": "CCCCCc1nc(Cl)c(C=O)n1>>CCCCCc1nc(Cl)c(CO)n1"
      },
      {
        "reaction_class": "N-Alkylation",
        "reagents": ["NaH", "DMF"],
        "conditions": "0°C to rt, 12h",
        "reaction_smiles": "CCCCCc1nc(Cl)c(CO)n1.CCc1ccc(-c2ccccc2-c2n[nH]nn2)cc1>>CCCCCc1nc(Cl)c(CO)n1CCc1ccc(-c2ccccc2-c2n[nH]nn2)cc1"
      },
      {
        "reaction_class": "Reduction",
        "reagents": ["H2", "Pd/C", "MeOH"],
        "conditions": "rt, 6h",
        "reaction_smiles": "CCCCCc1nc(Cl)c(CO)n1CCc1ccc(-c2ccccc2-c2n[nH]nn2)cc1>>CCCCCc1nc(Cl)c(CO)n1Cc1ccc(-c2ccccc2-c2n[nH]nn2)cc1"
      },
      {
        "reaction_class": "[3+2] Cycloaddition",
        "reagents": ["NaN3", "ZnCl2", "H2O"],
        "conditions": "100°C, 24h",
        "reaction_smiles": "CCc1ccc(-c2ccccc2C#N)cc1.NaN3>>CCc1ccc(-c2ccccc2-c2n[nH]nn2)cc1.Na"
      },
      {
        "reaction_class": "Suzuki Coupling",
        "reagents": ["Pd(PPh3)4", "Na2CO3"],
        "conditions": "85°C, toluene/EtOH/H2O, 18h",
        "reaction_smiles": "Brc1ccc(C#N)cc1.B(OH)2c2ccccc2>>CCc1ccc(-c2ccccc2C#N)cc1.PdBr2"
      }
    ]
  }
]

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--prompt-file", required=True)
    parser.add_argument("--output-file", required=True)
    parser.add_argument("--request-file", required=True)
    parser.add_argument("--stage", required=True)
    parser.add_argument("--round", type=int, required=True)
    args = parser.parse_args()

    # Write the routes to the output file
    dest = Path(args.output_file)
    dest.parent.mkdir(parents=True, exist_ok=True)

    with open(dest, "w", encoding="utf-8") as f:
        json.dump(ROUTES_DATA, f, indent=2, ensure_ascii=False)

    return 0

if __name__ == "__main__":
    import sys
    sys.exit(main())
