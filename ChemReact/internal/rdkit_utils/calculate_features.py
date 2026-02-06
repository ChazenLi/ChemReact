
import sys
import json
import argparse
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import QED
from rdkit.Chem import AllChem
from rdkit.Chem import GraphDescriptors

def calculate_features(smiles: str):
    result = {
        "input_smiles": smiles,
        "success": False, # This will be updated to True later if successful
        "error": None,
        "synthesis_critical": {},
        "metadata": {
            "synthesis_critical": "Key properties for synthesis and retrosynthesis.",
            "supplementary": "Additional descriptors (LogP, QED, Atomic Features) less critical for initial planning."
        },
        "supplementary": {} # Initialize supplementary here, it will be populated later
    }

    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            result["error"] = "Invalid SMILES string"
            return result

        # Sanitize for property calculation
        Chem.SanitizeMol(mol)

        # --- 1. Synthesis Critical Features ---
        
        # atomic_features loop removed (Moved to Module 8: atom_supplementary.py)
        
        result["synthesis_critical"] = {
            "MolWt": Descriptors.MolWt(mol),
            "TPSA": Descriptors.TPSA(mol),
            "H_Bond_Donors": Descriptors.NumHDonors(mol),
            "H_Bond_Acceptors": Descriptors.NumHAcceptors(mol)
        }

        # --- 2. Supplementary Features ---
        
        result["supplementary"] = {
             "LogP": Descriptors.MolLogP(mol),
             "Rotatable_Bonds": Descriptors.NumRotatableBonds(mol),
             "QED": QED.qed(mol),
             # "atomic_features": atoms, # Moved to Module 8
             "topological_descriptors": {}
        }

        try:
             result["supplementary"]["topological_descriptors"] = {
                "BertzCT": GraphDescriptors.BertzCT(mol),
                "BalabanJ": GraphDescriptors.BalabanJ(mol),
                "HallKierAlpha": GraphDescriptors.HallKierAlpha(mol),
                "Chi0": GraphDescriptors.Chi0(mol),
                "Number_of_Rings": rdMolDescriptors.CalcNumRings(mol)
             }
        except Exception as e:
             result["supplementary"]["topological_descriptors_error"] = str(e)

        result["success"] = True

    except Exception as e:
        result["error"] = str(e)

    return result

def main():
    parser = argparse.ArgumentParser(description="Calculate chemical features using RDKit.")
    parser.add_argument("smiles", help="Input SMILES string")
    args = parser.parse_args()

    features = calculate_features(args.smiles)
    print(json.dumps(features, indent=2))

if __name__ == "__main__":
    main()
