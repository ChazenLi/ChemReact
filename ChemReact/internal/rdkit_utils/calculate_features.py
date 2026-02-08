
import sys
import json
import argparse
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import QED
from rdkit.Chem import AllChem
from rdkit.Chem import GraphDescriptors
from rdkit.Chem.Scaffolds import MurckoScaffold


def _calculate_sa_score(mol) -> float:
    """
    Calculate Synthetic Accessibility Score (SAScore).
    Lower scores = easier to synthesize (1-10 scale).
    Uses fragment-based complexity estimation.
    """
    try:
        # Try to use RDKit's built-in sascorer if available
        from rdkit.Chem import RDConfig
        import os
        sascore_path = os.path.join(RDConfig.RDContribDir, 'SA_Score', 'sascorer.py')
        if os.path.exists(sascore_path):
            import importlib.util
            spec = importlib.util.spec_from_file_location("sascorer", sascore_path)
            sascorer = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(sascorer)
            return round(sascorer.calculateScore(mol), 3)
    except Exception:
        pass
    
    # Fallback: simplified SA estimation based on molecular features
    # Based on: Ertl & Schuffenhauer, J. Cheminform. 2009
    num_rings = rdMolDescriptors.CalcNumRings(mol)
    num_stereo = len(Chem.FindMolChiralCenters(mol, includeUnassigned=True))
    num_rotatable = Descriptors.NumRotatableBonds(mol)
    num_hetero = rdMolDescriptors.CalcNumHeteroatoms(mol)
    mol_wt = Descriptors.MolWt(mol)
    
    # Complexity factors (higher = harder to synthesize)
    ring_factor = 0.3 * num_rings
    stereo_factor = 0.5 * num_stereo
    size_factor = 0.01 * max(0, mol_wt - 150)
    hetero_factor = 0.1 * num_hetero
    
    # Base score + complexity
    sa_score = 2.0 + ring_factor + stereo_factor + size_factor + hetero_factor
    return round(min(10.0, max(1.0, sa_score)), 3)


def _get_scaffold_info(mol) -> dict:
    """Extract Murcko scaffold information."""
    try:
        core = MurckoScaffold.GetScaffoldForMol(mol)
        core_smiles = Chem.MolToSmiles(core)
        framework = MurckoScaffold.MakeScaffoldGeneric(core)
        framework_smiles = Chem.MolToSmiles(framework)
        return {
            "core_scaffold": core_smiles,
            "generic_framework": framework_smiles,
            "scaffold_rings": rdMolDescriptors.CalcNumRings(core),
        }
    except Exception:
        return {"core_scaffold": None, "generic_framework": None, "scaffold_rings": 0}


def calculate_features(smiles: str):
    result = {
        "input_smiles": smiles,
        "success": False,
        "error": None,
        "synthesis_critical": {},
        "metadata": {
            "synthesis_critical": "Key properties for synthesis and retrosynthesis.",
            "supplementary": "Additional descriptors (LogP, QED, Atomic Features) less critical for initial planning."
        },
        "supplementary": {}
    }

    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            result["error"] = "Invalid SMILES string"
            return result

        # Sanitize for property calculation
        Chem.SanitizeMol(mol)

        # --- 1. Synthesis Critical Features ---
        sa_score = _calculate_sa_score(mol)
        scaffold_info = _get_scaffold_info(mol)
        
        result["synthesis_critical"] = {
            "MolWt": round(Descriptors.MolWt(mol), 3),
            "TPSA": round(Descriptors.TPSA(mol), 3),
            "H_Bond_Donors": Descriptors.NumHDonors(mol),
            "H_Bond_Acceptors": Descriptors.NumHAcceptors(mol),
            "SAScore": sa_score,  # Synthetic Accessibility Score
            "Scaffold": scaffold_info,
            "Stereocenters": len(Chem.FindMolChiralCenters(mol, includeUnassigned=True)),
            "NumRings": rdMolDescriptors.CalcNumRings(mol),
            "NumHeteroatoms": rdMolDescriptors.CalcNumHeteroatoms(mol),
        }

        # --- 2. Supplementary Features ---
        
        result["supplementary"] = {
             "LogP": round(Descriptors.MolLogP(mol), 3),
             "Rotatable_Bonds": Descriptors.NumRotatableBonds(mol),
             "QED": round(QED.qed(mol), 3),
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
