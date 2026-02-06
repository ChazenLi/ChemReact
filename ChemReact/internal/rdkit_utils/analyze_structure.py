
import sys
import json
import argparse
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors

def analyze_smiles(smiles: str, generate_3d: bool = False):
    """
    Analyzes a SMILES string using RDKit and returns a dictionary of results.
    """
    result = {
        "input_smiles": smiles,
        "success": False,
        "error": None,
        "properties": {}
    }

    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            result["error"] = "Invalid SMILES string: Could not parse molecule."
            return result

        # 1. Sanitization (fixes aromaticity, valence, etc.)
        # RDKit does this by default on MolFromSmiles, but explicit is good for catch
        sanitize_flags = Chem.SanitizeFlags.SANITIZE_ALL
        sanitization_result = Chem.SanitizeMol(mol, sanitizeOps=sanitize_flags, catchErrors=True)
        
        if sanitization_result != Chem.SanitizeFlags.SANITIZE_NONE:
             result["warnings"] = f"Sanitization issues found: {sanitization_result}"

        # 2. Stereochemistry
        Chem.AssignStereochemistry(mol, force=True, cleanIt=True)
        chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)

        # 3. Canonicalization
        canonical_smiles = Chem.MolToSmiles(mol, isomericSmiles=True, canonical=True)

        # 4. Basic Properties
        formula = rdMolDescriptors.CalcMolFormula(mol)
        mol_wt = Descriptors.MolWt(mol)
        formal_charge = Chem.GetFormalCharge(mol)
        num_atoms = mol.GetNumAtoms()
        num_bonds = mol.GetNumBonds()
        
        # 5. 3D Conformer (Optional)
        has_3d = False
        if generate_3d:
            try:
                mol_3d = Chem.AddHs(mol)
                AllChem.EmbedMolecule(mol_3d, AllChem.ETKDG())
                has_3d = True
            except Exception as e:
                result["warnings_3d"] = str(e)

        result["success"] = True
        result["properties"] = {
            "canonical_smiles": canonical_smiles,
            "formula": formula,
            "molecular_weight": mol_wt,
            "formal_charge": formal_charge,
            "num_atoms": num_atoms,
            "num_bonds": num_bonds,
            "chiral_centers": chiral_centers,
            "stereochemistry_valid": len(chiral_centers) > 0 or "No stereocenters"
        }
        
    except Exception as e:
        result["error"] = str(e)

    return result

def main():
    parser = argparse.ArgumentParser(description="Analyze a SMILES string using RDKit.")
    parser.add_argument("smiles", help="Input SMILES string")
    parser.add_argument("--3d", dest="generate_3d", action="store_true", help="Attempt to generate 3D conformer")
    args = parser.parse_args()

    analysis = analyze_smiles(args.smiles, generate_3d=args.generate_3d)
    print(json.dumps(analysis, indent=2))

if __name__ == "__main__":
    main()
