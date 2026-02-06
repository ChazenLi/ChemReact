
import sys
import json
import argparse
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdForceFieldHelpers

def generate_conformer(smiles, force_field="MMFF", constrained_smarts=None):
    result = {
        "success": False,
        "input_smiles": smiles,
        "conformer_data": {}
    }
    
    try:
        # 1. Prepare Molecule (Add Hs is critical for 3D)
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            result["error"] = "Invalid SMILES"
            return result
        
        mol_h = Chem.AddHs(mol)
        
        # 2. Embed 3D Conformer (ETKDG)
        # Using ETKDGv3 if available, else standard ETKDG
        params = AllChem.ETKDGv3() if hasattr(AllChem, "ETKDGv3") else AllChem.ETKDG()
        params.randomSeed = 0xf00d # Fix seed for reproducibility in tests
        
        embed_code = AllChem.EmbedMolecule(mol_h, params)
        if embed_code == -1:
            result["error"] = "Conformer embedding failed (likely steric clash or complex geometry)."
            return result

        # 3. Setup Force Field
        if force_field == "MMFF":
            # Check availability
            if not AllChem.MMFFHasAllMoleculeParams(mol_h):
                result["warning"] = "MMFF parameters missing for some atoms, falling back to UFF."
                force_field = "UFF"
        
        if force_field == "MMFF":
            props = AllChem.MMFFGetMoleculeProperties(mol_h)
            ff = AllChem.MMFFGetMoleculeForceField(mol_h, props)
        else:
            ff = AllChem.UFFGetMoleculeForceField(mol_h)
            
        if not ff:
            result["error"] = f"Could not create {force_field} force field."
            return result

        # Initial Energy
        initial_energy = ff.CalcEnergy()

        # 4. Apply Constraints (if requested)
        constrained_atoms = []
        if constrained_smarts:
            core = Chem.MolFromSmarts(constrained_smarts)
            if core:
                matches = mol_h.GetSubstructMatches(core)
                if matches:
                    # Use the first match
                    constrained_atoms = list(matches[0])
                    for idx in constrained_atoms:
                        ff.AddFixedPoint(idx)
        
        # 5. Optimize
        ff.Minimize()
        
        final_energy = ff.CalcEnergy()
        
        # 6. Output
        # Remove Hs for display/saving usually, but for 3D PDB usually keep them. 
        # We will return the PDB block of the explicit H molecule.
        pdb_block = Chem.MolToPDBBlock(mol_h)
        
        result["success"] = True
        result["conformer_data"] = {
            "force_field": force_field,
            "initial_energy": initial_energy,
            "final_energy": final_energy,
            "energy_delta": final_energy - initial_energy,
            "constrained_atoms_indices": constrained_atoms,
            "pdb_block": pdb_block
        }

    except Exception as e:
        result["error"] = str(e)
        import traceback
        result["traceback"] = traceback.format_exc()

    return result

def main():
    parser = argparse.ArgumentParser(description="Generate and Optimize 3D Conformer.")
    parser.add_argument("smiles", help="Input SMILES")
    parser.add_argument("--ff", choices=["UFF", "MMFF"], default="MMFF", help="Force Field choice")
    parser.add_argument("--constraint", help="SMARTS pattern to fix during optimization", default=None)
    
    args = parser.parse_args()
    
    output = generate_conformer(args.smiles, args.ff, args.constraint)
    print(json.dumps(output, indent=2))

if __name__ == "__main__":
    main()
