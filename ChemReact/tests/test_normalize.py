"""
Tests for SMILES Normalization Module
"""

import unittest
import sys
from pathlib import Path

# Add project root to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from internal.normalize.smiles_normalizer import (
    normalize_smiles,
    normalize_reaction_smiles,
    remove_atom_numbers,
    normalize_smiles_list
)


class TestRemoveAtomNumbers(unittest.TestCase):
    """Test atom number removal"""
    
    def test_remove_single_mapping(self):
        result = remove_atom_numbers("[CH3:1]O")
        self.assertNotIn(":1", result)
    
    def test_remove_multiple_mappings(self):
        result = remove_atom_numbers("[CH3:1][OH:2]")
        self.assertNotIn(":1", result)
        self.assertNotIn(":2", result)
    
    def test_no_mapping(self):
        result = remove_atom_numbers("CCO")
        self.assertEqual(result, "CCO")
    
    def test_preserve_charge(self):
        result = remove_atom_numbers("[NH4+:1]")
        self.assertIn("+", result)


class TestNormalizeSmiles(unittest.TestCase):
    """Test SMILES normalization"""
    
    def test_basic_smiles(self):
        result = normalize_smiles("CCO")
        self.assertEqual(result["original"], "CCO")
        self.assertIn("display", result)
        self.assertIn("canonical", result)
    
    def test_mapped_smiles(self):
        result = normalize_smiles("[CH3:1][OH:2]")
        self.assertEqual(result["original"], "[CH3:1][OH:2]")
        # display should not contain mapping
        self.assertNotIn(":1", result["display"])
        self.assertNotIn(":2", result["display"])
    
    def test_empty_smiles(self):
        result = normalize_smiles("")
        self.assertEqual(result["original"], "")
    
    def test_invalid_smiles(self):
        result = normalize_smiles("invalid_smiles_xxx")
        # Should return original for invalid
        self.assertEqual(result["original"], "invalid_smiles_xxx")
    
    def test_aspirin(self):
        aspirin = "CC(=O)Oc1ccccc1C(=O)O"
        result = normalize_smiles(aspirin)
        self.assertEqual(result["original"], aspirin)
        # Display should be valid for visualization
        self.assertNotIn(":", result["display"])


class TestNormalizeReactionSmiles(unittest.TestCase):
    """Test reaction SMILES normalization"""
    
    def test_basic_reaction(self):
        rxn = "CCO.CC(=O)O>>CC(=O)OCC"
        result = normalize_reaction_smiles(rxn)
        self.assertEqual(result["original"], rxn)
        self.assertIn(">>", result["display"])
    
    def test_mapped_reaction(self):
        rxn = "[CH3:1]Br.[OH:2]>>[CH3:1][OH:2]"
        result = normalize_reaction_smiles(rxn)
        # Display should not contain mappings
        self.assertNotIn(":1", result["display"])
        self.assertNotIn(":2", result["display"])
        # But should still have reaction arrow
        self.assertIn(">>", result["display"])
    
    def test_no_arrow(self):
        result = normalize_reaction_smiles("CCO.CC")
        # Should return original if no >> separator
        self.assertEqual(result["original"], "CCO.CC")
    
    def test_empty_reaction(self):
        result = normalize_reaction_smiles("")
        self.assertEqual(result["original"], "")


class TestNormalizeSmilesList(unittest.TestCase):
    """Test batch normalization"""
    
    def test_batch(self):
        smiles_list = ["CCO", "[CH3:1]O", "c1ccccc1"]
        results = normalize_smiles_list(smiles_list)
        self.assertEqual(len(results), 3)
        for r in results:
            self.assertIn("original", r)
            self.assertIn("display", r)


class TestVisualizationCompatibility(unittest.TestCase):
    """Test that display format works with RDKit visualization"""
    
    def test_display_parseable(self):
        """Ensure display SMILES can be parsed by RDKit"""
        from rdkit import Chem
        
        test_cases = [
            "[CH3:1][OH:2]",
            "[C:1](=[O:2])([O:3])[c:4]1[cH:5][cH:6][cH:7][cH:8][c:9]1[C:10](=[O:11])[OH:12]",
            "[CH3:1][Br:2].[Na+:3].[OH-:4]>>[CH3:1][OH:5].[Na+:3].[Br-:6]"
        ]
        
        for smiles in test_cases:
            if ">>" in smiles:
                result = normalize_reaction_smiles(smiles)
            else:
                result = normalize_smiles(smiles)
            
            display = result["display"]
            
            # For reactions, split and check each part
            if ">>" in display:
                parts = display.split(">>")
                for part in parts:
                    for frag in part.split("."):
                        if frag.strip():
                            mol = Chem.MolFromSmiles(frag.strip())
                            self.assertIsNotNone(mol, f"Failed to parse: {frag}")
            else:
                mol = Chem.MolFromSmiles(display)
                self.assertIsNotNone(mol, f"Failed to parse: {display}")


if __name__ == "__main__":
    unittest.main()
