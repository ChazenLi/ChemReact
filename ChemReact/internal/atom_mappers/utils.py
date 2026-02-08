"""
Utility Functions for Atom Mapping - LLM Focused

Provides helper functions to extract key mapping information for LLM analysis.
"""

import re
from typing import Dict, List, Set, Tuple, Any
from rdkit import Chem


def _extract_mapped_bonds(side_smiles: str) -> Dict[Tuple[int, int], str]:
    """Return mapped bond pairs and bond type token from one reaction side."""
    bond_map: Dict[Tuple[int, int], str] = {}
    for frag in side_smiles.split("."):
        frag = frag.strip()
        if not frag:
            continue
        mol = Chem.MolFromSmiles(frag)
        if mol is None:
            continue
        for bond in mol.GetBonds():
            a1 = bond.GetBeginAtom().GetAtomMapNum()
            a2 = bond.GetEndAtom().GetAtomMapNum()
            if a1 <= 0 or a2 <= 0:
                continue
            pair = (a1, a2) if a1 < a2 else (a2, a1)
            bond_map[pair] = str(bond.GetBondType())
    return bond_map


def extract_atom_correspondence(mapped_smiles: str) -> Dict[int, Dict[str, str]]:
    """
    Extract atom correspondence between reactants and products.

    Args:
        mapped_smiles: Mapped reaction SMILES (e.g., "[CH3:1][OH:2]>>[CH3:1][O:2]")

    Returns:
        {atom_number: {"reactant": atom_type, "product": atom_type}}
        Example: {1: {"reactant": "CH3", "product": "CH3"}, 2: {"reactant": "OH", "product": "O"}}

    LLM Usage:
        This tells the LLM exactly which atoms map from reactants to products,
        enabling it to track bond formation/breaking at the atomic level.
    """
    if ">>" not in mapped_smiles:
        return {}

    pattern = re.compile(r'\[([A-Za-z][A-Za-z0-9+@-]*):(\d+)\]')
    reactants, products = mapped_smiles.split(">>", 1)

    # Parse reactants
    reactant_atoms = {}
    for match in pattern.finditer(reactants):
        atom_type = match.group(1)
        atom_num = int(match.group(2))
        if atom_num not in reactant_atoms:
            reactant_atoms[atom_num] = atom_type

    # Parse products
    product_atoms = {}
    for match in pattern.finditer(products):
        atom_type = match.group(1)
        atom_num = int(match.group(2))
        if atom_num not in product_atoms:
            product_atoms[atom_num] = atom_type

    # Build correspondence
    correspondence = {}
    for atom_num, r_type in reactant_atoms.items():
        correspondence[atom_num] = {
            "reactant": r_type,
            "product": product_atoms.get(atom_num, "?")
        }

    return correspondence


def get_reaction_center(mapped_smiles: str) -> Set[int]:
    """
    Identify atoms at the reaction center (atoms that changed).

    Args:
        mapped_smiles: Mapped reaction SMILES

    Returns:
        Set of atom numbers that are part of the reaction center

    LLM Usage:
        The LLM can use this to identify which atoms are involved in
        bond formation/breaking, helping it understand the reaction mechanism.
    """
    correspondence = extract_atom_correspondence(mapped_smiles)
    reaction_center = set()

    for atom_num, types in correspondence.items():
        if types["reactant"] != types["product"]:
            reaction_center.add(atom_num)

    return reaction_center


def get_atom_changes(mapped_smiles: str) -> List[Dict[str, Any]]:
    """
    Identify atoms that changed between reactants and products.

    Args:
        mapped_smiles: Mapped reaction SMILES

    Returns:
        List of changes:
        [
            {"atom_number": 1, "reactant_type": "C", "product_type": "O", "change_type": "modified"},
            {"atom_number": 2, "reactant_type": "Br", "product_type": None, "change_type": "removed"},
            {"atom_number": 3, "reactant_type": None, "product_type": "O", "change_type": "added"}
        ]

    LLM Usage:
        Provides detailed change information for each atom, helping the LLM
        provide step-by-step mechanistic explanations.
    """
    changes = []
    correspondence = extract_atom_correspondence(mapped_smiles)

    # Track changes
    for atom_num, types in correspondence.items():
        r_type = types["reactant"]
        p_type = types["product"]

        if p_type == "?":
            changes.append({
                "atom_number": atom_num,
                "reactant_type": r_type,
                "product_type": None,
                "change_type": "removed"
            })
        elif r_type != p_type:
            changes.append({
                "atom_number": atom_num,
                "reactant_type": r_type,
                "product_type": p_type,
                "change_type": "modified"
            })

    # Check for added atoms (in products but not reactants)
    if ">>" in mapped_smiles:
        _, products = mapped_smiles.split(">>", 1)
        pattern = re.compile(r'\[([A-Za-z][A-Za-z0-9+@-]*):(\d+)\]')

        product_nums = set()
        product_types = {}
        for match in pattern.finditer(products):
            atom_type = match.group(1)
            atom_num = int(match.group(2))
            product_nums.add(atom_num)
            product_types[atom_num] = atom_type

        for atom_num in product_nums:
            if atom_num not in correspondence:
                changes.append({
                    "atom_number": atom_num,
                    "reactant_type": None,
                    "product_type": product_types[atom_num],
                    "change_type": "added"
                })

    return changes


def get_bond_changes(mapped_smiles: str) -> Dict[str, List[Dict[str, Any]]]:
    """
    Identify bond changes in the reaction.

    Args:
        mapped_smiles: Mapped reaction SMILES

    Returns:
        {
            "formed": [{"atom1": 1, "atom2": 2}, ...],  # new bonds in products
            "broken": [{"atom1": 1, "atom2": 2}, ...]   # bonds broken from reactants
        }

    LLM Usage:
        Helps LLM understand which bonds are formed and broken,
        crucial for mechanistic analysis.
    """
    if ">>" not in mapped_smiles:
        return {"formed": [], "broken": []}

    reactants, products = mapped_smiles.split(">>", 1)
    reactant_bonds = _extract_mapped_bonds(reactants)
    product_bonds = _extract_mapped_bonds(products)

    formed: List[Dict[str, Any]] = []
    broken: List[Dict[str, Any]] = []

    reactant_pairs = set(reactant_bonds.keys())
    product_pairs = set(product_bonds.keys())

    for atom1, atom2 in sorted(product_pairs - reactant_pairs):
        formed.append({
            "atom1": atom1,
            "atom2": atom2,
            "bond_type": product_bonds[(atom1, atom2)],
        })

    for atom1, atom2 in sorted(reactant_pairs - product_pairs):
        broken.append({
            "atom1": atom1,
            "atom2": atom2,
            "bond_type": reactant_bonds[(atom1, atom2)],
        })

    # Bond order changes are represented as one broken + one formed on same pair.
    for atom1, atom2 in sorted(reactant_pairs & product_pairs):
        reactant_type = reactant_bonds[(atom1, atom2)]
        product_type = product_bonds[(atom1, atom2)]
        if reactant_type != product_type:
            broken.append({
                "atom1": atom1,
                "atom2": atom2,
                "bond_type": reactant_type,
            })
            formed.append({
                "atom1": atom1,
                "atom2": atom2,
                "bond_type": product_type,
            })

    return {"formed": formed, "broken": broken}


def summarize_for_llm(mapped_smiles: str, confidence: float = 0.0) -> str:
    """
    Generate a text summary of the atom mapping for LLM consumption.

    Args:
        mapped_smiles: Mapped reaction SMILES
        confidence: Mapping confidence score

    Returns:
        Human-readable summary with key insights for LLM

    Example Output:
        "Atom mapping confidence: 95.6%
         Reaction center involves 3 atoms: [1, 2, 3]
         Key changes:
         - Atom 1: C -> O (oxidation)
         - Atom 2: Br -> H (substitution)
         - Bond formed: C1-O2
         - Bond broken: C2-Br3"
    """
    lines = []
    lines.append(f"Atom mapping confidence: {confidence:.1%}")

    # Get reaction center
    reaction_center = get_reaction_center(mapped_smiles)
    if reaction_center:
        lines.append(f"Reaction center atoms: {sorted(reaction_center)}")

    # Get changes
    changes = get_atom_changes(mapped_smiles)
    if changes:
        lines.append("Atom changes:")
        for change in changes:
            atom = change["atom_number"]
            r_type = change.get("reactant_type", "?")
            p_type = change.get("product_type", "?")
            c_type = change["change_type"]
            lines.append(f"  Atom {atom}: {r_type} -> {p_type} ({c_type})")

    bond_changes = get_bond_changes(mapped_smiles)
    if bond_changes["formed"]:
        lines.append("Bonds formed:")
        for bond in bond_changes["formed"]:
            lines.append(
                f"  {bond['atom1']}-{bond['atom2']} ({bond.get('bond_type', '?')})"
            )
    if bond_changes["broken"]:
        lines.append("Bonds broken:")
        for bond in bond_changes["broken"]:
            lines.append(
                f"  {bond['atom1']}-{bond['atom2']} ({bond.get('bond_type', '?')})"
            )

    return "\n".join(lines)


def remove_atom_numbers(smiles: str) -> str:
    """
    Remove atom mapping numbers from SMILES.

    Args:
        smiles: SMILES with or without atom numbers

    Returns:
        SMILES without atom numbers

    LLM Usage:
        Useful when the LLM needs the canonical form without mapping annotations.
    """
    return re.sub(r'\[([A-Za-z][A-Za-z0-9+@-]*):(\d+)\]', r'[\1]', smiles)


def validate_mapped_reaction(mapped_smiles: str) -> Dict[str, Any]:
    """
    Validate a mapped reaction and return key metrics.

    Returns:
        {
            "valid": bool,
            "total_atoms": int,
            "mapped_atoms": int,
            "completeness": float,  # fraction of atoms that are mapped
            "warnings": List[str]
        }
    """
    result = {
        "valid": True,
        "total_atoms": 0,
        "mapped_atoms": 0,
        "completeness": 0.0,
        "warnings": []
    }

    if ">>" not in mapped_smiles:
        result["valid"] = False
        result["warnings"].append("Missing '>>' separator")
        return result

    pattern = re.compile(r'\[([A-Za-z][A-Za-z0-9+@-]*):(\d+)\]')
    all_atoms = re.compile(r'\[([A-Za-z][A-Za-z0-9+@-]*)\]')

    reactants, products = mapped_smiles.split(">>", 1)

    # Count atoms
    total_reactants = len(all_atoms.findall(reactants))
    total_products = len(all_atoms.findall(products))
    result["total_atoms"] = total_reactants + total_products

    mapped_reactants = set()
    for match in pattern.finditer(reactants):
        mapped_reactants.add(int(match.group(2)))

    mapped_products = set()
    for match in pattern.finditer(products):
        mapped_products.add(int(match.group(2)))

    result["mapped_atoms"] = len(mapped_reactants | mapped_products)
    result["completeness"] = result["mapped_atoms"] / result["total_atoms"] if result["total_atoms"] > 0 else 0

    # Check for unmapped atoms
    unmapped_reactants = total_reactants - len(mapped_reactants)
    unmapped_products = total_products - len(mapped_products)

    if unmapped_reactants > 0:
        result["warnings"].append(f"{unmapped_reactants} unmapped reactant atom(s)")

    if unmapped_products > 0:
        result["warnings"].append(f"{unmapped_products} unmapped product atom(s)")

    return result
