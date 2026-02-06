import os
from typing import List, Tuple
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import rdMolDraw2D
try:
    from PIL import Image, ImageDraw
    PIL_AVAILABLE = True
except ImportError:
    PIL_AVAILABLE = False

def ensure_dir(file_path):
    directory = os.path.dirname(file_path)
    if directory and not os.path.exists(directory):
        os.makedirs(directory)

def _mol_to_pil(mol, size):
    """Helper to convert RDKit mol to PIL Image"""
    d = rdMolDraw2D.MolDraw2DCairo(size[0], size[1])
    d.drawOptions().addStereoAnnotation = True
    d.DrawMolecule(mol)
    d.FinishDrawing()
    png_data = d.GetDrawingText()
    import io
    return Image.open(io.BytesIO(png_data))

def generate_molecule_image(smiles: str, output_path: str, legend: str = "", highlight_atoms: List[int] = None, size: Tuple[int, int] = (400, 400)):
    """
    Generates a 2D image of a single molecule.
    """
    try:
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            raise ValueError(f"Invalid SMILES: {smiles}")
        
        AllChem.Compute2DCoords(mol)
        
        d = rdMolDraw2D.MolDraw2DCairo(size[0], size[1])
        d.drawOptions().addStereoAnnotation = True
        
        # Prepare highlights
        if highlight_atoms is None:
            highlight_atoms = []
            
        d.DrawMolecule(mol, highlightAtoms=highlight_atoms, legend=legend)
        d.FinishDrawing()
        
        ensure_dir(output_path)
        d.WriteDrawingText(output_path)
        return True
    except Exception as e:
        print(f"Error generating molecule image: {e}")
        return False

def generate_reaction_image(rxn_smiles: str, output_path: str, size: Tuple[int, int] = (800, 300)):
    """
    Generates an image for a reaction SMILES (Reactants>>Products).
    """
    try:
        rxn = AllChem.ReactionFromSmarts(rxn_smiles, useSmiles=True)
        if not rxn:
             # Fallback: try to manually draw reactants and products side-by-side if strict parsing fails
             # specific parsing/handling might be needed for mapping numbers
             raise ValueError(f"Invalid Reaction SMILES: {rxn_smiles}")

        d = rdMolDraw2D.MolDraw2DCairo(size[0], size[1])
        d.DrawReaction(rxn)
        d.FinishDrawing()
        
        ensure_dir(output_path)
        d.WriteDrawingText(output_path)
        return True
    except Exception as e:
        print(f"Error generating reaction image: {e}")
        return False

def generate_route_grid(smiles_list: List[str], legends: List[str], output_path: str, mols_per_row: int = 3, sub_img_size: Tuple[int, int] = (300, 300)):
    """
    Generates a grid image for a list of molecules (e.g. precursors).
    """
    try:
        mols = [Chem.MolFromSmiles(s) for s in smiles_list]
        if not legends:
            legends = [""] * len(smiles_list)
        elif len(legends) < len(smiles_list):
            legends = legends + [""] * (len(smiles_list) - len(legends))

        # Filter None
        valid_mols = []
        valid_legends = []
        for m, l in zip(mols, legends):
            if m:
                AllChem.Compute2DCoords(m)
                valid_mols.append(m)
                valid_legends.append(l)

        if not valid_mols:
            raise ValueError("No valid molecules available for grid rendering.")
        
        img = Draw.MolsToGridImage(valid_mols, molsPerRow=mols_per_row, subImgSize=sub_img_size, legends=valid_legends, returnPNG=True)
        
        ensure_dir(output_path)
        with open(output_path, 'wb') as f:
            f.write(img)
        return True
    except Exception as e:
        print(f"Error generating route grid: {e}")
        return False

def generate_reaction_tree_image(target_smiles: str, precursor_smiles: List[str], output_path: str, canvas_size: Tuple[int, int] = (1000, 800)):
    """
    Generates a Tree View image: Target (Top) -> Arrow -> Precursors (Bottom).
    Requires PIL (Pillow).
    """
    if not PIL_AVAILABLE:
        print("PIL not available, skipping tree generation.")
        return False

    try:
        # 1. Prepare Molecules
        target_mol = Chem.MolFromSmiles(target_smiles)
        precursor_mols = [Chem.MolFromSmiles(s) for s in precursor_smiles if s]
        
        if not target_mol or not precursor_mols:
            print("Invalid SMILES for tree generation.")
            return False
            
        AllChem.Compute2DCoords(target_mol)
        for m in precursor_mols: AllChem.Compute2DCoords(m)

        # 2. Generate Sub-Images
        mol_w, mol_h = 300, 300
        target_img = _mol_to_pil(target_mol, (mol_w, mol_h))
        prec_imgs = [_mol_to_pil(m, (mol_w, mol_h)) for m in precursor_mols]
        
        # 3. Layout Dimensions
        num_precs = len(prec_imgs)
        spacing = 50
        arrow_section_h = 100
        
        total_prec_w = num_precs * mol_w + (num_precs - 1) * spacing
        canvas_w = max(mol_w, total_prec_w) + 100
        canvas_h = mol_h + arrow_section_h + mol_h + 50
        canvas_w = max(canvas_w, canvas_size[0])
        canvas_h = max(canvas_h, canvas_size[1])
        
        # 4. Draw Canvas
        bg_color = (255, 255, 255)
        canvas = Image.new('RGB', (canvas_w, canvas_h), bg_color)
        draw = ImageDraw.Draw(canvas)
        
        # 5. Place Target (Top Center)
        target_x = (canvas_w - mol_w) // 2
        target_y = 20
        canvas.paste(target_img, (target_x, target_y))
        
        # 6. Place Precursors (Bottom Row)
        prec_start_x = (canvas_w - total_prec_w) // 2
        prec_y = target_y + mol_h + arrow_section_h
        
        prec_centers_x = []
        
        for i, p_img in enumerate(prec_imgs):
            x = prec_start_x + i * (mol_w + spacing)
            canvas.paste(p_img, (x, prec_y))
            prec_centers_x.append(x + mol_w // 2)
            
        # 7. Draw Orthogonal Retro Tree Connectors (Target -> Precursor1 + Precursor2 + ...)
        target_bottom_center = (target_x + mol_w // 2, target_y + mol_h)
        trunk_x = target_bottom_center[0]
        # Keep bus well above precursor panels to avoid visual overlap.
        bus_y = target_bottom_center[1] + max(25, arrow_section_h // 3)
        bus_y = min(bus_y, prec_y - 24)

        # Main trunk from target down to horizontal bus.
        draw.line([target_bottom_center, (trunk_x, bus_y)], fill=(0, 0, 0), width=3)

        # Horizontal bus from leftmost to rightmost precursor branch.
        min_x = min(prec_centers_x)
        max_x = max(prec_centers_x)
        draw.line([(min_x, bus_y), (max_x, bus_y)], fill=(0, 0, 0), width=3)

        # Vertical branches from bus to each precursor top center.
        for pc_x in prec_centers_x:
            draw.line([(pc_x, bus_y), (pc_x, prec_y)], fill=(0, 0, 0), width=3)

        # Add '+' markers between precursors to emphasize "pre1 + pre2" composition.
        if len(prec_centers_x) > 1:
            plus_y = prec_y + mol_h // 2
            for i in range(len(prec_centers_x) - 1):
                mid_x = (prec_centers_x[i] + prec_centers_x[i + 1]) // 2
                draw.text((mid_x - 6, plus_y - 12), "+", fill=(0, 0, 0))
            
        # Save
        ensure_dir(output_path)
        canvas.save(output_path)
        return True

    except Exception as e:
        print(f"Error generating tree image: {e}")
        import traceback
        traceback.print_exc()
        return False
