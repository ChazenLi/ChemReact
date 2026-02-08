import io
import os
from typing import Any, Dict, List, Tuple
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


def _canonical_smiles(smiles: str) -> str:
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return ""
    for atom in mol.GetAtoms():
        if atom.GetAtomMapNum():
            atom.SetAtomMapNum(0)
    return Chem.MolToSmiles(mol, canonical=True)


def _parse_reaction_smiles(rxn: str) -> Tuple[List[str], List[str]]:
    if not isinstance(rxn, str) or ">>" not in rxn:
        return [], []
    left, right = rxn.split(">>", 1)
    reactants = [x.strip() for x in left.split(".") if x.strip()]
    products = [x.strip() for x in right.split(".") if x.strip()]
    return reactants, products


def suggest_molecule_image_size(smiles: str) -> Tuple[int, int]:
    """
    Suggest a larger drawing canvas for complex molecules so full structures stay visible.
    """
    mol = Chem.MolFromSmiles(smiles) if isinstance(smiles, str) else None
    if mol is None:
        return (480, 360)
    atom_count = mol.GetNumAtoms()
    width = min(900, max(420, 300 + atom_count * 16))
    height = min(700, max(320, 220 + atom_count * 12))
    return (width, height)


def _render_stage_panel(smiles_list: List[str], title: str, panel_size: Tuple[int, int] = (340, 300)):
    valid_mols = []
    valid_legends = []
    for smi in smiles_list:
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            continue
        AllChem.Compute2DCoords(mol)
        valid_mols.append(mol)
        valid_legends.append(smi)

    if not valid_mols:
        return None

    grid_png = Draw.MolsToGridImage(
        valid_mols,
        molsPerRow=min(3, max(1, len(valid_mols))),
        subImgSize=(220, 180),
        legends=valid_legends,
        returnPNG=True,
    )
    grid = Image.open(io.BytesIO(grid_png)).convert("RGB")
    panel = Image.new("RGB", panel_size, (255, 255, 255))
    draw = ImageDraw.Draw(panel)
    draw.text((10, 8), title, fill=(0, 0, 0))
    grid_x = max(0, (panel_size[0] - grid.width) // 2)
    grid_y = max(28, (panel_size[1] - grid.height) // 2 + 12)
    panel.paste(grid, (grid_x, grid_y))
    draw.rectangle([0, 0, panel_size[0] - 1, panel_size[1] - 1], outline=(180, 180, 180), width=1)
    return panel


def generate_rsps_route_pathway_image(target_smiles: str, route: Dict[str, Any], output_path: str):
    """
    Generate a precursor->target pathway based on step-level RS>>PS reactions.
    Stages: first-step RS (precursors) -> each-step PS (intermediates) -> target.
    """
    if not PIL_AVAILABLE:
        print("PIL not available, skipping RS>>PS pathway generation.")
        return False

    try:
        steps = route.get("steps", []) if isinstance(route, dict) else []
        if not isinstance(steps, list) or not steps:
            return False

        stages: List[Tuple[str, List[str]]] = []
        for idx, step in enumerate(steps, start=1):
            if not isinstance(step, dict):
                continue
            rxn = step.get("reaction_smiles") or step.get("rxn_smiles") or step.get("smirks")
            reactants, products = _parse_reaction_smiles(rxn)
            if idx == 1 and reactants:
                stages.append(("Precursors (RS)", reactants))
            if products:
                stages.append((f"Step {idx} Product (PS)", products))

        if not stages:
            return False

        target_can = _canonical_smiles(target_smiles)
        last_stage_smiles = stages[-1][1] if stages else []
        last_stage_can = {_canonical_smiles(s) for s in last_stage_smiles if _canonical_smiles(s)}
        if target_can and target_can not in last_stage_can:
            stages.append(("Target", [target_smiles]))
        elif target_can and stages[-1][0] != "Target":
            stages[-1] = ("Target", stages[-1][1])

        panels = []
        for title, smiles_list in stages:
            panel = _render_stage_panel(smiles_list, title)
            if panel is not None:
                panels.append(panel)
        if not panels:
            return False

        arrow_w = 70
        margin = 20
        total_w = sum(p.width for p in panels) + arrow_w * (len(panels) - 1) + margin * 2
        total_h = max(p.height for p in panels) + margin * 2
        canvas = Image.new("RGB", (total_w, total_h), (255, 255, 255))
        draw = ImageDraw.Draw(canvas)

        x = margin
        y_base = margin
        for i, panel in enumerate(panels):
            canvas.paste(panel, (x, y_base))
            if i < len(panels) - 1:
                y_mid = y_base + panel.height // 2
                arrow_start = (x + panel.width + 8, y_mid)
                arrow_end = (x + panel.width + arrow_w - 8, y_mid)
                draw.line([arrow_start, arrow_end], fill=(0, 0, 0), width=3)
                draw.polygon(
                    [
                        (arrow_end[0], arrow_end[1]),
                        (arrow_end[0] - 10, arrow_end[1] - 6),
                        (arrow_end[0] - 10, arrow_end[1] + 6),
                    ],
                    fill=(0, 0, 0),
                )
            x += panel.width + arrow_w

        ensure_dir(output_path)
        canvas.save(output_path)
        return True
    except Exception as e:
        print(f"Error generating RS>>PS pathway image: {e}")
        return False

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


def generate_precursor_target_orthogonal_tree_image(
    target_smiles: str,
    precursor_smiles: List[str],
    output_path: str,
    precursor_image_paths: List[str] = None,
):
    """
    Render full RS precursors as individual panels and connect each precursor to target
    using orthogonal (right-angle) guide lines.
    """
    if not PIL_AVAILABLE:
        print("PIL not available, skipping orthogonal precursor tree generation.")
        return False

    try:
        valid_precursors = [s for s in precursor_smiles if isinstance(s, str) and s.strip()]
        if not valid_precursors:
            return False

        target_mol = Chem.MolFromSmiles(target_smiles)
        if target_mol is None:
            return False
        AllChem.Compute2DCoords(target_mol)
        target_size = suggest_molecule_image_size(target_smiles)
        target_img = _mol_to_pil(target_mol, target_size).convert("RGB")

        precursor_imgs: List[Image.Image] = []
        for idx, smi in enumerate(valid_precursors):
            panel = None
            if isinstance(precursor_image_paths, list) and idx < len(precursor_image_paths):
                p = precursor_image_paths[idx]
                if isinstance(p, str) and os.path.exists(p):
                    panel = Image.open(p).convert("RGB")
            if panel is None:
                p_mol = Chem.MolFromSmiles(smi)
                if p_mol is not None:
                    AllChem.Compute2DCoords(p_mol)
                    p_size = suggest_molecule_image_size(smi)
                    panel = _mol_to_pil(p_mol, p_size).convert("RGB")
            if panel is not None:
                precursor_imgs.append(panel)

        if not precursor_imgs:
            return False

        spacing = 30
        margin = 30
        bottom_max_h = max(img.height for img in precursor_imgs)
        total_prec_w = sum(img.width for img in precursor_imgs) + spacing * (len(precursor_imgs) - 1)
        top_w = target_img.width
        canvas_w = max(top_w + 2 * margin, total_prec_w + 2 * margin)
        bus_gap = 80
        section_gap = 90
        canvas_h = margin + target_img.height + section_gap + bus_gap + bottom_max_h + margin

        canvas = Image.new("RGB", (canvas_w, canvas_h), (255, 255, 255))
        draw = ImageDraw.Draw(canvas)

        target_x = (canvas_w - target_img.width) // 2
        target_y = margin
        canvas.paste(target_img, (target_x, target_y))

        prec_start_x = (canvas_w - total_prec_w) // 2
        prec_y = canvas_h - margin - bottom_max_h
        prec_centers = []
        x = prec_start_x
        for panel in precursor_imgs:
            y = prec_y + (bottom_max_h - panel.height)
            canvas.paste(panel, (x, y))
            prec_centers.append((x + panel.width // 2, y))
            x += panel.width + spacing

        target_bottom = (target_x + target_img.width // 2, target_y + target_img.height)
        bus_y = target_bottom[1] + section_gap
        trunk_x = target_bottom[0]
        draw.line([target_bottom, (trunk_x, bus_y)], fill=(0, 0, 0), width=3)

        min_x = min(px for px, _ in prec_centers)
        max_x = max(px for px, _ in prec_centers)
        draw.line([(min_x, bus_y), (max_x, bus_y)], fill=(0, 0, 0), width=3)

        for pc_x, pc_top_y in prec_centers:
            elbow_y = min(pc_top_y - 20, bus_y + 120)
            draw.line([(pc_x, bus_y), (pc_x, elbow_y)], fill=(0, 0, 0), width=3)
            draw.line([(pc_x, elbow_y), (pc_x, pc_top_y)], fill=(0, 0, 0), width=3)

        ensure_dir(output_path)
        canvas.save(output_path)
        return True
    except Exception as e:
        print(f"Error generating orthogonal precursor tree image: {e}")
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
