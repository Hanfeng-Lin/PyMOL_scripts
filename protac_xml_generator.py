import tkinter as tk
from tkinter import filedialog, messagebox
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolTransforms, rdFMCS
from rdkit.Chem import Draw, rdDetermineBonds
from rdkit.Chem.SaltRemover import SaltRemover
from rdkit.Chem.Draw import rdMolDraw2D
from io import BytesIO
from PIL import Image, ImageTk
from lxml import etree
import re, io, os
import numpy as np

filename = ""
xml_file_path = "PROTAC.xml"
protac_xml = None
drawOptions = rdMolDraw2D.MolDrawOptions()
drawOptions.annotationFontScale = 1


class AtomSelectorApp:
    def __init__(self, master):
        self.master = master
        self.master.title("Atom Selector")
        self.master.columnconfigure(0, minsize=300)
        self.master.columnconfigure(1, minsize=500)

        self.file_path_label = tk.Label(master, text="Select linker conformer PDB file:")
        self.file_path_label.grid(row=0, column=0)

        self.file_path_entry = tk.Entry(master)
        self.file_path_entry.grid(row=0, column=1, sticky="we")

        self.browse_button = tk.Button(master, text="Browse", command=self.browse_file)
        self.browse_button.grid(row=0, column=2)

        # Canvas for linker
        self.canvas = tk.Canvas(master, width=500, height=500)
        self.canvas.grid(row=2, rowspan=2, column=1, columnspan=2)
        # Canvas for warhead
        self.warhead_canvas = tk.Canvas(master, width=400, height=400)
        self.warhead_canvas.grid(row=2, column=5)
        # Canvas for E3 ligand image
        self.e3ligand_canvas = tk.Canvas(master, width=400, height=400)
        self.e3ligand_canvas.grid(row=2, column=6)

        # Text widget to display ligand atoms
        self.warhead_ligand_output_label = tk.Label(master, text="Warhead ligand")
        self.warhead_ligand_output_label.grid(row=3, column=5)
        self.warhead_ligand_output_text = tk.Text(master, width=40, height=10)
        self.warhead_ligand_output_text.grid(row=4, column=5, rowspan=5, padx=5, pady=5, sticky="we")

        self.e3_ligand_output_label = tk.Label(master, text="E3 ligand")
        self.e3_ligand_output_label.grid(row=3, column=6)
        self.e3_ligand_output_text = tk.Text(master, width=40, height=10)
        self.e3_ligand_output_text.grid(row=4, column=6, rowspan=5, padx=5, pady=5, sticky="we")

        # Align button
        self.rotate_button = tk.Button(master, text="Try alignment", command=self.align_molecules)
        self.rotate_button.grid(row=3, column=0, sticky="e", padx=5)

        # Rotate button
        self.rotate_button = tk.Button(master, text="Rotate 90", command=self.rotate_image)
        self.rotate_button.grid(row=3, column=1, sticky="e", padx=5)

        # Flip button
        self.flip_button = tk.Button(master, text="Flip", command=self.flip_image)
        self.flip_button.grid(row=3, column=2, sticky="w", padx=5)

        # Label and entry for warhead atoms
        self.warhead_atoms_label = tk.Label(master, text="Warhead Atoms for alignment:")
        self.warhead_atoms_label.grid(row=4, column=0, sticky="e")
        self.warhead_atoms_entry = tk.Entry(master)
        self.warhead_atoms_entry.grid(row=4, column=1, padx=5, pady=5, sticky="we")

        # Label and entry for bridge point warhead
        self.bridge_point_warhead_label = tk.Label(master, text="Bridge Point for Warhead:")
        self.bridge_point_warhead_label.grid(row=5, column=0, sticky="e")
        self.bridge_point_warhead_entry = tk.Entry(master)
        self.bridge_point_warhead_entry.grid(row=5, column=1, padx=5, pady=5, sticky="we")

        # Label and entry for E3 atoms
        self.E3_atoms_label = tk.Label(master, text="E3 Atoms for alignment:")
        self.E3_atoms_label.grid(row=6, column=0, sticky="e")
        self.E3_atoms_entry = tk.Entry(master)
        self.E3_atoms_entry.grid(row=6, column=1, padx=5, pady=5, sticky="we")

        # Label and entry for bridge point E3 ligand
        self.bridge_point_E3ligand_label = tk.Label(master, text="Bridge Point for E3 Ligand:")
        self.bridge_point_E3ligand_label.grid(row=7, column=0, sticky="e")
        self.bridge_point_E3ligand_entry = tk.Entry(master)
        self.bridge_point_E3ligand_entry.grid(row=7, column=1, padx=5, pady=5, sticky="we")

        # Label and entry for useless atoms
        self.useless_atoms_label = tk.Label(master, text="Useless Atoms to be deleted when bridging:")
        self.useless_atoms_label.grid(row=8, column=0, sticky="e")
        self.useless_atoms_entry = tk.Entry(master)
        self.useless_atoms_entry.grid(row=8, column=1, padx=5, pady=5, sticky="we")

        # Button to generate XML
        self.generate_button = tk.Button(master, text="Generate XML", command=self.generate_xml)
        self.generate_button.grid(row=9, column=2)
        # Button to refresh image
        self.refresh_button = tk.Button(master, text="Refresh Image", command=self.refresh_image)
        self.refresh_button.grid(row=9, column=1, pady=10)
        # Button to clear input box
        self.clear_button = tk.Button(master, text="Clear", command=self.clear_input)
        self.clear_button.grid(row=9, column=0, sticky="e")

        # Text widget to display XML output
        self.xml_output_label = tk.Label(master, text="XML Output:")
        self.xml_output_label.grid(row=10, column=0, sticky="e")
        self.xml_output_text = tk.Text(master, width=60, height=10)
        self.xml_output_text.grid(row=10, column=1, columnspan=2, padx=5, pady=5, sticky="we")

        # Create dropdown menu for warhead ligands
        self.warhead_var = tk.StringVar()
        self.warhead_var.set("Select Warhead")
        self.warhead_var.trace("w", self.handle_selection)
        self.warhead_dropdown = tk.OptionMenu(master, self.warhead_var, "Select Warhead")
        self.warhead_dropdown.grid(row=0, column=5, padx=5, pady=5)

        # Create dropdown menu for E3 ligands
        self.e3ligand_var = tk.StringVar()
        self.e3ligand_var.set("Select E3 Ligand")
        self.e3ligand_var.trace("w", self.handle_selection)
        self.e3ligand_dropdown = tk.OptionMenu(master, self.e3ligand_var, "Select E3 Ligand")
        self.e3ligand_dropdown.grid(row=0, column=6, padx=5, pady=5)

        # Populate dropdown menus
        self.populate_dropdown()

        self.mol = None
        self.warhead_mol = None
        self.warhead_bridge_point_id, self.warhead_align_atom_ids, self.warhead_overlapping_atom_ids = None, None, None
        self.e3_mol = None
        self.e3_bridge_point_id, self.e3_align_atom_ids, self.e3_overlapping_atom_ids = None, None, None
        self.linker_highlight_map = {}
        self.linker_highlight_radii = {}

    def browse_file(self):
        file_path = filedialog.askopenfilename(filetypes=[("PDB Files", "*.pdb")])
        global filename
        filename = file_path.split('/')[-1][:-4]
        if file_path:
            self.file_path_entry.delete(0, tk.END)
            self.file_path_entry.insert(0, file_path)
            self.draw_molecule(file_path)

    def draw_molecule(self, pdb_file_path):
        """
        mol = Chem.MolFromPDBFile(pdb_file_path, sanitize=True, removeHs=False)
        rdDetermineBonds.DetermineBonds(mol, charge=0)
        """
        mol = Chem.MolFromPDBFile(pdb_file_path, removeHs=False)
        rdDetermineBonds.DetermineBonds(mol, charge=0)
        Chem.SanitizeMol(mol)
        mol = Chem.RemoveHs(mol)

        print(Chem.MolToSmiles(mol))
        mol.SetProp("_displayOptions", '{"kekulize": true, "addStereoAnnotation": true}')
        for atom in mol.GetAtoms():
            atom.SetProp('atomNote', str(atom.GetIdx() + 1))
        AllChem.Compute2DCoords(mol)
        if mol is not None:
            self.mol = mol
            img = self.get_molecule_image(mol)
            if img:
                print("Image size:", img.width(), "x", img.height())
                self.display_image(img, self.canvas)
            else:
                print("Failed to create image")
        else:
            messagebox.showerror("Error", "Failed to load molecule from PDB file.")

    def get_molecule_image(self, mol, size=500, highlight_atom_map=None, highlight_radii=None):
        drawer = rdMolDraw2D.MolDraw2DCairo(int(size), int(size))  # Adjust the size as needed
        drawer.SetDrawOptions(drawOptions)

        # Draw the molecule
        if highlight_atom_map:
            drawer.DrawMoleculeWithHighlights(mol, "", highlight_atom_map, {}, highlight_radii, {})
        else:
            drawer.DrawMolecule(mol)
        drawer.FinishDrawing()
        # Convert the drawing to a PIL image
        png_img = drawer.GetDrawingText()
        # Open the PNG image using PIL
        pil_img = Image.open(io.BytesIO(png_img))
        # Convert PIL image to Tkinter PhotoImage
        photo_image = ImageTk.PhotoImage(pil_img)
        return photo_image

    def display_image(self, img, canvas):
        self.canvas.delete("all")
        # Draw the image
        self.image_label = tk.Label(canvas, image=img)
        self.image_label.image = img  # Keep a reference to avoid garbage collection
        self.image_label.place(x=0, y=0, anchor=tk.NW)

    def refresh_image(self):
        # Get values from entry widgets
        pdb_file_path = self.file_path_entry.get()
        warhead_atoms = self.warhead_atoms_entry.get().split()
        bridge_point_warhead = self.bridge_point_warhead_entry.get().split()
        E3_atoms = self.E3_atoms_entry.get().split()
        bridge_point_E3ligand = self.bridge_point_E3ligand_entry.get().split()
        useless_atoms = self.useless_atoms_entry.get().split()

        mol = self.mol

        mol.SetProp("_displayOptions", '{"kekulize": true, "addStereoAnnotation": true}')
        for atom in mol.GetAtoms():
            atom.SetProp('atomNote', str(atom.GetIdx() + 1))
        if mol is None:
            messagebox.showerror("Error", "Failed to read molecule from PDB file.")
            return

        # Highlight atoms
        highlight_atom_map = {}
        highlight_radii = {}

        def add_atoms_to_highlight_map(atom_indices, rgb_tuple):
            for atom in atom_indices:
                atom_index = int(''.join(re.findall(r'\d+', atom))) - 1
                if 0 <= atom_index < mol.GetNumAtoms():
                    if atom_index in highlight_atom_map:
                        highlight_atom_map[atom_index].append(rgb_tuple)
                    else:
                        highlight_atom_map[atom_index] = [rgb_tuple]

        add_atoms_to_highlight_map(warhead_atoms, (1.0, 0.5, 0.5))
        add_atoms_to_highlight_map(bridge_point_warhead, (1.0, 0.5, 1.0))
        add_atoms_to_highlight_map(E3_atoms, (0.5, 1.0, 0.5))
        add_atoms_to_highlight_map(bridge_point_E3ligand, (1.0, 0.5, 1.0))
        add_atoms_to_highlight_map(useless_atoms, (1.0, 1.0, 0.8))

        self.linker_highlight_map = highlight_atom_map
        # Define highlight radii
        for key in highlight_atom_map:
            highlight_radii[key] = 0.5
        self.linker_highlight_radii = highlight_radii

        # Generate molecule image with highlighted atoms
        img = self.get_molecule_image(mol, highlight_atom_map=highlight_atom_map, highlight_radii=highlight_radii)
        self.display_image(img, self.canvas)

    def rotate_image(self):
        # Rotate the RDKit molecule
        Chem.rdMolTransforms.TransformConformer(self.mol.GetConformer(), np.array(
            [[0., 1., 0., 0.], [-1., 0., 0., 0.], [0., 0., 1., 0.], [0., 0., 0., 1.]]))
        img = self.get_molecule_image(self.mol, highlight_atom_map=self.linker_highlight_map, highlight_radii=self.linker_highlight_radii)
        self.display_image(img, self.canvas)

    def flip_image(self):
        # Flip the RDKit molecule
        Chem.rdMolTransforms.TransformConformer(self.mol.GetConformer(), np.array(
            [[-1., 0., 0., 0.], [0., 1., 0., 0.], [0., 0., -1., 0.], [0., 0., 0., 1.]]))
        img = self.get_molecule_image(self.mol, highlight_atom_map=self.linker_highlight_map, highlight_radii=self.linker_highlight_radii)
        self.display_image(img, self.canvas)

    def generate_xml(self):
        if self.mol is not None:
            linker_title = filename.split("_")[0]
            warhead_atoms = self.warhead_atoms_entry.get().split()
            bridge_point_warhead = self.bridge_point_warhead_entry.get()
            E3_atoms = self.E3_atoms_entry.get().split()
            bridge_point_E3ligand = self.bridge_point_E3ligand_entry.get()
            useless_atoms = self.useless_atoms_entry.get().split()

            # Generate XML
            self.generate_xml_file(linker_title, warhead_atoms, bridge_point_warhead, E3_atoms, bridge_point_E3ligand,
                                   useless_atoms)
        else:
            messagebox.showerror("Error", "Please select a PDB file first.")

    def generate_xml_file(self, linker_title, warhead_atoms, bridge_point_warhead, E3_atoms, bridge_point_E3ligand,
                          useless_atoms):
        # Create XML structure
        root = etree.Element("linker", title=linker_title)

        # Warhead alignment atoms
        warhead_align_atoms = etree.SubElement(root, "warhead_align_atoms")
        for atom in warhead_atoms:
            atom = ''.join(re.findall(r'\d+', atom))  # extract atom index 12 from C12
            atom_name = self.mol.GetAtomWithIdx(int(atom) - 1).GetSymbol() + atom
            etree.SubElement(warhead_align_atoms, "atom").text = atom_name

        # Bridge point for warhead
        atom = ''.join(re.findall(r'\d+', bridge_point_warhead))
        atom_name = self.mol.GetAtomWithIdx(int(atom) - 1).GetSymbol() + bridge_point_warhead
        etree.SubElement(root, "bridge_point_warhead").text = atom_name

        # E3 alignment atoms
        E3_align_atoms = etree.SubElement(root, "E3_align_atoms")
        for atom in E3_atoms:
            atom = ''.join(re.findall(r'\d+', atom))  # extract atom index 12 from C12
            atom_name = self.mol.GetAtomWithIdx(int(atom) - 1).GetSymbol() + atom
            etree.SubElement(E3_align_atoms, "atom").text = atom_name

        # Bridge point for E3 ligand
        atom = ''.join(re.findall(r'\d+', bridge_point_E3ligand))
        atom_name = self.mol.GetAtomWithIdx(int(atom) - 1).GetSymbol() + bridge_point_E3ligand
        etree.SubElement(root, "bridge_point_E3ligand").text = atom_name

        # Useless atoms
        useless_atoms_elem = etree.SubElement(root, "useless_atoms")
        for atom in useless_atoms:
            atom = ''.join(re.findall(r'\d+', atom))  # extract atom index 12 from C12
            atom_name = self.mol.GetAtomWithIdx(int(atom) - 1).GetSymbol() + atom
            etree.SubElement(useless_atoms_elem, "atom").text = atom_name

        # Output XML
        xml_string = etree.tostring(root, pretty_print=True, encoding='unicode')
        self.xml_output_text.delete(1.0, tk.END)
        self.xml_output_text.insert(tk.END, xml_string)

    def populate_dropdown(self):
        if not os.path.exists(xml_file_path):
            messagebox.showerror("Error", "XML file 'PROTAC.xml' not found in the current directory.")
            return
        try:
            global protac_xml
            protac_xml = etree.parse(xml_file_path)
            print("Parsing successful")
        except etree.XMLSyntaxError:
            messagebox.showerror("Error", "Failed to parse XML file.")
            return

        warhead_options = []
        e3ligand_options = []
        # Extract warhead ligand titles
        warhead_proteins = protac_xml.xpath('//target_protein/protein')
        for protein in warhead_proteins:
            protein_title = protein.get('title')
            ligands = protein.xpath('ligand')
            for ligand in ligands:
                ligand_title = ligand.get('title')
                warhead_options.append(protein_title + '-' + ligand_title)

        # Extract E3 ligand titles
        e3ligand_proteins = protac_xml.xpath('//E3_protein/protein')
        for protein in e3ligand_proteins:
            protein_title = protein.get('title')
            ligands = protein.xpath('ligand')
            for ligand in ligands:
                ligand_title = ligand.get('title')
                e3ligand_options.append(protein_title + '-' + ligand_title)

        print(warhead_options, e3ligand_options)

        self.warhead_dropdown['menu'].delete(0, tk.END)
        self.e3ligand_dropdown['menu'].delete(0, tk.END)
        for ligand in warhead_options:
            self.warhead_dropdown['menu'].add_command(label=ligand, command=tk._setit(self.warhead_var, ligand))
        for ligand in e3ligand_options:
            self.e3ligand_dropdown['menu'].add_command(label=ligand, command=tk._setit(self.e3ligand_var, ligand))

    def process_ligand(self, selected_ligand, pdb_path, canvas, output_textbox, ligandIsWarhead):
        if selected_ligand.startswith("Select"):
            return None

        hetatm_lines = []
        pdb_atom_names = []

        with open(pdb_path + selected_ligand + ".pdb", 'r') as f:
            for line in f:
                if line.startswith('HETATM'):
                    if line[12:16].strip() != "ZN":
                        hetatm_lines.append(line)
                        # Extract the atom name from columns 13-16
                        pdb_atom_names.append(line[12:16].strip())
        print(pdb_atom_names)

        pdb_block = ''.join(hetatm_lines)
        mol = Chem.MolFromPDBBlock(pdb_block, removeHs=False)

        if mol:
            remover = SaltRemover(defnData="[Zn]")
            mol = remover.StripMol(mol, dontRemoveEverything=True)

            # Label atoms with PDB atom names
            for i, atom in enumerate(mol.GetAtoms()):
                atom.SetProp("atomNote", pdb_atom_names[i])

            # Process the ligand and assign bond orders from SMILES template or guess if SMILES is not available
            ligand_title = selected_ligand.split('-')[1]
            SMILES = protac_xml.xpath(f"//ligand[@title='{ligand_title}']/SMILES")
            if SMILES:
                SMILES_template = Chem.MolFromSmiles(SMILES[0].text)
                mol = AllChem.AssignBondOrdersFromTemplate(SMILES_template, mol)
            else:
                rdDetermineBonds.DetermineBonds(mol, charge=0)
                Chem.SanitizeMol(mol)
                mol = Chem.RemoveHs(mol)

            # Compute 2D coordinates
            AllChem.Compute2DCoords(mol)

            # Get xml atom names
            bridge_point = protac_xml.xpath(f"//ligand[@title='{ligand_title}']/bridge_point")
            bridge_point_name = [atom.text for atom in bridge_point]
            align_atoms = protac_xml.xpath(f"//ligand[@title='{ligand_title}']/align_atoms/atom")
            align_atom_names = [atom.text for atom in align_atoms]
            overlapping_atoms = protac_xml.xpath(f"//ligand[@title='{ligand_title}']/overlapping_atoms/atom")
            overlapping_atoms_names = [atom.text for atom in overlapping_atoms]
            print(bridge_point_name, align_atom_names, overlapping_atoms_names)

            # Convert atom name back to atom ID

            def atom_name_to_id(mol, atom_names):
                atom_ids = []
                for atom_name in atom_names:
                    found = False
                    for atom in mol.GetAtoms():
                        if atom.GetProp("atomNote") == atom_name:
                            atom_ids.append(atom.GetIdx())
                            found = True
                            break
                    if not found:
                        print(f"Atom with name '{atom_name}' not found in the molecule")
                return atom_ids

            bridge_point_id = atom_name_to_id(mol, bridge_point_name)
            align_atom_ids = atom_name_to_id(mol, align_atom_names)
            overlapping_atom_ids = atom_name_to_id(mol, overlapping_atoms_names)
            print(bridge_point_id, align_atom_ids, overlapping_atom_ids)

            output_textbox.delete(1.0, tk.END)
            output_textbox.insert(tk.END, f"bridge_point: {str(bridge_point_name)}\n", "pink")
            output_textbox.insert(tk.END, f"align_atoms: {str(align_atom_names)}\n", "red_or_green")
            output_textbox.insert(tk.END, f"overlapping_atoms: {str(overlapping_atoms_names)}", "yellow")

            output_textbox.tag_configure("pink", background="#ff7fff")
            if ligandIsWarhead:
                output_textbox.tag_configure("red_or_green", background="#ff7f7f")
            else:
                output_textbox.tag_configure("red_or_green", background="#7fff7f")
            output_textbox.tag_configure("yellow", background="#ffffcc")

            highlight_atom_map = {}

            def add_atoms_to_highlight_map(atom_indices, rgb_tuple):
                for atom in atom_indices:
                    if atom in highlight_atom_map:
                        highlight_atom_map[atom].append(rgb_tuple)
                    else:
                        highlight_atom_map[atom] = [rgb_tuple]

            add_atoms_to_highlight_map(bridge_point_id, (1.0, 0.5, 1.0))
            if ligandIsWarhead:
                add_atoms_to_highlight_map(align_atom_ids, (1.0, 0.5, 0.5))
            else:
                add_atoms_to_highlight_map(align_atom_ids, (0.5, 1.0, 0.5))
            add_atoms_to_highlight_map(overlapping_atom_ids, (1.0, 1.0, 0.8))

            highlight_radii = {}
            for key in highlight_atom_map:
                highlight_radii[key] = 0.5
            self.linker_highlight_radii = highlight_radii

                # Display the image on the canvas
            img = self.get_molecule_image(mol, size=400, highlight_atom_map=highlight_atom_map,
                                          highlight_radii=highlight_radii)
            self.display_image(img, canvas)

            return mol, bridge_point_id, align_atom_ids, overlapping_atom_ids
        else:
            messagebox.showerror("Error", f"Failed to load ligand from {pdb_path}")
            return None, None, None, None

    def handle_selection(self, *args):
        selected_warhead = self.warhead_var.get()
        print(selected_warhead)
        selected_e3ligand = self.e3ligand_var.get()
        if selected_warhead != "Select Warhead":
            self.warhead_mol, self.warhead_bridge_point_id, self.warhead_align_atom_ids, self.warhead_overlapping_atom_ids = \
               self.process_ligand(selected_warhead, "./warhead_ligands/", self.warhead_canvas, self.warhead_ligand_output_text, ligandIsWarhead=True)
        if selected_e3ligand != "Select E3 ligand":
            self.e3_mol, self.e3_bridge_point_id, self.e3_align_atom_ids, self.e3_overlapping_atom_ids = \
                self.process_ligand(selected_e3ligand, "./E3_ligands/", self.e3ligand_canvas, self.e3_ligand_output_text, ligandIsWarhead=False)

    def align_molecules(self):
        # Load warhead,e3, and linker molecules
        warhead_mol = self.warhead_mol
        e3_mol = self.e3_mol
        linker_mol = self.mol

        # Set MCS parameters
        # params = rdFMCS.MCSParameters()
        # params.StoreAll = True

        # Perform MCS alignment
        if warhead_mol and linker_mol is not None:
            warhead_linker_mcs_result = rdFMCS.FindMCS([warhead_mol, linker_mol], ringMatchesRingOnly=True)
            match_warhead = warhead_mol.GetSubstructMatches(warhead_linker_mcs_result.queryMol, useQueryQueryMatches=True)
            match_linker = linker_mol.GetSubstructMatches(warhead_linker_mcs_result.queryMol, useQueryQueryMatches=True)
            print("warhead match linker: " + str(match_linker))
            warhead_linker_atom_map = dict(zip(match_warhead[0], match_linker[0]))
            print("warhead - linker pairs: " + str(warhead_linker_atom_map))
        else:
            warhead_linker_atom_map = {}

        if e3_mol and linker_mol is not None:
            e3_linker_mcs_result = rdFMCS.FindMCS([e3_mol, linker_mol],  bondCompare=rdFMCS.BondCompare.CompareOrderExact)  # Optional: ringMatchesRingOnly=True
            match_e3 = e3_mol.GetSubstructMatches(e3_linker_mcs_result.queryMol, useQueryQueryMatches=True)
            match_linker = linker_mol.GetSubstructMatches(e3_linker_mcs_result.queryMol, useQueryQueryMatches=True)
            print("e3 match linker: " + str(match_linker))
            e3_linker_atom_map = dict(zip(match_e3[0], match_linker[0]))
            print("e3 - linker pairs: " + str(e3_linker_atom_map))
        else:
            e3_linker_atom_map = {}

        linker_warhead_align_atom_id = []
        linker_e3_align_atom_id = []
        linker_overlapping_atom_id = []
        try:
            print("warhead align id: "+str(self.warhead_align_atom_ids))
            for atom in self.warhead_align_atom_ids:
                linker_warhead_align_atom_id.append(warhead_linker_atom_map[atom])
        except KeyError:
            messagebox.showerror("Error", f"Unable to find warhead align atoms. Maybe the maximum common structure doesn't align well.\n"
                                          f"warhead - linker pairs: {str(warhead_linker_atom_map)}")

        try:
            for atom in self.e3_align_atom_ids:
                linker_e3_align_atom_id.append(e3_linker_atom_map[atom])
        except KeyError:
            messagebox.showerror("Error", f"Unable to find E3 ligand align atoms. Maybe the maximum common structure doesn't align well.\n"
                                          f"E3 - linker pairs: {str(e3_linker_atom_map)}")

        linker_overlapping_atom_id = list(warhead_linker_atom_map.values()) + list(e3_linker_atom_map.values())
        linker_overlapping_atom_id_to_remove = []
        for atom in self.warhead_overlapping_atom_ids:
            try:
                linker_overlapping_atom_id_to_remove.append(warhead_linker_atom_map[atom])
            except KeyError:
                pass
        for atom in self.e3_overlapping_atom_ids:
            try:
                linker_overlapping_atom_id_to_remove.append(e3_linker_atom_map[atom])
            except KeyError:
                pass
        linker_overlapping_atom_id = [x for x in linker_overlapping_atom_id if x not in linker_overlapping_atom_id_to_remove]

        # output to textbox
        self.warhead_atoms_entry.delete(0, tk.END)
        self.warhead_atoms_entry.insert(tk.END, " ".join(str(atom+1) for atom in linker_warhead_align_atom_id))
        self.E3_atoms_entry.delete(0, tk.END)
        self.E3_atoms_entry.insert(tk.END, " ".join(str(atom + 1) for atom in linker_e3_align_atom_id))
        self.useless_atoms_entry.delete(0, tk.END)
        self.useless_atoms_entry.insert(tk.END, " ".join(str(atom+1) for atom in linker_overlapping_atom_id))
        self.bridge_point_warhead_entry.delete(0, tk.END)
        self.bridge_point_E3ligand_entry.delete(0, tk.END)

        highlight_atom_map = {}

        def add_atoms_to_highlight_map(atom_indices, rgb_tuple):
            for atom in atom_indices:
                if atom in highlight_atom_map:
                    highlight_atom_map[atom].append(rgb_tuple)
                else:
                    highlight_atom_map[atom] = [rgb_tuple]

        add_atoms_to_highlight_map(linker_warhead_align_atom_id, (1.0, 0.5, 0.5))
        add_atoms_to_highlight_map(linker_e3_align_atom_id, (0.5, 1.0, 0.5))
        add_atoms_to_highlight_map(linker_overlapping_atom_id, (1.0, 1.0, 0.8))

        self.linker_highlight_map = highlight_atom_map

        print(highlight_atom_map)

        highlight_radii = {}
        for key in highlight_atom_map:
            highlight_radii[key] = 0.5

        self.linker_highlight_radii = highlight_radii

        # Generate molecule image with highlighted atoms
        img = self.get_molecule_image(linker_mol, highlight_atom_map=highlight_atom_map, highlight_radii=highlight_radii)
        self.display_image(img, self.canvas)

    def clear_input(self):
        self.warhead_atoms_entry.delete(0, tk.END)
        self.E3_atoms_entry.delete(0, tk.END)
        self.useless_atoms_entry.delete(0, tk.END)
        self.bridge_point_warhead_entry.delete(0, tk.END)
        self.bridge_point_E3ligand_entry.delete(0, tk.END)


def main():
    root = tk.Tk()
    app = AtomSelectorApp(root)
    root.mainloop()


if __name__ == "__main__":
    main()
