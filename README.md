# rdeditor

Simple RDKit molecule editor GUI using PySide6

![rdeditor, the RDKit molecule editor](https://github.com/EBjerrum/rdeditor/blob/master/Screenshots/Main_window.png?raw=true)

## Installation

- requirements

RDKit, NumPy, PySide6 and pyqtdarktheme should be automatically installed when installed via pip.

- installation

```bash
pip install rdeditor

```

A launch script will also be added so that it can be started from the command line via the `rdEditor` command.

## Usage

Can be started from the command line with `rdEditor` or `rdEditor your_molecule.mol` to start edit an existing molecule.
Interactions with the molecule are done via clicking and dragging on the canvas, atoms or bonds. A choice of tools is available.

To edit a molecule, select the pen tool, and an atom, bond or template type and click on the canvas to add it.

Clicking:

- When clicking existing atoms or bonds, the clicked atom or bond will be modified, depending on the atom or bond type selected.
- If a bondtype is selected, the bond will be added to the clicked atom with a carbon atom at the other end.
- If a bondtype is selected and a bond is clicked, the bond will be changed to that type if different.
- If you click multiple times on a bond with an atomtype selected, the bondtype will cycle between single, double and triple bond.

Dragging:

- If you click and drag from an atom, the selected atomtype will be added at the end of a single bond.
- Dragging between two atoms will add a bond between them. If a bondtype is selected, the bondorder will correspond to the bondtype selected.
- Dragging on the canvas with an atomtype will add two atoms of that type with a single bond between them.
- Dragging on the canvas with a bondtype will simply add that bond with carbon atoms at both ends.

Templates:

- Templates work kind of like atoms, so if you click on an atom, the template will be added directly to that atom.
- If a template is selected and dragged from an atom, the template will be added with a single bond to the clicked atom.
- Some templates can also be added to bonds by clicking on the middle of the bond.
- Dragging on a canvas with a template will add a carbon atom and a single bond to the template.

Other actions:

- most other actions (R/S, E/Z, Increase/Decrease charge, Adjust atom number) works by clikcing on existing atoms.

#### Top Menu:

![top menu of rdeditor, the RDKit molecule editor](https://github.com/EBjerrum/rdeditor/blob/master/Screenshots/Top_Menu.png?raw=true)

From left to right

- Open: Open a molfile
- Save: Save current molecule
- Save As: Save current molecule with a new name

- Arrow: Select tool. Click on an atom to select it, click on the canvas to deselect. Clicking on multiple atoms one after another will select them, but only the lastly clicked one will be highlighted in red and used for operations, such as bond creation to another existing atom.
- Pen: Add tool. Clicking on an existing atom will add the current selected atom type to that atom with a single bond. Clicking on the canvas will add a disconnected atom. Clicking on a bond will cycle through single, double and triple bonds.
- R/S: Change the stereo chemistry of the selected atom
- E/Z: Change E/Z stereo of double bonds
- Increase/Decrease charge: Will increase or decrease the charge of the atom clicked
- Set atommap or R-group number: Will set the atommap or R-group number of the atom clicked
- clean up coordinates: recalculate coordinates disregarding existing coordinates.
- clean up chemistry. Sanitize and/or Kekulize the molecule.
- Delete atom/bond:
- Clear Canvas
- Undo.

#### Side Bar:

![top menu of rdeditor, the RDKit molecule editor](https://github.com/EBjerrum/rdeditor/blob/master/Screenshots/Side_bar.png?raw=true)

Most commonly used bond types, and atom types can be selected. Templates and R-group (dummy atoms) are also accessible. A Periodic table is accessible for exotic atom types.

#### Dropdown menus

Access to all standard operations as well as less used atom types and bond-types.

#### Settings

Themes can be selected from the ones available on your platform (Mac/Linux/Windows).

The debug level can be selected

Cleanup settings can be selected if the molecule should be sanitized or kekulized during cleanup.

## Development

Instructions to set it up in editable modes and instructions for eventual contributions can be found in the [DEVELOPER.md](./DEVELOPER.md) file.
Please reach out first, there may be a relevant development branch.

## Additional Reading

I wrote a blog post with an overview of the structure of the code.
[https://www.wildcardconsulting.dk/rdeditor-an-open-source-molecular-editor-based-using-python-pyside2-and-rdkit/](https://www.wildcardconsulting.dk/rdeditor-an-open-source-molecular-editor-based-using-python-pyside2-and-rdkit/)

We also published a preprint on ChemRxiv: [https://chemrxiv.org/engage/chemrxiv/article-details/65e6dcfa9138d23161b2979c](https://chemrxiv.org/engage/chemrxiv/article-details/65e6dcfa9138d23161b2979c)

## ISSUES

Please report issues at GitHub, it's tough getting all corners of a GUI tested.
