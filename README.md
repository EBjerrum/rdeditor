# rdeditor

Simple RDKit molecule editor GUI using PySide6
![rdeditor, the RDKit molecule editor](./Screenshots/Main_window.png)

## Installation

- requirements

RDKit, NumPy, PySide6 and pyqtdarktheme should be automatically pip installed by the setup.py script.

- installation

```bash
pip install rdeditor

```

A launch script will also be added so that it can be started from the command line via the `rdEditor` command.

## Usage

Can be started with `rdEditor` or `rdEditor your_molecule.mol` to start edit an existing molecule.
Interactions with the molecule are done via clicking on the canvas, atoms or bonds. A choice of tools is available.

#### Top Menu:

![top menu of rdeditor, the RDKit molecule editor](./Screenshots/Top_Menu.png)

From left to right

- Open: Open a molfile
- Save: Save current molecule
- Save As: Save current molecule with a new name

- Arrow: Select tool. Click on an atom to select it, click on the canvas to deselect. Clicking on multiple atoms one after another will select them, but only the lastly clicked one will be highlighted in red and used for operations, such as bond creation to another existing atom.
- Pen: Add tool. Clicking on an existing atom will add the current selected atom type to that atom with a single bond. Clicking on the canvas will add a disconnected atom. Clicking on a bond will cycle through single, double and triple bonds.
- Add bond / Join atoms: Will add a single bond between a clicked atom (or a previously selected atom) and the next atom clicked.
- Change Atom: Will substitute the atom clicked, with the currently selected atom type
- R/S: Change the stereo chemistry of the selected atom (see issues below)
- E/Z: Change E/Z stereo of double bonds (see issues below)
- Increase/Decrease charge: Will increase or decrease the charge of the atom clicked
- Delete atom/bond:
- Clear Canvas
- Undo.

#### Side Bar:

![top menu of rdeditor, the RDKit molecule editor](./Screenshots/Side_bar.png)

Most commonly used bond types, and atom types can be selected. A Periodic table is accessible for exotic atom types.

#### Dropdown menus

Access to all standard operations as well as less used atom types and bond-types.

#### Settings

Themes can be selected from the ones available on your platform (Mac/Linux/Windows)

The debug level can be selected

## Development

Instructions to set it up in editable modes and instructions for eventual contributions can be found in the [DEVELOPER.md](./DEVELOPER.md) file.
Please reach out first, there may be a relevant development branch.

## Additional Reading

I wrote a blog post with an overview of the structure of the code.
[https://www.wildcardconsulting.dk/rdeditor-an-open-source-molecular-editor-based-using-python-pyside2-and-rdkit/](https://www.wildcardconsulting.dk/rdeditor-an-open-source-molecular-editor-based-using-python-pyside2-and-rdkit/)

We also published a preprint on ChemRxiv: [https://chemrxiv.org/engage/chemrxiv/article-details/65e6dcfa9138d23161b2979c](https://chemrxiv.org/engage/chemrxiv/article-details/65e6dcfa9138d23161b2979c)

## ISSUES

Please report issues at GitHub, it's tough getting all corners of a GUI tested.
