import webbrowser, copy

import pathlib

import pandas as pd
import numpy as np

from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.Draw import rdMolDraw2D
from IPython.display import SVG
from rdkit.Chem import rdDepictor
from rdkit.Geometry.rdGeometry import Point2D

import ipywidgets as widgets
from traitlets import Unicode, Int, validate
from IPython.display import display, HTML


class MolSVGWidget(widgets.DOMWidget):
    """Gregs MolSVGWidget that joins with the javascript"""
    _view_name = Unicode('MolSVGView').tag(sync=True)
    _view_module = Unicode('molsvg_widget').tag(sync=True)
    _view_module_version = Unicode('0.0.1').tag(sync=True)
    
    svg = Unicode('', help="svg to be rendered").tag(sync=True)
    clicked_atom_idx = Unicode('', help="The index of the atom that was just clicked").tag(sync=True)
    clicked_bond_idx = Unicode('', help="The index of the bond that was just clicked").tag(sync=True)
    
class EditMol(object):
    def __init__(self, mol = Chem.MolFromSmiles("c1c([NH3+])cccc1[C@H](C)C(=O)OCC=CC")):
        #Load the javascript needed for MolSVGwidget
        self.inject_javascript()

        style = {'description_width': 'initial'}
        #Create the outputs and widgets
        self.o_mol = widgets.Output()  #Molecule view
        self.o_log = widgets.Output(layout=widgets.Layout(width='100%', height='160px', overflow_y='auto')) #For debug
        self.o_status = widgets.Output(layout=widgets.Layout(width='50%', height='80px', overflow_y='auto')) #Status for Users
        
        self.o_atomclicked = widgets.Text(description="Index of clicked atom",
                                         #layout = widgets.Layout(width="100px"),                                       
                                         style=style)
        self.o_bondclicked = widgets.Text(description="Index of clicked bond",
                                         #layout = widgets.Layout(width="100px"),                                       
                                         style=style)
        
        

        tglb_style = widgets.ToggleButtonsStyle(button_width = '30px')

        self.tglb_atomtype = widgets.ToggleButtons(
                                        options=['C', 'N', 'O', 'S', 'P', 'F', 'Cl','Br','I'],
                                        description='AtomType',
                                        disabled=False,
                                        button_style='', # 'success', 'info', 'warning', 'danger' or ''
                                        tooltips=['Hotbar of commonly used atom types'],
                                    #     icons=['check'] * 
                                        style=tglb_style,
                                        #layout=,
                                        button_width='10px',
                                    )

        self.dd_atomtype = widgets.Dropdown(options=['H', 'Li','Be','B', 'C', 'N', 'O', 'F','Na','Mg','Al','Si','P','S', 'Cl',
                                                    'K','Ca','Br','I'], #TODO full list of Atoms
                                            value='C',
                                            description='AtomType')
        
        self.tglb_atomtype.observe(self.changeAtomType, names="value")
        self.dd_atomtype.observe(self.changeAtomType, names="value")

        self.tglb_bondtype = widgets.ToggleButtons(
                                        options=['-', '=', '#'],
                                        description='BondType',
                                        disabled=False,
                                        button_style='', # 'success', 'info', 'warning', 'danger' or ''
                                        tooltips=['Single','Double','Triple'],
                                        style=tglb_style,
                                    #     icons=['check'] * 3
                                    )

        self.dd_bondtype = widgets.Dropdown(options=Chem.rdchem.BondType.names.keys(),
                                            value='SINGLE',
                                            description='BondType')

        self.tglb_bondtype.observe(self.changeBondType, names="value")
        self.dd_bondtype.observe(self.changeBondType, names="value")

                                    

        #Icons can be from some of these: https://fontawesome.com/icons?d=gallery

        self.tglb_actions = widgets.ToggleButtons(
                                        options=['Select ', 'Add ', 'Replace ', 'Remove ',
                                                   'Add Bond ',
                                                   "Increase Charge ", "Decrease Charge ",
                                                  "RS-toggle", "EZ-toggle"],
                                        description='Actions',
                                        disabled=False,
                                        button_style='', # 'success', 'info', 'warning', 'danger' or ''
                                        tooltips=['Hotbar of actions'],
                                        #style=,
                                        icons=['mouse-pointer', 'pencil', 'repeat', 'times-circle',
                                                   'link',
                                                   "plus-circle", "minus-circle",
                                                  "RS-toggle","EZ-toggle"],
                                        value='Add ',
                                    )
        
        
        
        self.dd_action = widgets.Dropdown(options=['Select', 'Add', 'Replace', 'Remove',
                                                   'Add Bond',
                                                   "Increase Charge", "Decrease Charge",
                                                  "RS-toggle", "EZ-toggle"],
                                            value='Add',
                                            description='Action')
        self.btn_undo = widgets.Button(description="Undo ", icon="undo")
        self.btn_undo.on_click(self.undo)

        self.btn_about = widgets.Button(description="About")
        self.btn_about.on_click(self.about)
        
        self.btn_clean = widgets.Button(description="Clean ", icon="broom")
        self.btn_clean.on_click(self.clean_mol_coords)
        

        #Lookup tables for bondtypes
        self.bondtypes = copy.copy(Chem.rdchem.BondType.names) #Need to copy, or will change the existing ones
        self.bondtypes.update({ "-": Chem.rdchem.BondType.SINGLE,
                                "=": Chem.rdchem.BondType.DOUBLE,
                                "#": Chem.rdchem.BondType.TRIPLE,
                            })

        self.bondtypes_inv = {i:k for k,i in Chem.rdchem.BondType.names.items()} #Use original in inverse lookup


        self._selectedAtoms = []

        
        #Build the interface from the widgets
        #First Menu of buttons
        display(widgets.HBox([self.dd_atomtype, self.dd_bondtype, self.btn_undo, self.btn_clean, self.btn_about]))
        #Second menu with the quick actions
        display(self.tglb_actions)
        display(widgets.HBox([self.tglb_atomtype, self.tglb_bondtype]))
        #Main Window
        display(self.o_mol)
        #Status bar
        display(widgets.HBox([self.o_atomclicked, self.o_bondclicked, self.o_status]))
        #Debug and such
        display(widgets.HBox([self.o_log]))
        
        #Create the SVG widget
        self.create_widget()
        
        #Set the mol
        self._prevmol = None #For undo
        self._mol = None
        self.mol = mol

    def about(self, _):
        webbrowser.open("https://www.cheminformania.com/rdeditor-an-open-source-molecular-editor-based-using-python-pyside2-and-rdkit/") #TODO update with URL for blog-post

    def inject_javascript(self):
        filepath = pathlib.Path(__file__).parent.absolute()
        javascript = open(f"{filepath}/molsvgwidget.js", "r").read()
        script = f"""<script>{javascript}</script>"""
        display(HTML(script))

    def changeAtomType(self, e):
        """Function to update the quick action tgl_buttons if dropdown selected, and vice versa"""
        #TODO, can this be linked somehow using the widgets framework?
        newatomtype = e["new"]
        self.log(newatomtype)
        #Set the tgl_buttons to new type
        if newatomtype in self.tglb_atomtype.options:
            self.tglb_atomtype.value = newatomtype
        else:
            self.tglb_atomtype.value = None
        #Set dropdown
        if newatomtype in self.dd_atomtype.options:
            self.dd_atomtype.value = newatomtype

    def changeBondType(self, e):
        """Function to update the quick action tgl_buttons if dropdown selected, and vice versa"""
        newbondtype = e["new"]
        self.log(newbondtype)
        #Issue, - = # in tglbuttons, but RDKit names in dd menu TODO review code, selecting single from DD doesn't update tgl_button
        if newbondtype in self.tglb_bondtype.options:
            #self.log("Here!")
            rdbondtype = self.bondtypes[newbondtype]
            textbondtype = self.bondtypes_inv[rdbondtype]
            #self.log("TBtype: %s"%textbondtype)
            self.dd_bondtype.value = textbondtype
            self.tglb_bondtype.value = newbondtype
        else:
            self.tglb_bondtype.value = None
        if newbondtype in self.dd_bondtype.options:
            self.dd_bondtype.value = newbondtype
        

    @property
    def atomType(self):
        return self.dd_atomtype.value

    @property
    def bondType(self):
        return self.dd_bondtype.value    

    @property
    def selectedAtoms(self):
        return self._selectedAtoms

    @selectedAtoms.setter
    def selectedAtoms(self, atomlist):
        if atomlist != self._selectedAtoms:
            assert type(atomlist) == list, "selectedAtoms should be a list of integers"
            assert all(isinstance(item, int) for item in atomlist), "selectedAtoms should be a list of integers"
            self._selectedAtoms = atomlist
            self.selection_changed()

    def selectAtom(self, atomidx):
        if not atomidx in self._selectedAtoms:
            self._selectedAtoms.append(atomidx)
            self.selection_changed()

    def unselectAtom(self, atomidx):
        self._selectedAtoms.remove(atomidx)
        self.selection_changed()

    def clearAtomSelection(self):
        if self._selectedAtoms != []:
            self._selectedAtoms = []
            self.selection_changed()            

    def selection_changed(self):
        svg = self.generate_svg()
        self.molwidget.svg = svg
        
        
    def computeNewCoords(self, ignoreExisting=False, canonOrient=False):
        """Computes new coordinates for the molecule taking into account all
        existing positions (feeding these to the rdkit coordinate generation as
        prev_coords).
        """
        # This code is buggy when you are not using the CoordGen coordinate
        # generation system, so we enable it here
        rdDepictor.SetPreferCoordGen(True)
        prev_coords = {}
        if self._mol.GetNumConformers() == 0:
            self.log("No Conformers found, computing all 2D coords")
        elif ignoreExisting:
            self.log("Ignoring existing conformers, computing all 2D coords")
        else:
            assert self._mol.GetNumConformers() == 1
            self.log("1 Conformer found, computing 2D coords not in found conformer")
            conf = self._mol.GetConformer(0)
            for a in self._mol.GetAtoms():
                pos3d = conf.GetAtomPosition(a.GetIdx()) #If an atom has the coordinate 0,0, we assume it is a new one
                if (pos3d.x, pos3d.y) == (0, 0):
                    continue
                prev_coords[a.GetIdx()] = Point2D(pos3d.x, pos3d.y)
        self.log("Coordmap %s"%prev_coords)
        self.log("canonOrient %s"%canonOrient)
        rdDepictor.Compute2DCoords(self._mol, coordMap=prev_coords, canonOrient=canonOrient)
        
    def clean_mol_coords(self, evt):
        self.computeNewCoords(ignoreExisting=True, canonOrient=True)
        self.sanitizeMol()
        svg = self.generate_svg()
        self.molwidget.svg = svg        


    def generate_svg(self):
        """Generate the svg needed for updating the molwidget"""
        
        self.drawer = rdMolDraw2D.MolDraw2DSVG(600,300)
        #TODO, what if self._drawmol doesn't exist?
        if self._drawmol != None:
            #Chiral tags on R/S
            chiraltags = Chem.FindMolChiralCenters(self._drawmol)
            opts = self.drawer.drawOptions()
            for tag in chiraltags:
                idx = tag[0]
                opts.atomLabels[idx]= self._drawmol.GetAtomWithIdx(idx).GetSymbol() + ':' + tag[1]
            self.log(len(self._selectedAtoms))
            if len(self._selectedAtoms) > 0: #TODO, if a selected atom is deleted, the draw function stops working!
                self.log("Highlight Atoms")
                #Bond colors can't be made a different color?
                #colors = {atomidx:(0.8,0.8,0.8) for atomidx in self._selectedAtoms}
                #colors[self._selectedAtoms[-1]] = (0.6,0.6,0.6)
                #bondcolors = {bondidx:(0.8,0.8,0.8) for bondidx in range(len(self._drawmol.GetBonds()))}
                colors={self._selectedAtoms[-1]:(1,0.2,0.2)} #Color lastly selected a different color
                self.drawer.DrawMolecule(self._drawmol, highlightAtoms=self._selectedAtoms, highlightAtomColors=colors )#, highlightBondColors=bondcolors)
            else:
                self.drawer.DrawMolecule(self._drawmol)
        self.drawer.TagAtoms(self._drawmol)
        self.drawer.FinishDrawing()
        svg = self.drawer.GetDrawingText()
        return svg

    def create_widget(self):
        """Create the interactive SVG mol widget"""
        self.molwidget = MolSVGWidget(svg='')
        self.draw_widget()
        self.create_observer()
        
    def draw_widget(self):
        """Display the mol widget"""
        self.o_mol.clear_output()
        with self.o_mol:
            display(self.molwidget)

    def create_observer(self):
        """Create the observers that should react to the clicked event"""
        self.molwidget.observe(self.on_atom_clicked, names="clicked_atom_idx")
        self.molwidget.observe(self.on_bond_clicked, names="clicked_bond_idx")

    def update_widget(self):
        """Update SVG of the widget"""
        self.computeNewCoords()
        svg = self.generate_svg()
        self.molwidget.svg = svg

    @property
    def mol(self):
        """Return the private mol"""
        return self._mol

    @mol.setter
    def mol(self, mol):
        """Set the private mol and initalize interactive SVG and update output widgets"""
        if mol == None:
            mol=Chem.MolFromSmiles('') #Create an empty mol
        if mol != self._mol:
            #TODO assert that this is a RDKit mol
            if self._mol != None:
                self._prevmol = Chem.Mol(self._mol.ToBinary()) #Copy mol and assign
            self._mol = mol
        self.computeNewCoords()
        self.sanitizeMol()
        self.update_widget()
        
    def log(self, message):
        with self.o_log:
            print(message)
            
    def status(self, message):
        self.o_status.clear_output()
        with self.o_status:
            print(message)

    def sanitizeMol(self, kekulize=True, drawkekulize=False):
        #self.computeNewCoords()
        #self._mol.UpdatePropertyCache() Updatepropertycache leads to sanitazion errors (did it always do that?)
        self._drawmol = Chem.Mol(self._mol.ToBinary()) #Copy mol to make changes
        try:
            Chem.SanitizeMol(self._drawmol)
            self.status("Sanitizable")
        except:
            #self.sanitizeSignal.emit("UNSANITIZABLE")
            self.status("UNSANITIZABLE")
            #self.logger.warning("Unsanitizable")
            try:
                self._drawmol.UpdatePropertyCache(strict=False)
            except:
                #self.sanitizeSignal.emit("UpdatePropertyCache FAIL")
                #self.logger.error("Update Property Cache failed")
                self.log("Update Property Cache failed")
        #Kekulize
        if kekulize:
            try:
                Chem.Kekulize(self._drawmol)
            except:
                
                self.log("Unkekulizable")
        try:
            self._drawmol = rdMolDraw2D.PrepareMolForDrawing(self._drawmol, kekulize=drawkekulize)
        except ValueError:  # Can happen on a kekulization failure
            self._drawmol = rdMolDraw2D.PrepareMolForDrawing(self._drawmol, kekulize=False)

    def undo(self, e):
        self.mol = self._prevmol

    def on_atom_clicked(self, e):
        """Callback for reacting to atom clicked"""
        with self.o_log: #Error messages will go here
            
            if e["new"] == "Event":#Bogus event. If value is not changed, event is not triggered, thus atom can't be clicked twice in a row
                self.log('Non event')
                return
            self.o_atomclicked.value = e["new"]
            self.log("Atom Clicked!")
            atomidx = int(e["new"])

            #Get action and delegate to function
            action = self.tglb_actions.value.strip()

            if action == "Add":
                self.add_atom_to(atomidx)
            elif action == "Remove":
                self.remove_atom(atomidx)
            elif action == "Select":
                self.select_atom_click(atomidx)
            elif action == "Replace":
                self.replace_atom(atomidx)
            elif action == "Add Bond":
                self.add_bond(atomidx)
            elif action == "Increase Charge":
                self.increase_charge(atomidx)
            elif action == "Decrease Charge":
                self.decrease_charge(atomidx)
            elif action == "RS-toggle":
                self.toogle_RS(atomidx)
            else:
                self.log("The combination of Atom click and Action %s undefined"%self.action)

    
    def add_atom_to(self, atomidx):
        rwmol = Chem.rdchem.RWMol(self.mol)
        newatom = Chem.rdchem.Atom(self.tglb_atomtype.value)
        newidx = rwmol.AddAtom(newatom)
        bondtype = self.bondtypes[self.tglb_bondtype.value]
        newbond = rwmol.AddBond(atomidx, newidx, order=bondtype)
        self.mol = rwmol
        
    def remove_atom(self, atomidx):
        rwmol = Chem.rdchem.RWMol(self.mol)
        rwmol.RemoveAtom(atomidx)
        #self.clearAtomSelection() # Removing atoms updates Idx'es
        self.mol = rwmol
        
    def replace_atom(self, atomidx):
        #Update atom properties with the text values from the widgets
        atomtype = self.tglb_atomtype.value
        prevatom = self.mol.GetAtomWithIdx(atomidx)
        rwmol = Chem.rdchem.RWMol(self.mol)
        newatom = Chem.rdchem.Atom(atomtype)
        newatom.SetIsAromatic(prevatom.GetIsAromatic())
        
        rwmol.ReplaceAtom(atomidx, newatom)
        #rwmol.UpdatePropertyCache()
        self.mol = rwmol
        
    def increase_charge(self, atomidx):
        #self.backupMol()
        atom = self.mol.GetAtomWithIdx(atomidx)
        atom.SetFormalCharge(atom.GetFormalCharge()+1)
        self.mol = self.mol
        
    def decrease_charge(self, atomidx):
            #self.backupMol()
            atom = self.mol.GetAtomWithIdx(atomidx)
            atom.SetFormalCharge(atom.GetFormalCharge()-1)
            self.mol = self.mol
            
 #   def select_atom_bak(self, atomidx):
 #       if atomidx in self._selectedAtoms:
 #           self.unselectAtom(atomidx)
 #       else:
 #           self.selectAtomAdd(atomidx)
    
            
    def select_atom_click(self, atomidx):
        if atomidx in self._selectedAtoms:
            self.unselectAtom(atomidx)
        else:
            self.selectAtom(atomidx)    
        
    #Double step action
    def add_bond(self, atomidx):
        if len(self.selectedAtoms) > 0:
            selected = self.selectedAtoms[-1]
            rwmol = Chem.rdchem.RWMol(self.mol)
            rwmol.AddBond(selected, atomidx, order=self.bondtypes[self.tglb_bondtype.value])
            self.mol = rwmol
        else:
            self.selectAtom(atomidx)

        
    def toogle_RS(self, atomidx):
        #self.backupMol()
        atom = self._mol.GetAtomWithIdx(atomidx)
        stereotype = atom.GetChiralTag()
        #self.log("Current stereotype of clicked atom %s"%stereotype)
        stereotypes = [Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CCW,
#                        Chem.rdchem.ChiralType.CHI_OTHER, this one doesn't show a wiggly bond
                        Chem.rdchem.ChiralType.CHI_UNSPECIFIED,
                        Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CW,
                        Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CCW]
        newidx = np.argmax(np.array(stereotypes) == stereotype)+1
        atom.SetChiralTag(stereotypes[newidx])
        #self.log("New stereotype set to %s"%atom.GetChiralTag())
        #rdDepictor.Compute2DCoords(self._mol)
        #self._mol.ClearComputedProps()
        self._mol.UpdatePropertyCache()
        #rdDepictor.Compute2DCoords(self._mol)
        #self.molChanged.emit()
        self.mol = self._mol


    def on_bond_clicked(self, e):
        with self.o_log: #Error messages will go here
            """Callback for reacting to a bond click"""
            if e["new"] == "Event":#Bogus event, to trigger change. Default behaviour is, that if value is not changed, event is not triggered, thus atom can't be clicked twice in a row
                self.log('Non event')
                return
            self.o_bondclicked.value = e["new"]
            self.log("Bond Clicked!")
            bondidx = int(e["new"])

            #Get action and delegate to function
            action = self.tglb_actions.value.strip()

            if action == "Add":
                self.toggle_bond(bondidx)
            elif action == "Add Bond":
                self.replace_bond(bondidx)
            elif action == "Remove":
                self.remove_bond(bondidx)
            elif action == "Select":
                self.select_bond_click(bondidx)
            elif action == "Replace":
                self.replace_bond(bondidx)
            elif action == "EZ-toggle":
                self.toggle_EZ(bondidx)
            else:
                self.log("The combination of bond click and Action %s undefined"%self.action)

    def replace_bond(self, bondidx):
        editmol = Chem.Mol(self._mol.ToBinary()) #Copy mol to make changes
        bond = editmol.GetBondWithIdx(bondidx)
        bondtype = self.bondtypes[self.dd_bondtype.value]
        bond.SetBondType(bondtype)
        self.mol = editmol
        
    def toggle_bond(self, bondidx):
        bondtype = self.mol.GetBondWithIdx(bondidx).GetBondType()
        bondtypes = [Chem.rdchem.BondType.TRIPLE,
                    Chem.rdchem.BondType.SINGLE,
                    Chem.rdchem.BondType.DOUBLE,
                    Chem.rdchem.BondType.TRIPLE]
        #Find the next type in the list based on current
        #If current is not in list => Then it selects the first and add 1 => SINGLE
        newidx = np.argmax(np.array(bondtypes) == bondtype)+1
        newtype = bondtypes[newidx]       
        editmol = Chem.Mol(self._mol.ToBinary()) #Copy mol to make changes
        bond = editmol.GetBondWithIdx(bondidx)
        bond.SetBondType(newtype)
        self.mol = editmol
    
    def remove_bond(self, bondidx):
        rwmol = Chem.rdchem.RWMol(self.mol)
        bond = rwmol.GetBondWithIdx(bondidx)
        rwmol.RemoveBond(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())
        self.mol = rwmol
        
    def select_bond_click(self, bondidx):
        self.log("Bond selection not supported")
        

    def assert_stereo_atoms(self, bond):
        if len(bond.GetStereoAtoms()) ==0:
            #get atoms and idx's of bond
            bondatoms = [bond.GetBeginAtom(), bond.GetEndAtom()]
            bondidx = [atom.GetIdx() for atom in bondatoms]
            
            #Figure out the atom idx's of the neigbor atoms, that are NOT the other end of the bond
            stereoatoms = []
            for bondatom in bondatoms:
                neighboridxs = [atom.GetIdx() for atom in bondatom.GetNeighbors()]
                neighboridx = [idx for idx in neighboridxs if idx not in bondidx][0]
                stereoatoms.append(neighboridx)
            #Set the bondstereoatoms
            bond.SetStereoAtoms(stereoatoms[0], stereoatoms[1])
        else:
            pass
        
       
    
    def toggle_EZ(self, bondidx):
        #self.backupMol()
        #Chem.rdmolops.AssignStereochemistry(self._mol,cleanIt=True,force=False)
        editmol = Chem.Mol(self._mol.ToBinary())
        bond = editmol.GetBondWithIdx(bondidx)
        stereotype = bond.GetStereo()
        self.assert_stereo_atoms(bond)
        #self.logger.debug("Current stereotype of clicked atom %s"%stereotype)
        #TODO: what if molecule already contain STEREOE or STEREOZ
        stereotypes = [Chem.rdchem.BondStereo.STEREONONE,
                        Chem.rdchem.BondStereo.STEREOCIS,
                        Chem.rdchem.BondStereo.STEREOTRANS,
                        #Chem.rdchem.BondStereo.STEREOANY, TODO, this should be crossed, but is not
                        Chem.rdchem.BondStereo.STEREONONE,]
        newidx = np.argmax(np.array(stereotypes) == stereotype)+1
        
        bond.SetStereo(stereotypes[newidx])
        editmol.ClearComputedProps()
        #Chem.rdmolops.AssignStereochemistry(self._mol,cleanIt=True,force=False)
        editmol.UpdatePropertyCache()
        self.mol = editmol
  