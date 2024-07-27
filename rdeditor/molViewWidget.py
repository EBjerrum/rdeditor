#!/usr/bin/python
# Import required modules
from __future__ import print_function
from PySide6 import QtCore, QtGui, QtSvg, QtWidgets, QtSvgWidgets
import sys
import copy

# from types import *
import logging

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Geometry.rdGeometry import Point2D

from rdeditor.utilities import validate_rgb


# The Viewer Class
class MolWidget(QtSvgWidgets.QSvgWidget):
    def __init__(self, mol=None, parent=None):
        # Also init the super class
        super(MolWidget, self).__init__(parent)

        # logging
        logging.basicConfig()
        self.logger = logging.getLogger()
        self.loglevel = logging.WARNING

        # This sets the window to delete itself when its closed, so it doesn't keep lingering in the background
        self.setAttribute(QtCore.Qt.WA_DeleteOnClose)
        # Private Properties
        self._mol = None  # The molecule
        self._drawmol = None  # Molecule for drawing
        self.drawer = None  # drawing object for producing SVG
        self._selectedAtoms = []  # List of selected atoms
        self._darkmode = False

        # Color settings
        self._unsanitizable_background_colour = None  # (1, 0.75, 0.75)
        self._last_selected_highlight_colour = (1, 0.2, 0.2)
        self._selected_highlight_colour = (1, 0.5, 0.5)

        # Sanitization Settings
        self._kekulize = False
        self._sanitize = False
        self._updatepropertycache = False

        # Bind signales to slots for automatic actions
        self.molChanged.connect(self.sanitize_draw)
        self.selectionChanged.connect(self.draw)
        self.drawSettingsChanged.connect(self.draw)
        self.sanitizeSignal.connect(self.changeSanitizeStatus)

        # Initialize class with the mol passed
        self.mol = mol

    ##Properties and their wrappers
    @property
    def loglevel(self):
        return self.logger.level

    @loglevel.setter
    def loglevel(self, loglvl):
        self.logger.setLevel(loglvl)

    @property
    def darkmode(self):
        return self._darkmode

    @darkmode.setter
    def darkmode(self, value: bool):
        self._darkmode = bool(value)
        self.draw()

    # Getter and setter for mol
    molChanged = QtCore.Signal(name="molChanged")

    @property
    def mol(self):
        return self._mol

    @mol.setter
    def mol(self, mol):
        if mol is None:
            mol = Chem.MolFromSmiles("")
        if mol != self._mol:
            # assert isinstance(mol, Chem.Mol)
            if self._mol is not None:
                self._prevmol = copy.deepcopy(self._mol)  # Chem.Mol(self._mol.ToBinary())  # Copy

            # # TODO make this failsafe
            # if self._updatepropertycache:
            #     try:
            #         mol.UpdatePropertyCache(strict=False)
            #     except Exception as e:
            #         self.sanitizeSignal.emit("UpdatePropertyCache FAIL")
            #         self.logger.error("Update Property Cache failed")

            # if self._sanitize:
            #     Chem.SanitizeMol(mol)
            # if self._kekulize:
            #     Chem.Kekulize(mol)
            self._mol = mol
            self.molChanged.emit()

    def setMol(self, mol):
        self.mol = mol

    # Handling of selections
    selectionChanged = QtCore.Signal(name="selectionChanged")

    def selectAtomAdd(self, atomidx):
        if atomidx not in self._selectedAtoms:
            self._selectedAtoms.append(atomidx)
            self.selectionChanged.emit()

    def selectAtom(self, atomidx):
        self._selectedAtoms = [atomidx]
        self.selectionChanged.emit()

    def unselectAtom(self, atomidx):
        self.selectedAtoms.remove(atomidx)
        self.selectionChanged.emit()

    def clearAtomSelection(self):
        if self._selectedAtoms != []:
            self._selectedAtoms = []
            self.selectionChanged.emit()

    @property
    def selectedAtoms(self):
        return self._selectedAtoms

    @selectedAtoms.setter
    def selectedAtoms(self, atomlist):
        if atomlist != self._selectedAtoms:
            assert isinstance(atomlist, list), "selectedAtoms should be a list of integers"
            assert all(isinstance(item, int) for item in atomlist), "selectedAtoms should be a list of integers"
            self._selectedAtoms = atomlist
            self.selectionChanged.emit()

    def setSelectedAtoms(self, atomlist):
        self.selectedAtoms = atomlist

    drawSettingsChanged = QtCore.Signal(name="drawSettingsChanged")

    @property
    def unsanitizable_background_colour(self):
        return self._unsanitizable_background_colour

    @unsanitizable_background_colour.setter
    def unsanitizable_background_colour(self, colour):
        if colour != self._unsanitizable_background_colour:
            if colour is not None:
                assert validate_rgb(colour)
            self._unsanitizable_background_colour = colour
            self.drawSettingsChanged.emit()

    @property
    def last_selected_highlight_colour(self):
        return self._last_selected_highlight_colour

    @last_selected_highlight_colour.setter
    def last_selected_highlight_colour(self, colour):
        assert validate_rgb(colour)
        if colour != self._last_selected_highlight_colour:
            self._last_selected_highlight_colour = colour
            self.drawSettingsChanged.emit()

    @property
    def selected_highlight_colour(self):
        return self._selected_highlight_colour

    @selected_highlight_colour.setter
    def selected_highlight_colour(self, colour):
        if colour != self._selected_highlight_colour:
            assert validate_rgb(colour)
            self._selected_highlight_colour = colour
            self.drawSettingsChanged.emit()

    # Actions and functions
    @QtCore.Slot()
    def draw(self):
        self.logger.debug("Updating SVG")
        svg = self.getMolSvg()
        self.load(QtCore.QByteArray(svg.encode("utf-8")))

    @QtCore.Slot()
    def sanitize_draw(self):
        # self.computeNewCoords()
        self.sanitizeDrawMol()
        self.draw()

    @QtCore.Slot()
    def changeSanitizeStatus(self, value):
        self.logger.debug(f"changeBorder called with value {value}")
        if value.upper() == "SANITIZABLE":
            self.molecule_sanitizable = True
        else:
            self.molecule_sanitizable = False

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
            self.logger.debug("No Conformers found, computing all 2D coords")
        elif ignoreExisting:
            self.logger.debug("Ignoring existing conformers, computing all 2D coords")
        else:
            assert self._mol.GetNumConformers() == 1
            self.logger.debug("1 Conformer found, computing 2D coords not in found conformer")
            conf = self._mol.GetConformer(0)
            for a in self._mol.GetAtoms():
                pos3d = conf.GetAtomPosition(a.GetIdx())
                if (pos3d.x, pos3d.y) == (0, 0):
                    continue
                prev_coords[a.GetIdx()] = Point2D(pos3d.x, pos3d.y)
        self.logger.debug("Coordmap %s" % prev_coords)
        self.logger.debug("canonOrient %s" % canonOrient)
        rdDepictor.Compute2DCoords(self._mol, coordMap=prev_coords, canonOrient=canonOrient)

    def canon_coords_and_draw(self):
        self.logger.debug("Recalculating coordinates")
        self.computeNewCoords(canonOrient=True, ignoreExisting=True)
        self._drawmol = copy.deepcopy(self._mol)  # Chem.Mol(self._mol.ToBinary())
        self.draw()

    sanitizeSignal = QtCore.Signal(str, name="sanitizeSignal")

    @QtCore.Slot()
    def sanitizeDrawMol(self, kekulize=False, drawkekulize=False):
        self.computeNewCoords()
        # self._drawmol_test = Chem.Mol(self._mol.ToBinary())  # Is this necessary?
        # self._drawmol = Chem.Mol(self._mol.ToBinary())  # Is this necessary?
        self._drawmol_test = copy.deepcopy(self._mol)  # Is this necessary?
        self._drawmol = copy.deepcopy(self._mol)  # Is this necessary?
        try:
            Chem.SanitizeMol(self._drawmol_test)
            self.sanitizeSignal.emit("Sanitizable")
        except Exception as e:
            self.sanitizeSignal.emit("UNSANITIZABLE")
            self.logger.warning("Unsanitizable")
            # try:
            #     self._drawmol.UpdatePropertyCache(strict=False)
            # except Exception as e:
            #     self.sanitizeSignal.emit("UpdatePropertyCache FAIL")
            #     self.logger.error("Update Property Cache failed")
        # # Kekulize
        # if kekulize:
        #     try:
        #         Chem.Kekulize(self._drawmol)
        #     except Exception as e:
        #         self.logger.warning("Unkekulizable")
        # try:
        #     self._drawmol = rdMolDraw2D.PrepareMolForDrawing(self._drawmol, kekulize=drawkekulize)
        # except ValueError:  # <- can happen on a kekulization failure
        #     self._drawmol = rdMolDraw2D.PrepareMolForDrawing(self._drawmol, kekulize=False)
        self._drawmol = rdMolDraw2D.PrepareMolForDrawing(self._drawmol, kekulize=False)

    finishedDrawing = QtCore.Signal(name="finishedDrawing")

    def getMolSvg(self):
        self.drawer = rdMolDraw2D.MolDraw2DSVG(300, 300)
        # TODO, what if self._drawmol doesn't exist?
        if self._drawmol is not None:
            # Chiral tags on R/S
            # chiraltags = Chem.FindMolChiralCenters(self._drawmol)
            opts = self.drawer.drawOptions()
            if self._darkmode:
                rdMolDraw2D.SetDarkMode(opts)
            if (not self.molecule_sanitizable) and self.unsanitizable_background_colour:
                opts.setBackgroundColour(self.unsanitizable_background_colour)
            opts.prepareMolsBeforeDrawing = False
            opts.addStereoAnnotation = True  # Show R/S and E/Z
            # for tag in chiraltags:
            #     idx = tag[0]
            #     opts.atomLabels[idx] = self._drawmol.GetAtomWithIdx(idx).GetSymbol() + ":" + tag[1]

            if len(self._selectedAtoms) > 0:
                colors = {atom_idx: self.selected_highlight_colour for atom_idx in self._selectedAtoms}
                colors[self._selectedAtoms[-1]] = (
                    self.last_selected_highlight_colour
                )  # Color lastly selected an optionally different color
                self.drawer.DrawMolecule(
                    self._drawmol,
                    highlightAtoms=self._selectedAtoms,
                    highlightAtomColors=colors,
                )
            else:
                self.drawer.DrawMolecule(self._drawmol)
        self.drawer.FinishDrawing()
        self.finishedDrawing.emit()  # Signal that drawer has finished
        svg = self.drawer.GetDrawingText().replace("svg:", "")
        return svg


if __name__ == "__main__":
    #    model = SDmodel()
    #    model.loadSDfile('dhfr_3d.sd')
    mol = Chem.MolFromSmiles("CCN(C)c1ccccc1S")
    # rdDepictor.Compute2DCoords(mol)
    myApp = QtWidgets.QApplication(sys.argv)
    molview = MolWidget(mol)
    molview.selectAtom(1)
    molview.selectedAtoms = [1, 2, 3]
    molview.show()
    myApp.exec()
