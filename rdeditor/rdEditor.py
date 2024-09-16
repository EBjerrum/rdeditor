#!/usr/bin/env python

# Import required modules
import sys
import os
import copy

from PySide6.QtWidgets import QMenu, QApplication, QStatusBar, QMessageBox, QFileDialog
from PySide6.QtCore import QSettings
from PySide6 import QtCore, QtGui, QtWidgets
from PySide6.QtCore import QUrl
from PySide6.QtGui import QDesktopServices, QIcon, QAction, QKeySequence

import qdarktheme

# Import model
from . import __version__
from .molEditWidget import MolEditWidget
from .ptable_widget import PTable

from rdkit import Chem
import qdarktheme


# The main window class
class MainWindow(QtWidgets.QMainWindow):
    # Constructor function
    def __init__(self, fileName=None, loglevel="WARNING"):
        super(MainWindow, self).__init__()
        self.pixmappath = os.path.abspath(os.path.dirname(__file__)) + "/pixmaps/"
        QtGui.QIcon.setThemeSearchPaths(
            # QtGui.QIcon.themeSearchPaths() +
            [os.path.abspath(os.path.dirname(__file__)) + "/icon_themes/"]
        )
        self.loglevels = ["Critical", "Error", "Warning", "Info", "Debug", "Notset"]
        self.editor = MolEditWidget()
        self.chemEntityActionGroup = QtGui.QActionGroup(self, exclusive=True)
        self.ptable = PTable(self.chemEntityActionGroup)
        self._fileName = None
        self.initGUI(fileName=fileName)
        self.applySettings()
        self.ptable.atomtypeChanged.connect(self.setAtomTypeName)
        # self.singleBondAction.trigger()
        self.ptable.atomActions["C"].trigger()

    # Properties
    @property
    def fileName(self):
        return self._fileName

    @fileName.setter
    def fileName(self, filename):
        if filename != self._fileName:
            self._fileName = filename
            self.setWindowTitle(str(filename))

    def initGUI(self, fileName=None):
        self.setWindowTitle("rdEditor")
        self.setWindowIcon(QIcon.fromTheme("appicon"))
        self.setGeometry(100, 100, 200, 150)

        self.center = self.editor
        self.center.setFixedSize(650, 650)
        self.setCentralWidget(self.center)
        self.fileName = fileName

        self.filters = "MOL Files (*.mol *.mol);;SMILES Files (*.smi *.smi);;Any File (*)"

        self.SetupComponents()

        self.infobar = QtWidgets.QLabel("")
        self.myStatusBar.addPermanentWidget(self.infobar, 0)

        if self.fileName is not None:
            self.editor.logger.info("Loading molecule from %s" % self.fileName)
            self.loadFile()

        self.editor.sanitizeSignal.connect(self.infobar.setText)

        self.show()

    def getAllActionsInMenu(self, qmenu: QMenu):
        all_actions = []

        # Iterate through actions in the current menu
        for action in qmenu.actions():
            if isinstance(action, QAction):
                if action.icon():
                    all_actions.append(action)
            elif isinstance(action, QMenu):  # If the action is a submenu, recursively get its actions
                all_actions.extend(self.getAllActionsInMenu(action))

        return all_actions

    def getAllIconActions(self, qapp: QApplication):
        all_actions = []

        # Iterate through all top-level widgets in the application
        for widget in qapp.topLevelWidgets():
            # Find all menus in the widget
            menus = widget.findChildren(QMenu)
            for menu in menus:
                # Recursively get all actions from each menu
                all_actions.extend(self.getAllActionsInMenu(menu))

        return all_actions

    def resetActionIcons(self):
        actions_with_icons = list(set(self.getAllIconActions(QApplication)))
        for action in actions_with_icons:
            icon_name = action.icon().name()
            self.editor.logger.debug(f"reset icon {icon_name}")
            action.setIcon(QIcon.fromTheme(icon_name))

    def applySettings(self):
        self.settings = QSettings("Cheminformania.com", "rdEditor")
        theme_name = self.settings.value("theme_name", "Fusion")

        self.applyTheme(theme_name)
        self.themeActions[theme_name].setChecked(True)

        loglevel = self.settings.value("loglevel", "Error")

        action = self.loglevelactions.get(loglevel, None)
        if action:
            action.trigger()

        sanitize_on_cleanup = self.settings.value("sanitize_on_cleanup", True, type=bool)
        self.editor.sanitize_on_cleanup = sanitize_on_cleanup
        self.cleanupSettingActions["sanitize_on_cleanup"].setChecked(sanitize_on_cleanup)

        kekulize_on_cleanup = self.settings.value("kekulize_on_cleanup", True, type=bool)
        self.editor.kekulize_on_cleanup = kekulize_on_cleanup
        self.cleanupSettingActions["kekulize_on_cleanup"].setChecked(kekulize_on_cleanup)

    # Function to setup status bar, central widget, menu bar, tool bar
    def SetupComponents(self):
        self.myStatusBar = QStatusBar()
        #        self.molcounter = QLabel("-/-")
        #        self.myStatusBar.addPermanentWidget(self.molcounter, 0)
        self.setStatusBar(self.myStatusBar)
        self.myStatusBar.showMessage("Ready", 10000)

        self.CreateActions()
        self.CreateMenus()
        self.CreateToolBars()

    # Actual menu bar item creation
    def CreateMenus(self):
        self.fileMenu = self.menuBar().addMenu("&File")
        # self.edit_menu = self.menuBar().addMenu("&Edit")

        self.toolMenu = self.menuBar().addMenu("&Tools")
        self.atomtypeMenu = self.menuBar().addMenu("&AtomTypes")
        self.bondtypeMenu = self.menuBar().addMenu("&BondTypes")
        self.templateMenu = self.menuBar().addMenu("Tem&plates")
        self.settingsMenu = self.menuBar().addMenu("&Settings")
        self.helpMenu = self.menuBar().addMenu("&Help")

        self.fileMenu.addAction(self.openAction)
        self.fileMenu.addAction(self.saveAction)
        self.fileMenu.addAction(self.saveAsAction)
        self.fileMenu.addSeparator()
        self.fileMenu.addAction(self.copyAction)
        self.fileMenu.addAction(self.pasteAction)
        self.fileMenu.addSeparator()
        self.fileMenu.addAction(self.exitAction)

        self.toolMenu.addAction(self.selectAction)
        self.toolMenu.addAction(self.addAction)
        # self.toolMenu.addAction(self.addBondAction)
        # self.toolMenu.addAction(self.replaceAction)
        self.toolMenu.addAction(
            self.rsAction
        )  # TODO, R/S and E/Z could be changed for a single action? it really depends if an atom or a bond is clicked!
        self.toolMenu.addAction(self.ezAction)
        self.toolMenu.addAction(self.increaseChargeAction)
        self.toolMenu.addAction(self.decreaseChargeAction)
        self.toolMenu.addAction(self.numberAtom)

        self.toolMenu.addSeparator()
        self.toolMenu.addAction(self.cleanCoordinatesAction)
        self.toolMenu.addAction(self.cleanupMolAction)
        self.toolMenu.addSeparator()
        self.toolMenu.addAction(self.undoAction)
        self.toolMenu.addSeparator()
        self.toolMenu.addAction(self.removeAction)
        self.toolMenu.addAction(self.clearCanvasAction)

        # Atomtype menu
        for action in self.atomActions:
            self.atomtypeMenu.addAction(action)
        self.specialatommenu = self.atomtypeMenu.addMenu("All Atoms")
        for atomnumber in self.ptable.ptable.keys():
            atomname = self.ptable.ptable[atomnumber]["Symbol"]
            self.specialatommenu.addAction(self.ptable.atomActions[atomname])

        # Bondtype Menu
        self.bondtypeMenu.addAction(self.singleBondAction)
        self.bondtypeMenu.addAction(self.doubleBondAction)
        self.bondtypeMenu.addAction(self.tripleBondAction)
        self.bondtypeMenu.addSeparator()
        # Bondtype Special types
        self.specialbondMenu = self.bondtypeMenu.addMenu("Special Bonds")
        for key in self.bondActions.keys():
            self.specialbondMenu.addAction(self.bondActions[key])

        # Templates menu
        for key in self.templateActions.keys():
            self.templateMenu.addAction(self.templateActions[key])

        # Settings menu
        self.themeMenu = self.settingsMenu.addMenu("Theme")
        self.populateThemeActions(self.themeMenu)
        self.loglevelMenu = self.settingsMenu.addMenu("Logging Level")
        for loglevel in self.loglevels:
            self.loglevelMenu.addAction(self.loglevelactions[loglevel])
        self.cleanupMenu = self.settingsMenu.addMenu("Cleanup")
        for key, action in self.cleanupSettingActions.items():
            self.cleanupMenu.addAction(action)

        # Help menu
        self.helpMenu.addAction(self.aboutAction)
        self.helpMenu.addSeparator()
        self.helpMenu.addAction(self.openChemRxiv)
        self.helpMenu.addAction(self.openRepository)
        self.helpMenu.addSeparator()
        self.helpMenu.addAction(self.aboutQtAction)

        # actionListAction = QAction(
        #     "List Actions", self, triggered=lambda: print(set(self.get_all_icon_actions_in_application(QApplication)))
        # )
        # self.helpMenu.addAction(actionListAction)

        # Debug level sub menu

    def populateThemeActions(self, menu: QMenu):
        stylelist = QtWidgets.QStyleFactory.keys() + ["Qdt light", "Qdt dark"]
        self.themeActionGroup = QtGui.QActionGroup(self, exclusive=True)
        self.themeActions = {}
        for style_name in stylelist:
            action = QAction(
                style_name,
                self,
                objectName=style_name,
                triggered=self.setTheme,
                checkable=True,
            )
            self.themeActionGroup.addAction(action)
            self.themeActions[style_name] = action
            menu.addAction(action)

    def CreateToolBars(self):
        self.mainToolBar = self.addToolBar("Main")
        # Main action bar
        self.mainToolBar.addAction(self.openAction)
        self.mainToolBar.addAction(self.saveAction)
        self.mainToolBar.addAction(self.saveAsAction)
        self.mainToolBar.addSeparator()
        self.mainToolBar.addAction(self.selectAction)
        self.mainToolBar.addAction(self.addAction)
        # self.mainToolBar.addAction(self.addBondAction)
        # self.mainToolBar.addAction(self.replaceAction)
        self.mainToolBar.addAction(self.rsAction)
        self.mainToolBar.addAction(self.ezAction)
        self.mainToolBar.addAction(self.increaseChargeAction)
        self.mainToolBar.addAction(self.decreaseChargeAction)
        self.mainToolBar.addAction(self.numberAtom)
        self.mainToolBar.addSeparator()
        self.mainToolBar.addAction(self.cleanCoordinatesAction)
        self.mainToolBar.addAction(self.cleanupMolAction)

        self.mainToolBar.addSeparator()
        self.mainToolBar.addAction(self.removeAction)
        self.mainToolBar.addAction(self.clearCanvasAction)
        # Bond types TODO are they necessary as can be toggled??
        self.mainToolBar.addSeparator()
        self.mainToolBar.addAction(self.undoAction)
        # Side Toolbar
        self.sideToolBar = QtWidgets.QToolBar(self)
        self.addToolBar(QtCore.Qt.LeftToolBarArea, self.sideToolBar)
        self.sideToolBar.addAction(self.singleBondAction)
        self.sideToolBar.addAction(self.doubleBondAction)
        self.sideToolBar.addAction(self.tripleBondAction)
        self.sideToolBar.addSeparator()
        self.sideToolBar.addAction(self.templateActions["benzene"])
        self.sideToolBar.addAction(self.templateActions["cyclohexane"])
        self.sideToolBar.addSeparator()
        self.sideToolBar.addAction(self.ptable.atomActions["R"])
        self.sideToolBar.addSeparator()
        for action in self.atomActions:
            self.sideToolBar.addAction(action)
        self.sideToolBar.addAction(self.openPtableAction)

    def loadSmilesFile(self, filename):
        self.fileName = filename
        with open(self.fileName, "r") as file:
            lines = file.readlines()
            if len(lines) > 1:
                self.editor.logger.warning("The SMILES file contains more than one line.")
                self.statusBar().showMessage("The SMILES file contains more than one line.")
                return None
            smiles = lines[0].strip()
            mol = Chem.MolFromSmiles(smiles)
            self.editor.mol = mol
            self.statusBar().showMessage(f"SMILES file {filename} opened")

    def loadMolFile(self, filename):
        self.fileName = filename
        mol = Chem.MolFromMolFile(str(self.fileName), sanitize=False, strictParsing=False)
        self.editor.mol = mol
        self.statusBar().showMessage(f"Mol file {filename} opened")

    def openFile(self):
        self.fileName, _ = QFileDialog.getOpenFileName(self, caption="Open file", filter=self.filters)
        return self.loadFile()

    def loadFile(self):
        if not self.fileName:
            self.editor.logger.warning("No file selected.")
            self.statusBar().showMessage("No file selected.")
            return
        if self.fileName.lower().endswith(".mol"):
            self.loadMolFile(self.fileName)
        elif self.fileName.lower().endswith(".smi"):
            self.loadSmilesFile(self.fileName)
        else:
            self.editor.logger.warning("Unknown file format. Assuming file as .mol format.")
            self.statusBar().showMessage("Unknown file format. Assuming file as .mol format.")
            self.loadMolFile(self.fileName)
            self.fileName += ".mol"

    def saveFile(self):
        if self.fileName is not None:
            Chem.MolToMolFile(self.editor.mol, str(self.fileName))
        else:
            self.saveAsFile()

    def saveAsFile(self):
        self.fileName, self.filterName = QFileDialog.getSaveFileName(self, filter=self.filters)
        if self.fileName != "":
            if self.filterName == "MOL Files (*.mol *.mol)":
                if not self.fileName.lower().endswith(".mol"):
                    self.fileName = self.fileName + ".mol"
                Chem.MolToMolFile(self.editor.mol, str(self.fileName))
                self.statusBar().showMessage("File saved as MolFile", 2000)
            elif self.filterName == "SMILES Files (*.smi *.smi)":
                if not self.fileName.lower().endswith(".smi"):
                    self.fileName = self.fileName + ".smi"
                smiles = Chem.MolToSmiles(self.editor.mol)
                with open(self.fileName, "w") as file:
                    file.write(smiles + "\n")
                self.statusBar().showMessage("File saved as SMILES", 2000)
            else:
                self.statusBar().showMessage("Invalid file format", 2000)

    def copy(self):
        selected_text = Chem.MolToSmiles(self.editor.mol, isomericSmiles=True)
        clipboard = QApplication.clipboard()
        clipboard.setText(selected_text)

    def paste(self):
        clipboard = QApplication.clipboard()
        text = clipboard.text()
        mol = Chem.MolFromSmiles(text, sanitize=False)
        if mol:
            try:
                Chem.SanitizeMol(copy.deepcopy(mol))  # ).ToBinary()))
            except Exception as e:
                self.editor.logger.warning(f"Pasted SMILES is not sanitizable: {e}")

            self.editor.assign_stereo_atoms(mol)
            Chem.rdmolops.SetBondStereoFromDirections(mol)

            self.editor.mol = mol
        else:
            self.editor.logger.warning(f"Failed to parse the content of the clipboard as a SMILES: {repr(text)}")

    def clearCanvas(self):
        self.editor.clearAtomSelection()
        self.editor.mol = None
        self.fileName = None
        self.statusBar().showMessage("Canvas Cleared")

    def closeEvent(self, event):
        self.editor.logger.debug("closeEvent triggered")
        self.exitFile()
        event.ignore()

    def exitFile(self):
        response = self.msgApp("Confirmation", "This will quit the application. Do you want to Continue?")
        if response == "Y":
            self.ptable.close()
            exit(0)
        else:
            self.editor.logger.debug("Abort closing")

    # Function to show Diaglog box with provided Title and Message
    def msgApp(self, title, msg):
        userInfo = QMessageBox.question(self, title, msg)
        if userInfo == QMessageBox.Yes:
            return "Y"
        if userInfo == QMessageBox.No:
            return "N"
        self.close()

    def aboutHelp(self):
        QMessageBox.about(
            self,
            "About Simple Molecule Editor",
            f"""A Simple Molecule Editor where you can edit molecules\n
Based on RDKit! http://www.rdkit.org/\n
Some icons from http://icons8.com\n
Source code: https://github.com/EBjerrum/rdeditor\n
Version: {__version__}
            """,
        )

    def setAction(self):
        sender = self.sender()
        self.editor.setAction(sender.objectName())
        self.myStatusBar.showMessage("Action %s selected" % sender.objectName())

    # TODO, the various setTypes could be unified, as the editor now understands a single chementity
    def setRingType(self):
        sender = self.sender()
        self.editor.setChemEntity(sender.objectName())
        self.myStatusBar.showMessage("Ringtype %s selected" % sender.objectName())

    def setBondType(self):
        sender = self.sender()
        self.editor.setChemEntity(sender.objectName())
        self.myStatusBar.showMessage("Bondtype %s selected" % sender.objectName())

    def setAtomType(self):
        sender = self.sender()
        self.editor.setChemEntity(sender.objectName())
        # self.editor.setRingType(None)
        self.myStatusBar.showMessage("Atomtype %s selected" % sender.objectName())

    def setAtomTypeName(self, atomname):
        self.editor.setChemEntity(str(atomname))
        self.myStatusBar.showMessage("Atomtype %s selected" % atomname)

    def openPtable(self):
        self.ptable.show()

    def setLogLevel(self):
        loglevel = self.sender().objectName().split(":")[-1]  # .upper()
        self.editor.logger.setLevel(loglevel.upper())
        self.editor.logger.log(self.editor.logger.getEffectiveLevel(), f"loglevel set to {loglevel}")
        self.settings.setValue("loglevel", loglevel)
        self.settings.sync()

    def setTheme(self):
        sender = self.sender()
        theme_name = sender.objectName()
        self.myStatusBar.showMessage(f"Setting theme or style to {theme_name}")
        self.applyTheme(theme_name)
        self.settings.setValue("theme_name", theme_name)
        self.settings.sync()

    def is_dark_mode(self):
        """Hack to detect if we have a dark mode running"""
        app = QApplication.instance()
        palette = app.palette()
        # Get the color of the window background
        background_color = palette.color(QtGui.QPalette.Window)
        # Calculate the luminance (brightness) of the color
        luminance = (
            0.299 * background_color.red() + 0.587 * background_color.green() + 0.114 * background_color.blue()
        ) / 255
        # If the luminance is below a certain threshold, it's considered dark mode
        return luminance < 0.5

    def applyTheme(self, theme_name):
        if "dark" in theme_name:
            self.set_dark()
        elif "light" in theme_name:
            self.set_light()
        elif self.is_dark_mode():
            self.set_dark()
        else:
            self.set_light()

        app = QApplication.instance()
        app.setStyleSheet("")  # resets style
        if theme_name in QtWidgets.QStyleFactory.keys():
            app.setStyle(theme_name)
        else:
            if theme_name == "Qdt light":
                qdarktheme.setup_theme("light")
            elif theme_name == "Qdt dark":
                qdarktheme.setup_theme("dark")

        self.resetActionIcons()

    def set_light(self):
        QIcon.setThemeName("light")
        self.editor.darkmode = False
        self.editor.logger.info("Resetting theme for light theme")

    def set_dark(self):
        QIcon.setThemeName("dark")
        self.editor.darkmode = True
        self.editor.logger.info("Resetting theme for dark theme")

    def openUrl(self):
        url = self.sender().data()
        QDesktopServices.openUrl(QUrl(url))

    def set_setting(self):
        action = self.sender()
        if isinstance(action, QAction):
            setting_name = action.objectName()
            if hasattr(self.editor, setting_name):
                if getattr(self.editor, setting_name) != action.isChecked():
                    setattr(self.editor, setting_name, action.isChecked())
                    self.editor.logger.error(f"Changed editor setting {setting_name} to {action.isChecked()}")
                    self.settings.setValue(setting_name, action.isChecked())
                    self.settings.sync()
            else:
                self.editor.logger.error(f"Error, could not find setting, {setting_name}, on editor object!")

    # Function to create actions for menus and toolbars
    def CreateActions(self):
        self.openAction = QAction(
            QIcon.fromTheme("open"),
            "O&pen",
            self,
            shortcut=QKeySequence.Open,
            statusTip="Open an existing file",
            triggered=self.openFile,
        )

        self.saveAction = QAction(
            QIcon.fromTheme("icons8-Save"),
            "S&ave",
            self,
            shortcut=QKeySequence.Save,
            statusTip="Save file",
            triggered=self.saveFile,
        )

        self.saveAsAction = QAction(
            QIcon.fromTheme("icons8-Save as"),
            "Save As",
            self,
            shortcut=QKeySequence.SaveAs,
            statusTip="Save file as ..",
            triggered=self.saveAsFile,
        )

        self.exitAction = QAction(
            QIcon.fromTheme("icons8-Shutdown"),
            "E&xit",
            self,
            shortcut="Ctrl+Q",
            statusTip="Exit the Application",
            triggered=self.exitFile,
        )

        self.aboutAction = QAction(
            QIcon.fromTheme("about"),
            "A&bout",
            self,
            statusTip="Displays info about text editor",
            triggered=self.aboutHelp,
        )

        self.aboutQtAction = QAction(
            "About &Qt",
            self,
            statusTip="Show the Qt library's About box",
            triggered=QApplication.aboutQt,
        )

        self.openPtableAction = QAction(
            QIcon.fromTheme("ptable"),
            "O&pen Periodic Table",
            self,
            shortcut=QKeySequence.Open,
            statusTip="Open the periodic table for atom type selection",
            triggered=self.openPtable,
        )

        # Copy-Paste actions
        self.copyAction = QAction(
            QIcon.fromTheme("icons8-copy-96"),
            "Copy SMILES",
            self,
            shortcut=QKeySequence.Copy,
            statusTip="Copy the current molecule as a SMILES string",
            triggered=self.copy,
        )

        self.pasteAction = QAction(
            QIcon.fromTheme("icons8-paste-100"),
            "Paste SMILES",
            self,
            shortcut=QKeySequence.Paste,
            statusTip="Paste the clipboard and parse assuming it is a SMILES string",
            triggered=self.paste,
        )

        # Edit actions
        self.actionActionGroup = QtGui.QActionGroup(self, exclusive=True)
        self.selectAction = QAction(
            QIcon.fromTheme("icons8-Cursor"),
            "Se&lect",
            self,
            shortcut="Ctrl+L",
            statusTip="Select Atoms",
            triggered=self.setAction,
            objectName="Select",
            checkable=True,
        )
        self.actionActionGroup.addAction(self.selectAction)

        self.addAction = QAction(
            QIcon.fromTheme("icons8-Edit"),
            "&Add",
            self,
            shortcut="Ctrl+A",
            statusTip="Add Atoms",
            triggered=self.setAction,
            objectName="Add",
            checkable=True,
        )
        self.actionActionGroup.addAction(self.addAction)

        self.addBondAction = QAction(
            QIcon.fromTheme("icons8-Pinch"),
            "Add &Bond",
            self,
            shortcut="Ctrl+B",
            statusTip="Add Bond",
            triggered=self.setAction,
            objectName="Add Bond",
            checkable=True,
        )
        self.actionActionGroup.addAction(self.addBondAction)

        self.replaceAction = QAction(
            QIcon.fromTheme("icons8-Replace Atom"),
            "&Replace",
            self,
            shortcut="Ctrl+R",
            statusTip="Replace Atom/Bond",
            triggered=self.setAction,
            objectName="Replace",
            checkable=True,
        )
        self.actionActionGroup.addAction(self.replaceAction)

        self.rsAction = QAction(
            QIcon.fromTheme("Change_R_S"),
            "To&ggle R/S",
            self,
            shortcut="Ctrl+G",
            statusTip="Toggle Stereo Chemistry",
            triggered=self.setAction,
            objectName="RStoggle",
            checkable=True,
        )
        self.actionActionGroup.addAction(self.rsAction)

        self.ezAction = QAction(
            QIcon.fromTheme("Change_E_Z"),
            "Toggle &E/Z",
            self,
            shortcut="Ctrl+E",
            statusTip="Toggle Bond Stereo Chemistry",
            triggered=self.setAction,
            objectName="EZtoggle",
            checkable=True,
        )
        self.actionActionGroup.addAction(self.ezAction)

        self.removeAction = QAction(
            QIcon.fromTheme("icons8-Cancel"),
            "D&elete",
            self,
            shortcut="Ctrl+D",
            statusTip="Delete Atom or Bond",
            triggered=self.setAction,
            objectName="Remove",
            checkable=True,
        )
        self.actionActionGroup.addAction(self.removeAction)

        self.increaseChargeAction = QAction(
            QIcon.fromTheme("icons8-Increase Font"),
            "I&ncrease Charge",
            self,
            shortcut="Ctrl++",
            statusTip="Increase Atom Charge",
            triggered=self.setAction,
            objectName="Increase Charge",
            checkable=True,
        )
        self.actionActionGroup.addAction(self.increaseChargeAction)

        self.decreaseChargeAction = QAction(
            QIcon.fromTheme("icons8-Decrease Font"),
            "D&ecrease Charge",
            self,
            shortcut="Ctrl+-",
            statusTip="Decrease Atom Charge",
            triggered=self.setAction,
            objectName="Decrease Charge",
            checkable=True,
        )
        self.actionActionGroup.addAction(self.decreaseChargeAction)

        self.numberAtom = QAction(
            QIcon.fromTheme("atommapnumber"),
            "Set atommap or R-group number",
            self,
            statusTip="Set atommap or R-group number",
            triggered=self.setAction,
            objectName="Number Atom",
            checkable=True,
        )
        self.actionActionGroup.addAction(self.numberAtom)
        self.addAction.setChecked(True)

        # BondTypeActions
        self.bondtypeActionGroup = QtGui.QActionGroup(self, exclusive=True)

        self.singleBondAction = QAction(
            QIcon.fromTheme("icons8-Single"),
            "S&ingle Bond",
            self,
            shortcut="Ctrl+1",
            statusTip="Set bondtype to SINGLE",
            triggered=self.setBondType,
            objectName="SINGLE",
            checkable=True,
        )
        self.chemEntityActionGroup.addAction(self.singleBondAction)

        self.doubleBondAction = QAction(
            QIcon.fromTheme("icons8-Double"),
            "Double Bond",
            self,
            shortcut="Ctrl+2",
            statusTip="Set bondtype to DOUBLE",
            triggered=self.setBondType,
            objectName="DOUBLE",
            checkable=True,
        )
        self.chemEntityActionGroup.addAction(self.doubleBondAction)

        self.tripleBondAction = QAction(
            QIcon.fromTheme("icons8-Triple"),
            "Triple Bond",
            self,
            shortcut="Ctrl+3",
            statusTip="Set bondtype to TRIPLE",
            triggered=self.setBondType,
            objectName="TRIPLE",
            checkable=True,
        )
        self.chemEntityActionGroup.addAction(self.tripleBondAction)

        # Build dictionary of ALL available bondtypes in RDKit
        self.bondActions = {}
        for key in self.editor.bondtypes.keys():
            action = QAction(
                "%s" % key,
                self,
                statusTip="Set bondtype to %s" % key,
                triggered=self.setBondType,
                objectName=key,
                checkable=True,
            )
            self.chemEntityActionGroup.addAction(action)
            self.bondActions[key] = action
        # Replace defined actions

        self.bondActions["SINGLE"] = self.singleBondAction
        self.bondActions["DOUBLE"] = self.doubleBondAction
        self.bondActions["TRIPLE"] = self.tripleBondAction

        # self.singleBondAction.setChecked(True)

        # Template Actions
        self.templateActions = {}

        # TODO can we add these automatically, i.e. if theres a similar named icon available?
        self.templateActions["benzene"] = QAction(
            QIcon.fromTheme("benzene"),
            "Benzene Ring",
            self,
            shortcut="Ctrl+4",
            statusTip="Draw Benzene",
            triggered=self.setRingType,
            objectName="benzene",
            checkable=True,
        )

        self.templateActions["cyclohexane"] = QAction(
            QIcon.fromTheme("cyclohexane"),
            "Cyclohexane",
            self,
            shortcut="Ctrl+5",
            statusTip="Draw Cyclohexane",
            triggered=self.setRingType,
            objectName="cyclohexane",
            checkable=True,
        )

        for key in self.editor.available_rings:
            if key not in self.templateActions:
                action = QAction(
                    key,
                    self,
                    statusTip=f"Set template to {key}",
                    triggered=self.setRingType,
                    objectName=key,
                    checkable=True,
                )
                self.templateActions[key] = action

        for action in self.templateActions.values():
            self.chemEntityActionGroup.addAction(action)

        # Misc Actions
        self.undoAction = QAction(
            QIcon.fromTheme("prev"),
            "U&ndo",
            self,
            shortcut="Ctrl+Z",
            statusTip="Undo/Redo changes to molecule Ctrl+Z",
            triggered=self.editor.undo,
            objectName="undo",
        )

        self.clearCanvasAction = QAction(
            QIcon.fromTheme("icons8-Trash"),
            "C&lear Canvas",
            self,
            shortcut="Ctrl+X",
            statusTip="Clear Canvas (no warning)",
            triggered=self.clearCanvas,
            objectName="Clear Canvas",
        )

        self.cleanCoordinatesAction = QAction(
            QIcon.fromTheme("RecalcCoord"),
            "Recalculate coordinates &F",
            self,
            shortcut="Ctrl+F",
            statusTip="Re-calculates coordinates and redraw",
            triggered=self.editor.canon_coords_and_draw,
            objectName="Recalculate Coordinates",
        )

        self.cleanupMolAction = QAction(
            QIcon.fromTheme("CleanupChem"),
            "Cleanup Chemistry",
            self,
            # shortcut="Ctrl+F",
            statusTip="Sanitizes and Kekulizes molecule according to settings",
            triggered=self.editor.cleanup_mol,
            objectName="Cleanup Mol",
        )

        self.cleanupSettingActions = {}

        self.sanitizeSettingAction = QAction(
            # QIcon.fromTheme("icons8-Broom"),
            "Sanitize molecule on cleanup",
            self,
            # shortcut="Ctrl+F",
            statusTip="Perform Sanitization during chemistry cleanup",
            triggered=self.set_setting,
            checkable=True,
            # checked=self.editor.sanitize_on_cleanup,
            objectName="sanitize_on_cleanup",
        )
        self.cleanupSettingActions["sanitize_on_cleanup"] = self.sanitizeSettingAction

        self.kekulizeSettingAction = QAction(
            # QIcon.fromTheme("icons8-Broom"),
            "Kekulize molecule on cleanup",
            self,
            # shortcut="Ctrl+F",
            statusTip="Perform Kekulization after chemistry cleanup",
            triggered=self.set_setting,
            checkable=True,
            # checked=self.editor.kekulize_on_cleanup,
            objectName="kekulize_on_cleanup",
        )
        self.cleanupSettingActions["kekulize_on_cleanup"] = self.kekulizeSettingAction

        # Atom Actions in actiongroup, reuse from ptable widget
        self.atomActions = []
        for atomname in ["H", "B", "C", "N", "O", "F", "P", "S", "Cl", "Br", "I"]:
            action = self.ptable.atomActions[atomname]
            self.atomActions.append(action)

        self.loglevelactions = {}
        self.loglevelActionGroup = QtGui.QActionGroup(self, exclusive=True)
        for key in self.loglevels:
            self.loglevelactions[key] = QAction(
                key,
                self,
                statusTip="Set logging level to %s" % key,
                triggered=self.setLogLevel,
                objectName="loglevel:%s" % key,
                checkable=True,
            )
            self.loglevelActionGroup.addAction(self.loglevelactions[key])

        self.openChemRxiv = QAction(
            QIcon.fromTheme("icons8-Exit"),
            "ChemRxiv Preprint",
            self,
            # shortcut="Ctrl+F",
            statusTip="Opens the ChemRxiv preprint",
            triggered=self.openUrl,
            data="https://doi.org/10.26434/chemrxiv-2024-jfhmw",
        )

        self.openRepository = QAction(
            QIcon.fromTheme("icons8-Exit"),
            "GitHub repository",
            self,
            # shortcut="Ctrl+F",
            statusTip="Opens the GitHub repository",
            triggered=self.openUrl,
            data="https://github.com/EBjerrum/rdeditor",
        )


def launch(loglevel="WARNING"):
    "Function that launches the mainWindow Application"
    # Exception Handling
    try:
        myApp = QApplication(sys.argv)
        if len(sys.argv) > 1:
            mainWindow = MainWindow(fileName=sys.argv[1], loglevel=loglevel)
        else:
            mainWindow = MainWindow(loglevel=loglevel)
        myApp.exec()
        sys.exit(0)
    except NameError:
        print("Name Error:", sys.exc_info()[1])
    except SystemExit:
        print("Closing Window...")
    except Exception:
        print(sys.exc_info()[1])


if __name__ == "__main__":
    launch(loglevel="DEBUG")
