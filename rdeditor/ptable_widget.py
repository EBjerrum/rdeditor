#!/usr/bin/python
# -*- coding: utf-8 -*-
import sys
import logging

from PySide6 import QtGui, QtCore, QtWidgets

from rdeditor.ptable import ptable


class PTable(QtWidgets.QWidget):
    def __init__(self, actionGroup):
        super(PTable, self).__init__()
        self.ptable = ptable
        self.initUI(actionGroup)
        # logging
        self.logger = logging.getLogger()

    def initUI(self, actionGroup):
        grid = QtWidgets.QGridLayout()
        # Create actions dictionary and group dictionary
        # self.atomActionGroup = QtWidgets.QActionGroup(self, exclusive=True)
        self.atomActions = {}
        # for atomname in self.editor.atomtypes.keys(): Gives unsorted list
        for key in self.ptable.keys():
            atomname = self.ptable[key]["Symbol"]
            action = QtGui.QAction(
                "%s" % atomname,
                self,
                statusTip="Set atomtype to %s" % atomname,
                triggered=self.atomtypePush,
                objectName=atomname,
                checkable=True,
            )
            # self.atomActionGroup.addAction(action)
            actionGroup.addAction(action)
            self.atomActions[atomname] = action
            button = QtWidgets.QToolButton()
            button.setDefaultAction(action)
            button.setFocusPolicy(QtCore.Qt.NoFocus)
            button.setMaximumWidth(40)

            if self.ptable[key]["Group"] is not None:
                grid.addWidget(button, self.ptable[key]["Period"], self.ptable[key]["Group"])
            else:
                if key < 72:
                    grid.addWidget(button, 9, key - 54)
                else:
                    grid.addWidget(button, 10, key - 86)
        # Ensure spacing between main table and actinides/lathanides
        grid.addWidget(QtWidgets.QLabel(""), 8, 1)

        self.setLayout(grid)

        self.move(300, 150)
        self.setWindowTitle("Periodic Table")

    atomtypeChanged = QtCore.Signal(str, name="atomtypeChanged")

    def atomtypePush(self):
        sender = self.sender()
        self.atomtypeChanged.emit(sender.objectName())

    # For setting the new atomtype
    def selectAtomtype(self, atomname):
        if atomname in self.atomActions.keys():
            self.atomActions[atomname].setChecked(True)
        else:
            self.debug.error("Unknown atomtype or key error: %s" % atomname)


def main():
    app = QtWidgets.QApplication(sys.argv)
    pt = PTable()
    pt.selectAtomtype("N")
    pt.show()
    sys.exit(app.exec())


if __name__ == "__main__":
    main()
