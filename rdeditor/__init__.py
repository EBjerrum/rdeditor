# from molViewWidget import molViewWidget
# from molEditWidget import MolEditWidget
# from ptable_widget import PTable
try:
    from rdeditor._version import __version__
except ImportError:  # pragma: no cover
    __version__ = "not-installed"

from .rdEditor import MainWindow
