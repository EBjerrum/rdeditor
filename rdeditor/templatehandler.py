from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Geometry.rdGeometry import Point2D, Point3D
from rdkit.Chem import rdDepictor

from .reaction import ProductTemplate


class TemplateHandler:
    reverse = False
    templates = {
        "benzene": {
            "canvas": "C1=CC=CC=C1",
            "atom": "[998*:1]>>[beginisotope*:1]1=C-C=C-C=C-1",
            "sp3": "[998*:1]-[999*:2]>>[beginisotope*:1]1-[endisotope*:2]=C-C=C-C=1",
            "sp2": "[998*:1]~[999*:2]>>[beginisotope*:1]1~[endisotope*:2]-C=C-C=C-1",
            "aromatic": "[998*:1]~[999*:2]>>[beginisotope*:1]1:[beginisotope*:2]:C:C:C:C:1",
        },
        "cyclohexane": {
            "canvas": "C1CCCCC1",
            "atom": "[998*:1]>>[beginisotope*:1]1CCCCC1",
            "sp3": "[998*:1]-[999*:2]>>[beginisotope*:1]1-[endisotope*:2]-C-C-C-C-1",
            "sp2": "[998*:1]~[999*:2]>>[beginisotope*:1]1~[endisotope*:2]-C-C-C-C-1",
            "aromatic": "[998*:1]~[999*:2]>>[beginisotope*:1]1:[beginisotope*:2]-C-C-C-C-1",
        },
        "cyclopentane": {
            "canvas": "C1CCCC1",
            "atom": "[998*:1]>>[beginisotope*:1]1CCCC1",
            "sp3": "[998*:1]-[999*:2]>>[beginisotope*:1]1-[endisotope*:2]-C-C-C-1",
            "sp2": "[998*:1]~[999*:2]>>[beginisotope*:1]1~[endisotope*:2]-C-C-C-1",
            "aromatic": "[998*:1]~[999*:2]>>[beginisotope*:1]1:[beginisotope*:2]-C-C-C-1",
        },
        "cyclobutane": {
            "canvas": "C1CCC1",
            "atom": "[998*:1]>>[beginisotope*:1]1CCC1",
            "sp3": "[998*:1]-[999*:2]>>[beginisotope*:1]1-[endisotope*:2]-C-C-1",
            "sp2": "[998*:1]~[999*:2]>>[beginisotope*:1]1~[endisotope*:2]-C-C-1",
            "aromatic": "[998*:1]~[999*:2]>>[beginisotope*:1]1:[beginisotope*:2]-C-C-1",
        },
        "cyclopropane": {
            "canvas": "C1CC1",
            "atom": "[998*:1]>>[beginisotope*:1]1CC1",
            "sp3": "[998*:1]-[999*:2]>>[beginisotope*:1]1-[endisotope*:2]-C-1",
            "sp2": "[998*:1]~[999*:2]>>[beginisotope*:1]1~[endisotope*:2]-C-1",
            "aromatic": "[998*:1]~[999*:2]>>[beginisotope*:1]1:[beginisotope*:2]-C-1",
        },
        "carboxylic acid": {
            "canvas": "C(=O)[O]",
            "atom": "[998*:1]>>[beginisotope*:1](=O)[O]",
        },
        # These types of templates need more work, i.e. if an NC bond is clicked, the addition can be non-sanitizable
        # due to the explicit H (or vice versa!)
        # "0-pyrrole": {
        #     "atom": "[998*:1]>>[beginisotope*:1]-[N]1-C=C-C=C-1",
        #     "sp3": "[998*:1]-[999*:2]>>[beginisotope*:1]1-[endisotopeN:2]-C=C-C=1",
        #     "sp2": "[998*:1]=[999*:2]>>[beginisotope*:1]1:[endisotopeN:2]:C:C:C:1",  # Not Kekulized
        #     "aromatic": "[998*:1]~[999*:2]>>[beginisotope*:1]1:[beginisotopeN:2]:c:c:c:1",
        # },
        # "1-pyrrole": {
        #     "atom": "[998*:1]>>[beginisotope*:1]-C1=C-C=C-[NH]1",
        #     Chem.BondType.SINGLE: "[998*:1]-[999*:2]>>[beginisotope*:1]1:[endisotope*:2]:[NH]:C:C:1",  # Not kekulized
        #     "sp2": "[998*:1]=[999*:2]>>[beginisotope*:1]1=[endisotope*:2]-[NH]-C=C-1",
        #     "aromatic": "[998*:1]~[999*:2]>>[beginisotope*:1]1:[beginisotope*:2]:[nH]:c:c:1",
        # },
        # "2-pyrrole": {
        #     "atom": "[998*:1]>>[beginisotope*:1]-C1-C=C-[NH]-C=1",
        #     Chem.BondType.SINGLE: "[998*:1]-[999*:2]>>[beginisotope*:1]1-[endisotope*:2]=C-[NH]-C=1",
        #     "sp2": "[998*:1]=[999*:2]>>[beginisotope*:1]1:[endisotope*:2]:C:[NH]:C:1",  # Not kekulized
        #     "aromatic": "[998*:1]~[999*:2]>>[beginisotope*:1]1:[beginisotope*:2]:c:[nH]:c:1",
        # },
    }

    def __init__(self):
        templates = {
            "benzene": "*[C:1]1=[C:2]C=CC=C1",
            "cyclohexane": "*[C:1]1[C:2]CCCC1",
            "cyclopentane": "*[C:1]1[C:2]CCC1",
            "cyclobutane": "*[C:1]1[C:2]CC1",
            "cyclopropane": "*[C:1]1[C:2]C1",
            "carboxylic acid": "*[C:1](=O)[O]",
        }
        self.productTemplates = {}

        for templateName, template in templates.items():
            self.productTemplates[templateName] = ProductTemplate(template)

    @property
    def templateslabels(self):
        return tuple(self.templates.keys())

    def apply_template_to_atom(self, beginatom: Chem.rdchem.Atom, templatelabel: str) -> Chem.Mol:
        """Apply to an atom"""
        productTemplate = self.productTemplates[templatelabel]

        newmol = productTemplate.apply(beginatom.GetOwningMol(), [beginatom], order=1)

        if newmol:
            return newmol
        else:
            raise RuntimeWarning(f"Applying template returned no molecule!")

    def apply_template_to_bond(self, bond: Chem.rdchem.Bond, templatelabel: str) -> Chem.Mol:
        """Apply to a bond"""
        # Reverse is for future usage for unsymmetric templates. Ctrl-Z and reapply to add the template reversed
        mol = bond.GetOwningMol()
        # mol.UpdatePropertyCache()
        Chem.SetHybridization(mol)

        beginatom = bond.GetEndAtom() if self.reverse else bond.GetBeginAtom()
        endatom = bond.GetBeginAtom() if self.reverse else bond.GetEndAtom()
        self.reverse = not self.reverse

        productTemplate = self.productTemplates[templatelabel]

        newmol = productTemplate.apply(mol, [beginatom, endatom], order=2)

        if newmol:
            return newmol
        else:
            raise RuntimeWarning(f"Applying template returned no molecule!")

    def apply_template_to_canvas(self, mol: Chem.Mol, point: Point2D, templatelabel: str) -> Chem.Mol:
        """Apply to canvas"""
        productTemplate = self.productTemplates[templatelabel]

        newmol = productTemplate.apply(mol, [], order=0)

        if newmol:
            return newmol
        else:
            raise RuntimeWarning(f"Applying template returned no molecule!")

    def react_and_keep_fragments(self, mol, rxn):
        """RDKit only returns the fragment of Mol object that is reacted,
        hence this function to keep all other fragments."""
        # Split the molecule into fragments
        fragments = list(Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=False))

        # Perform the reaction on the reacting fragment
        for i, reacting_fragment in enumerate(fragments):
            products = rxn.RunReactants((reacting_fragment,))
            if products:
                fragments.pop(i)
                result = products[0][0]
                for frag in fragments:
                    result = Chem.CombineMols(result, frag)

                return result
        else:
            return None
