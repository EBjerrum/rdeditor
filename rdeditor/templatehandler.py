from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Geometry.rdGeometry import Point2D, Point3D
from rdkit.Chem import rdDepictor


class TemplateHandler:
    reverse = False
    templates = {
        "benzene": {
            "canvas": "c1ccccc1",
            "atom": "[998*:1]>>[beginisotope*:1]-C1=C-C=C-C=C-1",
            Chem.BondType.SINGLE: "[998*:1]-[999*:2]>>[beginisotope*:1]1-[endisotope*:2]=C-C=C-C=1",
            Chem.BondType.DOUBLE: "[998*:1]=[999*:2]>>[beginisotope*:1]1=[endisotope*:2]-C=C-C=C-1",
            Chem.BondType.AROMATIC: "[998*:1]~[999*:2]>>[beginisotope*:1]1:[beginisotope*:2]:C:C:C:C:1",
        },
        "cyclohexane": {
            "canvas": "C1CCCCC1",
            "atom": "[998*:1]>>[beginisotope*:1]-C1CCCCC1",
            Chem.BondType.SINGLE: "[998*:1]-[999*:2]>>[beginisotope*:1]1-[endisotope*:2]-C-C-C-C-1",
            Chem.BondType.DOUBLE: "[998*:1]=[999*:2]>>[beginisotope*:1]1=[endisotope*:2]-C-C-C-C-1",
            Chem.BondType.AROMATIC: "[998*:1]~[999*:2]>>[beginisotope*:1]1:[beginisotope*:2]-C-C-C-C-1",
        },
        "cyclopentane": {
            "canvas": "C1CCCC1",
            "atom": "[998*:1]>>[beginisotope*:1]-C1CCCC1",
            Chem.BondType.SINGLE: "[998*:1]-[999*:2]>>[beginisotope*:1]1-[endisotope*:2]-C-C-C-1",
            Chem.BondType.DOUBLE: "[998*:1]=[999*:2]>>[beginisotope*:1]1=[endisotope*:2]-C-C-C-1",
            Chem.BondType.AROMATIC: "[998*:1]~[999*:2]>>[beginisotope*:1]1:[beginisotope*:2]-C-C-C-1",
        },
        "cyclobutane": {
            "canvas": "C1CCC1",
            "atom": "[998*:1]>>[beginisotope*:1]-C1CCC1",
            Chem.BondType.SINGLE: "[998*:1]-[999*:2]>>[beginisotope*:1]1-[endisotope*:2]-C-C-1",
            Chem.BondType.DOUBLE: "[998*:1]=[999*:2]>>[beginisotope*:1]1=[endisotope*:2]-C-C-1",
            Chem.BondType.AROMATIC: "[998*:1]~[999*:2]>>[beginisotope*:1]1:[beginisotope*:2]-C-C-1",
        },
        "cyclopropane": {
            "canvas": "C1CC1",
            "atom": "[998*:1]>>[beginisotope*:1]-C1CC1",
            Chem.BondType.SINGLE: "[998*:1]-[999*:2]>>[beginisotope*:1]1-[endisotope*:2]-C-1",
            Chem.BondType.DOUBLE: "[998*:1]=[999*:2]>>[beginisotope*:1]1=[endisotope*:2]-C-1",
            Chem.BondType.AROMATIC: "[998*:1]~[999*:2]>>[beginisotope*:1]1:[beginisotope*:2]-C-1",
        },
        # These types of templates need more work, if an NC bond is clicked, the addition can be non-sanitizable due to the explicit H (or vice versa!)
        # "0-pyrrole": {
        #     "atom": "[998*:1]>>[beginisotope*:1]-[N]1-C=C-C=C-1",
        #     Chem.BondType.SINGLE: "[998*:1]-[999*:2]>>[beginisotope*:1]1-[endisotopeN:2]-C=C-C=1",
        #     Chem.BondType.DOUBLE: "[998*:1]=[999*:2]>>[beginisotope*:1]1:[endisotopeN:2]:C:C:C:1",  # Not Kekulized
        #     Chem.BondType.AROMATIC: "[998*:1]~[999*:2]>>[beginisotope*:1]1:[beginisotopeN:2]:c:c:c:1",
        # },
        # "1-pyrrole": {
        #     "atom": "[998*:1]>>[beginisotope*:1]-C1=C-C=C-[NH]1",
        #     Chem.BondType.SINGLE: "[998*:1]-[999*:2]>>[beginisotope*:1]1:[endisotope*:2]:[NH]:C:C:1",  # Not kekulized
        #     Chem.BondType.DOUBLE: "[998*:1]=[999*:2]>>[beginisotope*:1]1=[endisotope*:2]-[NH]-C=C-1",
        #     Chem.BondType.AROMATIC: "[998*:1]~[999*:2]>>[beginisotope*:1]1:[beginisotope*:2]:[nH]:c:c:1",
        # },
        # "2-pyrrole": {
        #     "atom": "[998*:1]>>[beginisotope*:1]-C1-C=C-[NH]-C=1",
        #     Chem.BondType.SINGLE: "[998*:1]-[999*:2]>>[beginisotope*:1]1-[endisotope*:2]=C-[NH]-C=1",
        #     Chem.BondType.DOUBLE: "[998*:1]=[999*:2]>>[beginisotope*:1]1:[endisotope*:2]:C:[NH]:C:1",  # Not kekulized
        #     Chem.BondType.AROMATIC: "[998*:1]~[999*:2]>>[beginisotope*:1]1:[beginisotope*:2]:c:[nH]:c:1",
        # },
    }

    @property
    def templateslabels(self):
        return tuple(self.templates.keys())

    def apply_template_to_atom(self, beginatom: Chem.rdchem.Atom, templatelabel: str) -> Chem.Mol:
        """Apply to an atom"""
        beginisotope = beginatom.GetIsotope()

        # TODO assert that the molecule doesnt have this isotope numbers!
        beginatom.SetIsotope(998)

        template_set = self.templates[templatelabel]

        if "atom" not in template_set:
            raise ValueError(f"Template {templatelabel} not supported by atom click")

        template = template_set["atom"]

        template = template.replace("beginisotope", str(beginisotope))

        templateaddition = AllChem.ReactionFromSmarts(template)

        newmol = self.react_and_keep_fragments(beginatom.GetOwningMol(), templateaddition)

        # products = templateaddition.RunReactants((beginatom.GetOwningMol(),))

        # if len(products) != 1:
        #     print(f"Warning: number of products of template addition different than expected, was {len(products)}")

        # product = products[0][0]

        beginatom.SetIsotope(beginisotope)

        return newmol

    def apply_template_to_bond(self, bond: Chem.rdchem.Bond, templatelabel: str) -> Chem.Mol:
        """Apply to a bond"""
        # Reverse is for future usage for unsymmetric templates. Ctrl-Z and reapply to add the template reversed
        beginatom = bond.GetEndAtom() if self.reverse else bond.GetBeginAtom()
        endatom = bond.GetBeginAtom() if self.reverse else bond.GetEndAtom()
        self.reverse = not self.reverse

        beginisotope = (
            beginatom.GetIsotope()
        )  # TODO, hmm we are manipulating the parent molecule here. Can it be avoided?
        endisotope = endatom.GetIsotope()

        # TODO assert that the molecule doesnt have these isotope numbers!
        beginatom.SetIsotope(998)
        endatom.SetIsotope(999)

        bondtype = bond.GetBondType()

        template_set = self.templates[templatelabel]

        if bondtype not in template_set:
            raise ValueError(f"Template {templatelabel} not supported by clicked bondtype {bondtype}")

        template = template_set[bondtype]

        template = template.replace("beginisotope", str(beginisotope)).replace("endisotope", str(endisotope))

        templateaddition = AllChem.ReactionFromSmarts(template)

        mol = bond.GetOwningMol()

        # We need to split fragments, then find the one with the isotopes and react them, then recombine

        newmol = self.react_and_keep_fragments(mol, templateaddition)

        # templateaddition.RunReactantInPlace(mol, removeUnmatchedAtoms=False)

        # products = templateaddition.RunReactants()

        # if len(products) != 1:
        #     print(f"Warning: number of products of template addition different than expected, was {len(products)}")

        # product = products[0][0]

        beginatom.SetIsotope(beginisotope)
        endatom.SetIsotope(endisotope)

        return newmol

    def apply_template_to_canvas(self, mol: Chem.Mol, point: Point2D, templatelabel: str) -> Chem.Mol:
        """Apply to canvas"""
        template = Chem.MolFromSmiles(self.templates[templatelabel]["canvas"])

        if mol.GetNumAtoms() == 0:
            point.x = 0.0
            point.y = 0.0

        combined = Chem.rdchem.RWMol(Chem.CombineMols(mol, template))
        # This should only trigger if we have an empty canvas
        if not combined.GetNumConformers():
            rdDepictor.Compute2DCoords(combined)
        conf = combined.GetConformer(0)
        p3 = Point3D(point.x, point.y, 0)
        conf.SetAtomPosition(mol.GetNumAtoms(), p3)
        return combined

    def react_and_keep_fragments(self, mol, rxn):
        # Split the molecule into fragments
        fragments = Chem.GetMolFrags(mol, asMols=True)

        # Identify the fragment that will react (you may need to modify this logic)
        # reacting_fragment = fragments[0]  # Assume the first fragment is the one reacting

        # Perform the reaction on the reacting fragment
        for i, reacting_fragment in enumerate(fragments):
            products = rxn.RunReactants((reacting_fragment,))
            if products:
                fragments.pop(i)
                break

        if products:
            # Combine the reacted product with the unchanged fragments
            result = products[0][0]
            for frag in fragments:
                result = Chem.CombineMols(result, frag)

            return result
        else:
            return None

    def react_and_keep_fragments(self, mol, rxn):
        # Split the molecule into fragments
        fragments = list(Chem.GetMolFrags(mol, asMols=True))

        # Perform the reaction on the reacting fragment
        for i, reacting_fragment in enumerate(fragments):
            products = rxn.RunReactants((reacting_fragment,))
            if products:
                fragments.pop(i)
                break

        if products:
            # Combine the reacted product with the unchanged fragments
            result = products[0][0]
            for frag in fragments:
                result = Chem.CombineMols(result, frag)

            return result
        else:
            return None
