from rdkit import Chem
from typing import List

def isVirtualAtom(atom: Chem.Atom):
    if getAtomMapNumber(atom) == None:
        return False
    return True

def getAtomMapNumber(atom: Chem.Atom):
    if atom.HasProp("molAtomMapNumber"):
        return atom.GetIntProp("molAtomMapNumber")
    if atom.GetAtomicNum() == 0:
        return 0
    return None

def toRealAtom(atom: Chem.Atom):
    a = Chem.Atom(atom)

    if a.HasProp("molAtomMapNumber"):
        a.ClearProp("molAtomMapNumber")
    
    return a

def getOtherBonds(atom: Chem.Atom):
    bonds = []

    bond: Chem.Bond
    for bond in atom.GetBonds():
        otherAtom = bond.GetOtherAtom(atom)

        if atom.GetAtomicNum() != 0 and isVirtualAtom(otherAtom):
            continue

        bonds.append(BondTemplate(atom, bond))

    return bonds

def addTemplate(mol: Chem.RWMol, template: Chem.Mol):
    atomIdxMap = []

    atom: Chem.Atom
    for atom in template.GetAtoms():
        atomIdxMap.append(mol.AddAtom(atom))
    
    bond: Chem.Bond
    for bond in template.GetBonds():
        mol.AddBond(atomIdxMap[bond.GetBeginAtomIdx()], atomIdxMap[bond.GetEndAtomIdx()], bond.GetBondType())

    return atomIdxMap

class AtomRef:
    def __init__(self, atom: Chem.Atom):
        self.atom = atom
        self.atomMapNumber = getAtomMapNumber(atom)
        self.atomIsVirtual = isVirtualAtom(atom)
    
    def getAtomIdx(self, atomIdxMap: List[int], virtualAtomMap):
        if not self.atomIsVirtual:
            return atomIdxMap[self.atom.GetIdx()]
        
        return virtualAtomMap[self.atomMapNumber]

class BondTemplate:
    def __init__(self, atom: Chem.Atom, bond: Chem.Bond):
        self.atom = AtomRef(atom)
        self.other = AtomRef(bond.GetOtherAtom(atom))
        self.bondType = bond.GetBondType()
    
    def makeBond(self, mol: Chem.RWMol, atomIdxMap: List[int], virtualAtomMap):
        atomIdx = self.atom.getAtomIdx(atomIdxMap, virtualAtomMap)
        otherIdx = self.other.getAtomIdx(atomIdxMap, virtualAtomMap)
        mol.AddBond(atomIdx, otherIdx, self.bondType)

class ProductTemplate:
    def __init__(self, smarts):
        self.mol = Chem.RWMol(Chem.MolFromSmarts(smarts))
        self.virtualAtoms = {}
        self.virtualBonds = []
        atomsToRemove = []

        atom: Chem.Atom
        for atom in self.mol.GetAtoms():
            atomMapNumber = getAtomMapNumber(atom)
            if atomMapNumber != None:
                self.virtualAtoms[atomMapNumber] = (toRealAtom(atom), getOtherBonds(atom))
                atomsToRemove.append(atom)

        bond: Chem.Bond
        for bond in self.mol.GetBonds():
            beginAtom = bond.GetBeginAtom()
            beginAtomMapNumber = getAtomMapNumber(beginAtom)

            if beginAtomMapNumber == None or beginAtom.GetAtomicNum() == 0:
                continue

            endAtom = bond.GetEndAtom()
            endAtomMapNumber = getAtomMapNumber(endAtom)

            if endAtomMapNumber == None or endAtom.GetAtomicNum() == 0:
                continue

            self.virtualBonds.append(BondTemplate(beginAtom, bond))
        
        for atom in atomsToRemove:
            self.mol.RemoveAtom(atom.GetIdx())

        self.maxOrder = len(self.virtualAtoms)

        if self.maxOrder > 1 and 0 in self.virtualAtoms:
            self.maxOrder -= 1

    def apply(self, mol: Chem.Mol, atomRefs: List[Chem.Atom], order: int = 0):
        rwmol = Chem.RWMol(mol)
        return self.applyRW(rwmol, list(map(lambda atom: rwmol.GetAtomWithIdx(atom.GetIdx()), atomRefs)), order)
    
    def applyRW(self, mol: Chem.RWMol, atomRefs: List[Chem.Atom], order: int = 0):
        if order > self.maxOrder:
            return None

        startIdx = 1 if order > 1 or 0 not in self.virtualAtoms else 0
        
        atomIdxMap = addTemplate(mol, self.mol)

        virtualAtomIdxMap = {}
        for i, atom in enumerate(atomRefs[:order]):
            virtualAtomIdxMap[i + startIdx] = atomRefs[i].GetIdx()

        atom: Chem.Atom
        for i, (atom, bonds) in self.virtualAtoms.items():
            if atom.GetAtomicNum() == 0:
                continue

            if i in virtualAtomIdxMap:
                continue

            virtualAtomIdxMap[i] = mol.AddAtom(atom)
        
        for i, (atom, bonds) in self.virtualAtoms.items():
            if order != 1 and i == 0:
                continue

            for bond in bonds:
                bond.makeBond(mol, atomIdxMap, virtualAtomIdxMap)

        if order <= 1:
            for bond in self.virtualBonds:
                bond.makeBond(mol, atomIdxMap, virtualAtomIdxMap)

        return mol

if __name__ == '__main__':
    mol = Chem.MolFromSmarts("c1ccccc1")
    productTemplate = ProductTemplate("*[c:1]1[c:2]cccc1")

    for i in range(3):
        product = productTemplate.apply(mol, [mol.GetAtomWithIdx(0), mol.GetAtomWithIdx(1)], order=i)
        Chem.SanitizeMol(product)
        Chem.Kekulize(product)
        print(i, Chem.MolToSmiles(product))
