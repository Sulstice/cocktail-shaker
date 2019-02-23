#!/usr/bin/env python
#
# Runs R Group Converter for Library Creation
#
# ----------------------------------------------------------

from rdkit import Chem

class RGroupMolObject(object):
    """

    This class is used to enumerate an object
    """

    R_GROUP_LIST = ["Br", "Cl", "F", "O", "O=C", "OC=O", "N"]

    __version_parser__ = 1.0
    __allow_update__ = False

    def __init__(self, molecule=None):
        if molecule != None:
            self.molecule = molecule
            self.original_smiles = Chem.MolToSmiles(self.molecule)
            self.validate_molecule()

    def validate_molecule(self, smiles=None):
        """

        This method takes the smiles string and runs through the validation check of a smiles string.

        """

        # Use MolVS to validate the smiles to make sure enumeration and r group connections are correct
        # at least in the 1D Format.
        from molvs import validate_smiles

        try:
            if smiles:
                validate_smiles(smiles)
            else:
                validate_smiles(self.original_smiles)
        except RaiseMoleculeError:
            print ("Not a Valid Smiles, Please check the formatting: %s" % self.original_smiles)

    def r_group_enumerator(self):
        """

        This method will run through the R_Group List to generate smiles strings.

        """
        Chem.ReplaceSubstructs(self.molecule)



if __name__ == "__main__":
        start_mol = Chem.MolFromSmiles('c1cc(CCCO)ccc1')
        truncate = Chem.DeleteSubstructs(start_mol, Chem.MolFromSmiles('O'))
        molecule = RGroupMolObject(Chem.MolFromSmiles(i))
        molecule.validate_molecule()


        # combo = Chem.CombineMols(truncate, mod)
        # Chem.MolToSmiles(combo)
        # edcombo = Chem.EditableMol(combo)
        # edcombo.AddBond(3, 9, order=Chem.rdchem.BondType.SINGLE)
        # back = edcombo.GetMol()
        # print(Chem.MolToSmiles(back))
