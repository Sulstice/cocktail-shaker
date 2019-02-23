#!/usr/bin/env python
#
# Runs R Group Converter for Library Creation
#
# ----------------------------------------------------------

# imports
# ---------
from rdkit import Chem
import ruamel.yaml as yaml


# Load datasources
# -------------
with open("datasources/R_Groups.yaml") as stream:
    try:
        R_GROUP_DATASOURCE = yaml.safe_load(stream)
        R_GROUPS = R_GROUP_DATASOURCE['R_Groups']
    except yaml.YAMLError as exc:
        print(exc)


class RaiseMoleculeError(Exception):

    """

    Raise Molecule Error if for some reason we can't evaluate a SMILES, 2D, or 3D molecule.

    """
    def __init__(self, message, errors):
        super().__init__(message)
        self.errors = errors

class RGroupMolObject(object):
    """

    This class is used to take in a molecule and replace any R groups with a one of the groups from the R-Group.

    """

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

        Arguments:
            self (Object): Class RGroupMolObject
            smiles (string): Smiles string that needs to be evaluated
        Returns:
            N / A

        Exceptions:
            RaiseMoleculeError (Exception): MolVs Stacktrace and the smiles string that failed.

        """

        # Use MolVS to validate the smiles to make sure enumeration and r group connections are correct
        # at least in the 1D Format.
        from molvs import validate_smiles

        try:
            validate_smiles(self.original_smiles)
        except RaiseMoleculeError as RME:
            print ("Not a Valid Smiles, Please check the formatting: %s" % self.original_smiles)
            print ("MolVs Stacktrace %s" % RME)

    def r_group_enumerator(self):
        """

        This method will run through the R_Group List to generate smiles strings.

        Arguments:
            self (Object): Class RGroupMolObject

        """

        CHEMICAL_SMILES_R_GROUPS = set()

        for key, value in R_GROUPS.items():
            if value in self.original_smiles:
                r_group_molecule = RGroupMolObject(Chem.MolFromSmiles(value))
                # delete part of the molecule containing the structure.
                truncated_molecule = Chem.DeleteSubstructs(self.molecule, r_group_molecule.molecule)
                filtered_r_groups = {filtered_key: filtered_value for filtered_key, filtered_value \
                                     in R_GROUPS.items() if not filtered_value == value }

                for filtered_key, filtered_value in filtered_r_groups.items():
                    combination_molecule = Chem.CombineMols(
                                                truncated_molecule,
                                                RGroupMolObject(Chem.MolFromSmiles(filtered_value)).molecule)
                    editable_combination = Chem.EditableMol(combination_molecule)


                    for i in range(0, Chem.Mol.GetNumAtoms(truncated_molecule)):

                        editable_combination.AddBond(i, 9, order=Chem.rdchem.BondType.SINGLE)

                        molecule_returned = editable_combination.GetMol()
                        CHEMICAL_SMILES.add(Chem.MolToSmiles(molecule_returned))

        return CHEMICAL_SMILES_R_GROUPS

if __name__ == "__main__":
        scaffold_molecule = RGroupMolObject(Chem.MolFromSmiles('c1cc(CCCO)ccc1'))
        chemcials = scaffold_molecule.r_group_enumerator()
        from enumeration import enumerate_smiles