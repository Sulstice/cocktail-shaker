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

    __version_error_parser__ = 1.0
    __allow_update__ = False

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
        from molvs import Validator
        if molecule != None:
            self.molecule = molecule
            self.original_smiles = Chem.MolToSmiles(self.molecule)

            # Validation
            validator_format = '%(asctime)s - %(levelname)s - %(validation)s - %(message)s'
            self.validate = Validator(log_format=validator_format)

    def validate_smiles(self, smiles):
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
        from molvs import validate_smiles as vs

        try:
            vs(smiles)
        except RaiseMoleculeError as RME:
            print ("Not a Valid Smiles, Please check the formatting: %s" % self.original_smiles)
            print ("MolVs Stacktrace %s" % RME)

    def validate_molecule(self, molecule=None):

        """

        This function will be used to validate molecule objects

        Arguments:
            self (Object): Class RGroupMolObject
            molecule (RDKit Object): Molecule object we need to sanitize.
        Returns:
            N / A

        Exceptions:
            RaiseMoleculeError (Exception): Raise the Raise Molcule Error if the molecule is not valid.

        TODO: Verify Sanitize molcule that the validation works
        """

        if not molecule:
            try:
                Chem.rdmolops.SanitizeMol(molecule)
            except RaiseMoleculeError as RME:
                print ("Not a valid molecule: %s" % RME)
            finally:
                return molecule

    def find_r_groups(self):

        """

        Find functional groups that ligand library loader supports

        :return:

        """

        pattern_payload = {}

        for key, value in R_GROUPS.items():
            pattern = Chem.MolFromSmarts(value[1])
            if self.molecule.GetSubstructMatches(pattern,uniquify=False):
                print ("Found pattern: %s" % len(self.molecule.GetSubstructMatches(pattern,uniquify=False)))
                pattern_payload[value[1]] = len(self.molecule.GetSubstructMatches(pattern,uniquify=False))

    def r_group_enumerator(self):


        # find one R Group

        pattern = Chem.MolFromSmarts('[OX2H]')
        print ("number of matches:", len(self.molecule.GetSubstructMatches(pattern,uniquify=False)))

        # Enumerate through the R Groups stored in the system.
        for key, value in R_GROUPS.items():
            modified_molecule = Chem.ReplaceSubstructs(self.molecule, pattern, Chem.MolFromSmiles(value[0]),
                                                       replaceAll=True)

            # Validate the molecule.
            self.validate_molecule(molecule=modified_molecule[0])

            print (Chem.MolToSmiles(modified_molecule[0]))

if __name__ == "__main__":
        scaffold_molecule = RGroupMolObject(Chem.MolFromSmiles('c1cc(CCCO)ccc1'))
        chem = scaffold_molecule.find_r_groups()
        # chemcials = scaffold_molecule.r_group_enumerator()
