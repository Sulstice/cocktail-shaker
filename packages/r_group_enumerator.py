#!/usr/bin/env python
#
# Runs R Group Converter for Library Creation
#
# ----------------------------------------------------------

# imports
# ---------
from rdkit import Chem
import ruamel.yaml as yaml

from packages.file_handling import FileWriter, FileParser

# Load datasources
# -------------
def load_datasources():

    """

    Load all the datasources for running this package in local context.

    This might slow down performance later -- we can opt in to load sources of data dependent on the functional group.

    """
    from pathlib import Path
    datasource_location = Path(__file__).absolute().parent
    with open(str(datasource_location) + "/datasources/R_Groups.yaml") as stream:
        try:
            global R_GROUP_DATASOURCE
            R_GROUP_DATASOURCE = yaml.safe_load(stream)

            global R_GROUPS
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

    def __init__(self, molecules):
        from molvs import Validator

        # I will allow the user to pass a string but for easier sake down the road
        # I will reimplement it as a list.
        if type(molecules) is not list:
            self.molecules = [molecules]
        else:
            self.molecules = molecules

        for molecule in self.molecules:
            self.original_smiles = Chem.MolToSmiles(molecule)
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

    def validate_molecule(self, molecule):

        """

        This function will be used to validate molecule objects

        Arguments:
            self (Object): Class RGroupMolObject
            molecule (RDKit Object): Molecule object we need to sanitize.
        Returns:
            molecule (RDKit Object): The RDkit Molecule molecule object

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

        for molecule in self.molecules:
            for key, value in R_GROUPS.items():
                pattern = Chem.MolFromSmarts(value[1])
                if molecule.GetSubstructMatches(pattern,uniquify=False):
                    print ("Found Functional Group: %s | Pattern Count: %s" % (key,
                                                                               len(molecule.GetSubstructMatches(
                                                                                   pattern,uniquify=False))))
                    pattern_payload[key] = [value[0], value[1]]

        return pattern_payload

    def r_group_enumerator(self, patterns_found):

        """

        TODO: Do this faster than O(n)^3 as this algorithm is not the most efficient.
        """

        modified_molecules = []
        for molecule in self.molecules:
            for key, value in patterns_found.items():
                    smarts_mol = Chem.MolFromSmarts(value[1])
                    for r_functional_group, r_data in R_GROUPS.items():
                        # Skip redundacies if the r group is already matched.
                        if r_data[1] == value[1]:
                            continue

                        modified_molecule = Chem.ReplaceSubstructs(molecule, smarts_mol,
                                                                  Chem.MolFromSmiles(r_data[0]), replaceAll=True)

                        modified_molecules.append(modified_molecule[0])

        return modified_molecules

if __name__ == "__main__":

        import optparse
        load_datasources()

        # Have the option to parse in a file if needed.
        parser = optparse.OptionParser()
        parser.add_option('-f', '--file',
                  action="store", dest="file",
                  help="the file path",  nargs='?')

        # Have the option to enumerate later.
        parser.add_option('-e', '--enumerate',
                          action="store", dest="enumerate",
                          help="enumerate all representations of file pass in options: True | False",
                          nargs='?')
        options, args = parser.parse_args()

        if options.file:
            file = FileParser(options.file)

        scaffold_molecule = RGroupMolObject([Chem.MolFromSmiles('c1cc(CCCO)ccc1'), Chem.MolFromSmiles('c1cc(CCCBr)ccc1')])
        patterns_found = scaffold_molecule.find_r_groups()
        modified_molecules = scaffold_molecule.r_group_enumerator(patterns_found=patterns_found)
        FileWriter("test", modified_molecules, "sdf")
        # FileWriter("test", modified_molecules, "txt")

