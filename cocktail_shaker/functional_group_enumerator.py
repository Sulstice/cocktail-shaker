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
            print ("Datasources not loading correctly, Please contact lead developer")
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

class Cocktail(object):
    """

    This class is used to take in a molecule and replace any R groups with a one of the groups from the R-Group.

    """

    __version_parser__ = 1.0
    __allow_update__ = False

    def __init__(self, molecules):
        from molvs import Validator

        # I will allow the user to pass a string but for easier sake down the road
        # I will reimplement it as a list.

        load_datasources()
        self.dimensionality = '1D'
        self.modified_molecules = []
        rdkit_rendered_molecules = []
        for molecule in molecules:
            rdkit_rendered_molecules.append(Chem.MolFromSmiles(molecule))

        self.molecules = rdkit_rendered_molecules

        for molecule in self.molecules:
            self.original_smiles = Chem.MolToSmiles(molecule)
            # Validation
            validator_format = '%(asctime)s - %(levelname)s - %(validation)s - %(message)s'
            self.validate = Validator(log_format=validator_format)

    def validate_smiles(self, smiles):
        """

        This method takes the smiles string and runs through the validation check of a smiles string.

        Arguments:
            self (Object): Class Cocktail
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
            return True
        except RaiseMoleculeError as RME:
            print ("Not a Valid Smiles, Please check the formatting: %s" % self.original_smiles)
            print ("MolVs Stacktrace %s" % RME)

        return False

    def validate_molecule(self, molecule):

        """

        This function will be used to validate molecule objects

        Arguments:
            self (Object): Class Cocktail
            molecule (RDKit Object): Molecule object we need to sanitize.
        Returns:
            molecule (RDKit Object): The RDkit Molecule molecule object

        Exceptions:
            RaiseMoleculeError (Exception): Raise the Raise Molcule Error if the molecule is not valid.

        TODO: Verify Sanitize molecule that the validation works
        """

        if not molecule:
            try:
                Chem.rdmolops.SanitizeMol(molecule)
            except RaiseMoleculeError as RME:
                print ("Not a valid molecule: %s" % RME)
            finally:
                return molecule

    def detect_functional_groups(self):

        """

        Find functional groups that ligand library loader supports

        :return:

        """

        pattern_payload = {}
        load_datasources()

        print ("Detecting Functional Groups...")

        for molecule in self.molecules:
            for functonal_group, pattern in R_GROUPS.items():
                for key, value in pattern[0].items():
                    smart_pattern = Chem.MolFromSmarts(value[1])
                    if molecule.GetSubstructMatches(smart_pattern,uniquify=False):
                        print ("Found Functional Group: %s | Pattern Count: %s" % (key,
                                                                                   len(molecule.GetSubstructMatches(
                                                                                       smart_pattern,uniquify=False))))
                    pattern_payload[key] = [value[0], value[1]]

        return pattern_payload

    def shake(self, functional_groups=["all"], shape=None):

        """

        Used to swap out molecules based on the patterns found from the detection.

        Arguments:
            self (Object): Cocktail object of the list of molecules

        Return:
            modified_molecules (List): List of the RDKit molecule objects that have had their structures replaced.

        TODO: Do this faster than O(n)^3 as this algorithm is not the most efficient.

        """

        # Run detection first to see and validate what functional groups have been found.

        patterns_found = self.detect_functional_groups()
        print ("Shaking Compound....")

        if functional_groups[0] == 'all':
            modified_molecules = []
            for molecule in self.molecules:
                for key, value in patterns_found.items():
                        smarts_mol = Chem.MolFromSmarts(value[1])
                        for functonal_group, pattern in R_GROUPS.items():
                            for r_functional_group, r_data in pattern[0].items():
                                # Skip redundacies if the r group is already matched.
                                if r_data[1] == value[1]:
                                    continue
                                try:
                                    modified_molecule = Chem.ReplaceSubstructs(molecule, smarts_mol,
                                                                              Chem.MolFromSmiles(r_data[0]), replaceAll=True)
                                    modified_molecules.append(modified_molecule[0])
                                except RaiseMoleculeError:
                                    print ("Molecule Formed is not possible")
                                    # continue
        else:
            modified_molecules = []
            for molecule in self.molecules:
                for key, value in patterns_found.items():
                    smarts_mol = Chem.MolFromSmarts(value[1])
                    for functional_group, pattern in R_GROUPS.items():
                        if functional_group not in functional_groups:
                            continue
                        else:
                            for r_functional_group, r_data in pattern[0].items():
                                # Skip redundacies if the r group is already matched.
                                if r_data[1] == value[1]:
                                    continue
                                try:
                                    modified_molecule = Chem.ReplaceSubstructs(molecule, smarts_mol,
                                                                               Chem.MolFromSmiles(r_data[0]), replaceAll=True)
                                    modified_molecules.append(modified_molecule[0])
                                except RaiseMoleculeError:
                                    print ("Molecule Formed is not possible")
                                    # continue
        self.modified_molecules = modified_molecules

        print ("Molecules Generated: {}".format(len(modified_molecules)))

        return modified_molecules

    def enumerate(self, enumeration_complexity=None, dimensionality=None):

        """

        Enumerate the drug library based on dimension.

        Arguments:
            molecules (List): a list of molecules that the user would like enumerated.
            enumeration_complexity (String): Declares how many times we will want to discover another molecule
                                             configuration
            dimensionality (String): Enumerate based on dimensionality (1D, 2D, 3D)

        Returns:
            enumerated_molecules (List): Dependent on the dimensionality of the user it can be -> a list of smiles, or
                                         a list of RDKit Molecule Objects.

        """

        # Enumeration comes from the user iwatobipen
        # https://iwatobipen.wordpress.com/2018/11/15/generate-possible-list-of-smlies-with-rdkit-rdkit/

        print ("Enumerating Compunds....")

        if enumeration_complexity.lower() == 'low':
            complexity = 10
        elif enumeration_complexity.lower() == 'medium':
            complexity = 100
        elif enumeration_complexity.lower() == 'high':
            complexity = 1000
        else:
            complexity = 10

        enumerated_molecules = []
        for molecule in self.modified_molecules:
            for i in range(complexity):
                smiles_enumerated = Chem.MolToSmiles(molecule, doRandom=True)
                if dimensionality == '1D' and smiles_enumerated not in enumerated_molecules:
                    enumerated_molecules.append(smiles_enumerated)
                elif dimensionality == '2D':
                    if not smiles_enumerated in enumerated_molecules:
                        enumerated_molecules.append(Chem.MolFromSmiles(smiles_enumerated))
                elif dimensionality == '3D':
                    print ('Not supported yet')
                    return enumerated_molecules

        return enumerated_molecules



