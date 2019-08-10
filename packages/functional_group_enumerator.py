#!/usr/bin/env python
#
# Runs R Group Converter for Library Creation
#
# ----------------------------------------------------------

# imports
# ---------
from rdkit import Chem
import ruamel.yaml as yaml

from file_handler import FileWriter, FileParser

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

class RequestError(Exception):

    __version_error_parser__ = 1.0
    __allow_update__ = False

    """

    Raise a Request Handler Error if the response code is not 200. 

    """
    def __init__(self, message, errors):
        super().__init__(message)
        self.errors = errors


class RequestHandler(object):


    __version__ = 1.0
    __allow_update__ = False

    """
    
    Handler API requests to external libraries
    
    
    """

    def __init__(self, url):

        self.url = url

    def get(self):
        """

        Issues the request of the URL

        TODO: Handle requests that need a payload (assuming this will be the case for MolPort down the road)

        """

        # imports
        # -------
        from urllib.request import urlopen
        from lxml import etree

        try:
            request = urlopen(self.url)
        except RequestError:
            print ("Error handling request, please try again or contact Suliman Sharif")

        if request.getcode() != 200:
            raise RequestError(message="Error handling request, please try again or contact Suliman Sharif", errors='Request Error')

        response = etree.parse(request).getroot()
        for child in response.iter('item'):
            response = child.text

        return response

class CactusResolver(object):

    __version__ = 1.0
    __allow_update__ = False

    """
    
    Handle responses coming back from cactus. 
    
    """

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

        TODO: Verify Sanitize molcule that the validation works
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

        for molecule in self.molecules:
            for key, value in R_GROUPS.items():
                pattern = Chem.MolFromSmarts(value[1])
                if molecule.GetSubstructMatches(pattern,uniquify=False):
                    print ("Found Functional Group: %s | Pattern Count: %s" % (key,
                                                                               len(molecule.GetSubstructMatches(
                                                                                   pattern,uniquify=False))))
                    pattern_payload[key] = [value[0], value[1]]

        return pattern_payload

    def shake(self, shape=None):

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

        modified_molecules = []
        for molecule in self.molecules:
            for key, value in patterns_found.items():
                    smarts_mol = Chem.MolFromSmarts(value[1])
                    for r_functional_group, r_data in R_GROUPS.items():
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


        return modified_molecules

    def _retrieve_3D_rendering(self):

        """

        TODO: Incorporate TwirlyMol into a python version

        Arguments:
             molecules (List): List of molecules we need to retrieve the 3D rendering.

        Returns:
            three_d_molecules_text (List): A list of the 3D molecules representation in text.

        Exceptions:
            RequestError: Raises when the response code from cactus is not 200 (something wrong on their server side).

        """

        _API_BASE = 'https://cactus.nci.nih.gov/chemical/structure'
        _FILE_FORMATS = ['mol2']

        def _construct_api_url(molecule, get3d=True, **kwargs):
            """

            Return the URL for the desired API endpoint.

            Arguments:
                molecule (String): Smiles representation of the molecule we are interested in.
                get3d (Boolean): Parameter Cactus needs to handle to return the mol2 data.

            Returns:
                url (string): The full url used to request from cactus.

            """

            # imports
            # -------
            from urllib.parse import quote, urlencode

            kwargs['format'] = 'mol2'
            representation = 'file'
            tautomers = 'tautomers:%s' % molecule
            url = '{}/{}/{}/{}'.format(_API_BASE, quote(tautomers), representation, 'xml/')
            if get3d:
                kwargs['get3d'] = True

            url += '?%s' % urlencode(kwargs)

            return url

        three_d_molecules_text = []

        bulk = False

        if len(self.molecules) > 1:
            bulk = True
        for molecule in self.molecules:
            # Our request will use smiles.
            molecule = Chem.MolToSmiles(molecule)
            url = _construct_api_url(molecule)

            # Parse out to the request handler
            request = RequestHandler(url)
            response = request.get()

            three_d_molecules_text.append(response)

        return three_d_molecules_text

    @staticmethod
    def enumerate(molecules, enumeration_complexity=None, dimensionality=None):

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

        if enumeration_complexity.lower() == 'low':
            complexity = 10
        elif enumeration_complexity.lower() == 'medium':
            complexity = 100
        elif enumeration_complexity.lower() == 'high':
            complexity = 1000
        else:
            complexity = 10

        enumerated_molecules = []
        for molecule in molecules:
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

if __name__ == "__main__":

        starting_compound = Cocktail([Chem.MolFromSmiles('c1cc(CCCO)ccc1'), Chem.MolFromSmiles('c1cc(CCCBr)ccc1')])
        starting_compound._retrieve_3D_rendering()
        modified_molecules = starting_compound.shake()
        # smiles = Cocktail.enumerate(modified_molecules, enumeration_complexity='Low', dimensionality='2D')

        # FileWriter("test", modified_molecules, "sdf")
        # FileWriter("test", modified_molecules, "txt")
