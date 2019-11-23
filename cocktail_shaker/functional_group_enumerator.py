#!/usr/bin/env python
#
# Runs R Group Converter for Library Creation
#
# ----------------------------------------------------------

# Imports
# ---------
from rdkit import Chem
import ruamel.yaml as yaml
import itertools
import progressbar

# Cocktail Shaker Imports
# -----------------------
from validation import MoleculeValidator

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

class Cocktail(object):
    """

    This class is used to take in a molecule and replace any R groups with a one of the groups from the R-Group.

    """

    __version_parser__ = 1.0
    __allow_update__ = False

    def __init__(self, peptide_backbone, ligand_library = [], enumeration = False, dimensionality = '1D', enumeration_complexity='Low'):

        # imports
        # -------
        import re

        # I will allow the user to pass a string but for easier sake down the road
        # I will reimplement it as a list.

        load_datasources()
        self.peptide_backbone = str(peptide_backbone)
        self.ligand_library = ligand_library
        self.peptide_backbone_length = int(max(re.findall("\d", self.peptide_backbone)))

        if self.peptide_backbone_length > len(self.ligand_library):
            print ("Cocktail Shaker Error: Peptide Backbone Length needs to be less than or equal to for your library")
            raise IndexError

        if enumeration:
            self.enumerate(dimensionality, enumeration_complexity)

    def shake(self):

        """

        Generate all combinations of a molecule

        """

        results = []
        peptide = self.peptide_backbone
        combinations = list(itertools.permutations(self.ligand_library, self.peptide_backbone_length))

        print ("Generating Compounds...")

        for i in progressbar.progressbar(range(len(combinations))):
            combination = list(combinations[i])
            peptide_molecule = str(peptide)
            for j in range(0, len(combination)):

                modified_molecule = Chem.ReplaceSubstructs(Chem.MolFromSmiles(peptide_molecule),
                                                               Chem.MolFromSmiles('[*:'+ str(j+1)+']'),
                                                               Chem.MolFromSmiles(str(combination[j])))

                peptide_molecule = Chem.MolToSmiles(modified_molecule[0])

            results.append(str(peptide_molecule))

        # Remove Duplicates
        results = list(set(results))

        # Validate Smiles
        MoleculeValidator(results, smiles=True)

        return results

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
                    print ('3D Functionality is not supported yet!')
                    return enumerated_molecules

        return enumerated_molecules



