#!/usr/bin/env python
#
# Runs R Group Converter for Library Creation
#
# ----------------------------------------------------------

# Imports
# ---------
from rdkit import Chem
import ruamel.yaml as yaml

# Relative Imports
# ----------------
from validation import RaiseMoleculeError

# Load datasources
# ----------------
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


class LoadCustomLibrary(object):

    """

    Load a custom datasource per the user.

    Convert their SMILES to SMARTS

    """

    __version__ = "1.0.2"


    def __init__(self, molecules):

        # imports
        # -------

        from molvs import Validator
        from validation import MoleculeValidator

        self.molecules = molecules
        self.custom_datasources = False

        validator = MoleculeValidator(self.molecules, smiles=True)


class Cocktail(object):
    """

    This class is used to take in a molecule and replace any R groups with a one of the groups from the R-Group.

    """

    __version_parser__ = 1.0
    __allow_update__ = False

    def __init__(self, peptide_backbone, ligand_library = [], dimensionality = '1D', enumeration_complexity='Low'):

        # imports
        # -------
        import re
        from rdkit import Chem
        from molvs import Validator

        # I will allow the user to pass a string but for easier sake down the road
        # I will reimplement it as a list.

        load_datasources()
        self.peptide_backbone = peptide_backbone
        self.ligand_library = ligand_library

        # Guardrail for combinations
        # max_number_side_chains = [character for character in re.split("[^0-9]", self.peptide_backbone) if character != '']
        # if int(max(map(int, max_number_side_chains))) != len(self.ligand_library):
        #     print ("Ligand Library must match the amount of Side Chains")
        #     raise IndexError

        self.dimensionality = dimensionality
        self.enumeration_complexity = enumeration_complexity

    def shake(self):

        """

        Generate all combinations of a molecule

        """

        results = []
        peptide = self.peptide_backbone
        import itertools
        combinations = (list(itertools.permutations(self.ligand_library)))
        print (combinations)

        for combination in combinations:
            combination = list(combination)
            peptide_molecule = ''
            for j in range(0, len(self.ligand_library)):
                if j == 0:
                    peptide_molecule = str(peptide).replace('[*:' + str(j + 1) +']', combination[j])
                else:
                    peptide_molecule = str(peptide_molecule).replace('[*:' + str(j + 1) +']', combination[j])
            results.append(peptide_molecule)

        print (results)
        return results


    def _shake(self, functional_groups=["all"], shape=None):

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
        modified_molecules = []
        if functional_groups[0] == 'all':
            for molecule in self.molecules:
                for key, value in patterns_found.items():
                        smarts_mol = Chem.MolFromSmarts(value[1])
                        for functional_group, pattern in R_GROUPS.items():
                            if functional_group == 'R-Groups':
                                continue
                            else:
                                for i in range(0, len(pattern)):
                                    for r_functional_group, r_data in pattern[0].items():
                                        if r_data[1] == value[1]:
                                            continue
                                        try:
                                            if self.markush_structure:
                                                modified_molecule = Chem.ReplaceSubstructs(molecule, Chem.MolFromSmiles('[*:1]'),
                                                                                           Chem.MolFromSmiles(r_data[0]), replaceAll=True)
                                            else:
                                                modified_molecule = Chem.ReplaceSubstructs(molecule, smarts_mol,
                                                                                          Chem.MolFromSmiles(r_data[0]), replaceAll=True)
                                            modified_molecules.append(modified_molecule[0])
                                        except RaiseMoleculeError:
                                            print ("Molecule Formed is not possible")
        else:
            for molecule in self.molecules:
                for key, value in patterns_found.items():
                    smarts_mol = Chem.MolFromSmarts(value[1])
                    for functional_group, pattern in R_GROUPS.items():
                        if functional_group == 'R-Groups':
                            continue
                        elif functional_group not in functional_groups:
                            continue
                        else:
                            for r_functional_group, r_data in pattern[0].items():
                                # Skip redundacies if the r group is already matched.
                                if r_data[1] == value[1]:
                                    continue
                                try:
                                    if self.markush_structure:
                                        modified_molecule = Chem.ReplaceSubstructs(molecule, Chem.MolFromSmiles('[*:1]'),
                                                                                   Chem.MolFromSmiles(r_data[0]))
                                    else:
                                        modified_molecule = Chem.ReplaceSubstructs(molecule, smarts_mol,
                                                                                   Chem.MolFromSmiles(r_data[0]))
                                    modified_molecules.append(modified_molecule[0])
                                except RaiseMoleculeError:
                                    print ("Molecule Formed is not possible")
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
                    print ('3D Functionality is not supported yet!')
                    return enumerated_molecules

        return enumerated_molecules



