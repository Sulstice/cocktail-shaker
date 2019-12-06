#!/usr/bin/env python
#
# Runs R Group Converter for Library Creation
#
# ----------------------------------------------------------

# Imports
# ---------
from rdkit import Chem

# Suppress RDKit's StdOut and StdErr
from rdkit import RDLogger
lg = RDLogger.logger()
lg.setLevel(RDLogger.CRITICAL)


import ruamel.yaml as yaml
import itertools
import progressbar

# Cocktail Shaker Imports
# -----------------------
from .validation import MoleculeValidator

class Cocktail(object):
    """

    This class is used to take in a molecule and replace any R groups with a one of the groups from the R-Group.

    """

    __version_parser__ = "1.1.0"
    __allow_update__ = False

    def __init__(self, peptide_backbone, ligand_library = [], enable_isomers = False, include_amino_acids = False):

        """

        Initialize the Cocktail Object with optional parameters

        Arguments:
            peptide_backbone (String): Peptide Backbone with the slots allocated.
            ligand_library (List): The library of ligands the user would like for the combinations
            enable_isomers (Bool): Whether the user would like to include stereoisomers in the results
            include_amino_acids (Bool): Whether the user would like to include amino acids in their combinations

        """

        # imports
        # -------
        import re

        # I will allow the user to pass a string but for easier sake down the road
        # I will reimplement it as a list.

        self.peptide_backbone = str(peptide_backbone)
        print(peptide_backbone)
        self.ligand_library = ligand_library

        # Detect the proline amino acid on the N-terminus, set the max peptide length accordingly.
        if self.peptide_backbone[0:6] == 'N2CCCC2':
            self.peptide_backbone_length = int(max(map(int, re.findall(r"\d+", self.peptide_backbone[7:]))))

        self.peptide_backbone_length = int(max(map(int, re.findall(r"\d+", self.peptide_backbone))))
        self.enable_isomers = enable_isomers
        self.include_amino_acids = include_amino_acids

        self.combinations = []

        if (self.peptide_backbone_length > len(self.ligand_library) or not ligand_library) and self.include_amino_acids == False:
            print ("Cocktail Shaker Error: Peptide Backbone Length needs to be less than or equal to for your library")
            raise IndexError

    def shake(self):

        """

        Generate all combinations of a molecule

        """

        results = []
        peptide = self.peptide_backbone

        if self.include_amino_acids:
            if not self.ligand_library:
                self.ligand_library = self._load_amino_acids()
            else:
                self.ligand_library.extend(self._load_amino_acids())

        combinations = list(itertools.permutations(self.ligand_library, self.peptide_backbone_length))

        print ("Generating Compounds...")

        for i in progressbar.progressbar(range(len(combinations))):
            combination = list(combinations[i])
            peptide_molecule = str(peptide)
            for j in range(0, len(combination)):

                modified_molecule = Chem.ReplaceSubstructs(Chem.MolFromSmiles(peptide_molecule),
                                                               Chem.MolFromSmiles('[*:'+ str(j+1)+']'),
                                                               Chem.MolFromSmiles(str(combination[j])))

                peptide_molecule = Chem.MolToSmiles(modified_molecule[0], isomericSmiles=True)

            # Enable StereoChemistry
            if self.enable_isomers:

                molecule = Chem.MolFromSmiles(str(peptide_molecule))
                options = Chem.EnumerateStereoisomers.StereoEnumerationOptions(unique=True, tryEmbedding=True)
                isomers = tuple(Chem.EnumerateStereoisomers.EnumerateStereoisomers(molecule, options=options))
                for smile in sorted(Chem.MolToSmiles(isomer, isomericSmiles=True) for isomer in isomers):
                    results.append(str(smile))
            else:
                results.append(str(peptide_molecule))

        # Remove Duplicates
        results = list(set(results))

        # Validate Smiles
        MoleculeValidator(results, smiles=True)

        # Store in combinations if enumeration
        self.combinations = results

        return results

    def enumerate(self, dimensionality = '1D', enumeration_complexity='Low'):

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
            complexity = 100
        elif enumeration_complexity.lower() == 'medium':
            complexity = 1000
        elif enumeration_complexity.lower() == 'high':
            complexity = 10000
        else:
            complexity = 10

        enumerated_molecules = []
        for i in progressbar.progressbar(range(len(self.combinations))):
            for _ in range(complexity):
                molecule = Chem.MolFromSmiles(self.combinations[i])
                smiles_enumerated = Chem.MolToSmiles(molecule, doRandom=True)
                if dimensionality == '1D' and smiles_enumerated not in enumerated_molecules:
                    enumerated_molecules.append(smiles_enumerated)
                elif dimensionality == '2D':
                    if not smiles_enumerated in enumerated_molecules:
                        enumerated_molecules.append(Chem.MolFromSmiles(smiles_enumerated))
                elif dimensionality == '3D':
                    print ('3D Functionality is not supported yet!')
                    return enumerated_molecules

        # Validate Smiles
        MoleculeValidator(enumerated_molecules, smiles=True)

        return enumerated_molecules

    def _load_amino_acids(self):

        """

        Load the amino acids.

        Returns:
            natural_amino_acids (List): The natural amino acid list represented in smiles

        Reference:
            amino_acids = {
                "Alanine": "C",
                "Arginine": "CCCCNC(N)=N",
                "Asparagine": "CCC(N)=O",
                "Aspartic Acid": "CC(O)=O",
                "Cysteine": "CS",
                "Glutamic Acid": "CCC(O)=O",
                "Glutamine": "CCC(N)=O",
                "Glycine": "[H]",
                "Histidine":"CC1=CNC=N1",
                "Isoleucine": "C(CC)([H])C",
                "Leucine": "CC(C)C",
                "Lysine": "CCCCN",
                "Methionine": "CCSC",
                "PhenylAlanine": "CC1=CC=CC=C1",
                "Proline": "",
                "Serine": "CO",
                "Threonine": "C(C)([H])O",
                "Tryptophan": "CCC1=CNC2=C1C=CC=C2",
                "Tyrosine": "CC1=CC=C(O)C=C1",
                "Valine": "C(C)C"
            }

        """

        natural_amino_acids = ["C", "CCCNC(N)=N", "CCC(N)=O", "CC(O)=O", "CS", "CCC(O)=O", "CCC(O)=O", "CCC(N)=O", "[H]",
                               "CC1=CNC=N1", "C(CC)([H])C", "CC(C)C", "CCCCN", "CCSC", "CC1=CC=CC=C1", "CO", "C(C)([H])O",
                               "CCC1=CNC2=C1C=CC=C2", "CC1=CC=C(O)C=C1", "C(C)C"]

        return natural_amino_acids

        

