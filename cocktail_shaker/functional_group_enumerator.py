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
import functools

# Cocktail Shaker Imports
# ---
# --------------------
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
        self.ligand_library = ligand_library

        # Detect the proline amino acid on the N-terminus, set the max peptide length accordingly.
        if self.peptide_backbone[0:6] == 'N2CCCC2':
            self.peptide_backbone_length = int(max(map(int, re.findall(r"\d+", self.peptide_backbone[7:]))))

        self.peptide_backbone_length = int(max(map(int, re.findall(r"\d+", self.peptide_backbone))))
        self.enable_isomers = enable_isomers
        self.include_amino_acids = include_amino_acids

        # Experimental cache argument.
        self.cache = False

        self.combinations = []

        if (self.peptide_backbone_length > len(self.ligand_library) or not ligand_library) and self.include_amino_acids == False:
            print ("Cocktail Shaker Error: Peptide Backbone Length needs to be less than or equal to for your library")
            raise IndexError

    def shake(self, store_as_pickle = False, compound_filters = []):

        """

        Generate all combinations of a molecule

        Arguments:
            store_as_pickle (Bool): If the user wants to store the results into a pickle file.
            compound_filters (List): List of filters that the user would like to apply to their library.

        Returns:
            results (List): List of compounds generated into SMILES format.

        """

        results = []
        peptide = self.peptide_backbone

        if self.include_amino_acids:
            if not self.ligand_library:
                self.ligand_library = self._load_amino_acids()
            else:
                self.ligand_library.extend(self._load_amino_acids())

        combinations = list(itertools.permutations(self.ligand_library, self.peptide_backbone_length))

        print ("Generating %s Compounds..." % str(len(combinations)))

        for i in (range(len(combinations))):
            combination = list(combinations[i])
            peptide_molecule = str(peptide)
            for j in range(0, len(combination)):

                modified_molecule = Chem.ReplaceSubstructs(Chem.MolFromSmiles(peptide_molecule),
                                                               Chem.MolFromSmiles('[*:'+ str(j+1)+']'),
                                                               Chem.MolFromSmiles(str(combination[j])))

                peptide_molecule = Chem.MolToSmiles(modified_molecule[0], isomericSmiles=True, canonical=True)

            # Enable StereoChemistry
            if self.enable_isomers:

                from rdkit.Chem import EnumerateStereoisomers
                molecule = Chem.MolFromSmiles(str(peptide_molecule))
                options = EnumerateStereoisomers.StereoEnumerationOptions(unique=True, tryEmbedding=True)
                isomers = tuple(EnumerateStereoisomers.EnumerateStereoisomers(molecule, options=options))
                for smile in sorted(Chem.MolToSmiles(isomer, isomericSmiles=True,  canonical=True) for isomer in isomers):
                    results.append(str(smile))
            else:
                results.append(str(peptide_molecule))

        # Remove Duplicates
        results = list(set(results))

        # Validate Smiles
        MoleculeValidator(results, smiles=True)

        # Store in combinations if enumeration
        self.combinations = results

        # Handle storage into a pickle file.
        if store_as_pickle:

            import pickle
            with open("cocktail_shaker.pickle", 'wb') as cache_handle:
                print("saving result to cache '%s'" % "cocktail_shaker.pickle")
                pickle.dump(results, cache_handle)

        if compound_filters:
            results = self._apply_drug_filters(compound_filters)
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

        print ("Enumerating %s Compounds...." % len(self.combinations))

        if enumeration_complexity.lower() == 'low':
            complexity = 100
        elif enumeration_complexity.lower() == 'medium':
            complexity = 1000
        elif enumeration_complexity.lower() == 'high':
            complexity = 10000
        else:
            complexity = 10

        enumerated_molecules = []
        for i in range(len(self.combinations)):
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
        if dimensionality == '2D':
            MoleculeValidator(enumerated_molecules, smiles=False)
        else:
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

    def _apply_drug_filters(self, filters):

        """

        Applies a drug filter or more and returns the modified compound library.

        Arguments:
            filters (List): List of filters that you would like to apply.

        Reference:
            filters = [
                        "Lipinski",
                        "Ghose",
                        "Veber",
                        "Rule of 3",
                        "REOS",
                        "Drug-like",
                        "All",
            ]

        """

        from rdkit.Chem import Descriptors

        print ("Applying Drug Filters...")

        # Store our final list
        results = []

        # Boolean for filters
        lipinski = False
        rule_of_3 = False
        ghose_filter = False
        veber_filter = False
        reos_filter = False
        drug_like_filter = False
        apply_all_filters = False

        # Determine which filters to apply
        if "Lipinski" in filters:
            lipinski = True
        if "Rule of 3" in filters:
            rule_of_3 = True
        if "Ghose" in filters:
            ghose_filter = True
        if "Veber" in filters:
            veber_filter = True
        if "REOS" in filters:
            reos_filter = True
        if "Drug-Like"in filters:
            drug_like_filter = True
        if "All" in filters or lipinski and reos_filter and rule_of_3 and ghose_filter and veber_filter and drug_like_filter:
            apply_all_filters = True
            lipinski = True
            rule_of_3 = True
            ghose_filter = True
            veber_filter = True
            reos_filter = True
            drug_like_filter = True

        for i in range(len(self.combinations)):

            lipinski_pass = False
            ghose_pass = False
            veber_pass = False
            rule_of_3_pass = False
            reos_pass = False
            drug_like_pass = False


            molecule_smiles = self.combinations[i]
            molecule = Chem.MolFromSmiles(molecule_smiles)

            # Default Descriptors
            molecular_weight = Descriptors.ExactMolWt(molecule)
            logp = Descriptors.MolLogP(molecule)
            h_bond_donor = Descriptors.NumHDonors(molecule)
            h_bond_acceptors = Descriptors.NumHAcceptors(molecule)
            rotatable_bonds = Descriptors.NumRotatableBonds(molecule)
            number_of_atoms = Chem.rdchem.Mol.GetNumAtoms(molecule)
            molar_refractivity = Chem.Crippen.MolMR(molecule)
            topological_surface_area_mapping = Chem.QED.properties(molecule).PSA
            formal_charge = Chem.rdmolops.GetFormalCharge(molecule)
            heavy_atoms = Chem.rdchem.Mol.GetNumHeavyAtoms(molecule)
            num_of_rings = Chem.rdMolDescriptors.CalcNumRings(molecule)

            # Lipinski
            if lipinski:
                if molecular_weight <= 500 and logp <= 5 and h_bond_donor <= 5 and h_bond_acceptors <= 5 and rotatable_bonds <= 5:
                    if apply_all_filters:
                        lipinski_pass = True
                    else:
                        results.append(molecule_smiles)

            # Ghose Filter
            if ghose_filter:
                if molecular_weight >= 160 and molecular_weight <= 480 and logp >= 0.4 and logp <= 5.6 and number_of_atoms >= 20 and number_of_atoms <= 70 and molar_refractivity >= 40 and molar_refractivity <= 130:
                    if apply_all_filters:
                        ghose_pass = True
                    else:
                        results.append(molecule_smiles)

            # Veber Filter
            if veber_filter:
                if rotatable_bonds <= 10 and topological_surface_area_mapping <= 140:
                    if apply_all_filters:
                        veber_pass = True
                    else:
                        results.append(molecule_smiles)

            # Rule of 3
            if rule_of_3:
                if molecular_weight <= 300 and logp <= 3 and h_bond_donor <= 3 and h_bond_acceptors <= 3 and rotatable_bonds <= 3:
                    if apply_all_filters:
                        rule_of_3_passs = True
                    else:
                        results.append(molecule_smiles)

            # REOS Filter
            if reos_filter:
                if molecular_weight >= 200 and molecular_weight <= 500 and logp >= int(0 - 5) and logp <= 5 and h_bond_donor >= 0 and h_bond_donor <= 5 and h_bond_acceptors >= 0 and h_bond_acceptors <= 10 and formal_charge >= int(0-2) and formal_charge <= 2 and rotatable_bonds >= 0 and rotatable_bonds <= 8 and heavy_atoms >= 15 and heavy_atoms <= 50:
                    if apply_all_filters:
                        reos_pass = True
                    else:
                        results.append(molecule_smiles)

            #Drug Like Filter
            if drug_like_filter:
                if molecular_weight < 400 and num_of_rings > 0 and rotatable_bonds < 5 and h_bond_donor <= 5 and h_bond_acceptors <= 10 and logp < 5:
                    if apply_all_filters:
                        drug_like_pass = True
                    else:
                        results.append(molecule_smiles)

            if lipinski_pass and ghose_pass and veber_pass and rule_of_3_pass and reos_pass and drug_like_pass:
                results.append(molecule_smiles)

        return results