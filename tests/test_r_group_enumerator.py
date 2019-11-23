# -*- coding: utf-8 -*-
#
# Testing content definitions.
#
# ------------------------------------------------


# imports
# -------
import ruamel.yaml as yaml
from cocktail_shaker.functional_group_enumerator import Cocktail
from cocktail_shaker.peptide_builder import PeptideMolecule


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

def test_validation_smiles():

    """

    Validate that the validator package is working correctly,

    Come up with more a validation test for verification of molecules.

    """
    peptide_backbone = PeptideMolecule(1)
    cocktail = Cocktail(
        peptide_backbone,
        ligand_library = ['Br']
    )
    modified_molecules = cocktail.shake()
