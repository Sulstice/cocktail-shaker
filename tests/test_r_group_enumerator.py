# -*- coding: utf-8 -*-
#
# Testing content definitions.
#
# ------------------------------------------------


# imports
# -------
import yaml
from rdkit import Chem
from packages.r_group_enumerator import RGroupMolObject

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

    TODO: Key thing to note molecule validator passes strings that are incoherently incorrect when rendered into a
    molecule object.

    Come up with more a validation test for verification of molecules.

    """
    molecule_success = RGroupMolObject(Chem.MolFromSmiles("c1cc(CCCO)ccc1"))

    assert molecule_success.validate_smiles(smiles="c1cc(CCCO)ccc1") is True

def test_primary_finding_r_groups():

    """

    This test is associated with the RGroupMolObject being able to detect R Groups
    dependent on the datasource.

    If any test fails, it means we cannot support that R Group currently.

    TODO: Handle more sophisticated molecules and run tests on more complex strucutres.

    """

    # Alcohols
    # --------
    # primary_alcohol_molecule = ["c1cc(CCCO)ccc1", '[OX2H]']
    primary_alcohol_molecule = RGroupMolObject(Chem.MolFromSmiles("c1cc(CCCO)ccc1"))
    patterns = primary_alcohol_molecule.find_r_groups()
    assert patterns['Primary Alcohol'][0] == 'O'
    assert patterns['Primary Alcohol'][1] == '[OX2H]'

    # Carboxylic Acids
    # ----------------

    # included is the acid and the conjugate base
    # primary_carboxylic_acid = ['C1C(CCCC1)CCC(=O)O', '[CX3](=O)[OX1H0-,OX2H1]']
    primary_carboxylic_acid = RGroupMolObject(Chem.MolFromSmiles("C1C(CCCC1)CCC(=O)O"))
    patterns = primary_carboxylic_acid.find_r_groups()
    assert patterns['Carboxylic Acid'][0] == 'C(=O)O'
    assert patterns['Carboxylic Acid'][1] == '[CX3](=O)[OX1H0-,OX2H1]'

    # # Halogens
    # # --------
    # primary_bromine_molecule = ["c1cc(CCCBr)ccc1", '[Br]']
    # primary_chlorine_molecule = ["c1cc(CCCCl)ccc1", '[Cl]']
    # primary_fluorine_molecule = ["c1cc(CCCF)ccc1", '[F]']

    # Bromine Compounds
    primary_bromine_molecule = RGroupMolObject(Chem.MolFromSmiles("c1cc(CCCBr)ccc1"))
    patterns = primary_bromine_molecule.find_r_groups()
    assert patterns['Bromine'][0] == 'Br'
    assert patterns['Bromine'][1] == '[Br]'

    # Chlorine Compounds
    primary_chlorine_molecule = RGroupMolObject(Chem.MolFromSmiles("c1cc(CCCCl)ccc1"))
    patterns = primary_chlorine_molecule.find_r_groups()
    assert patterns['Chlorine'][0] == 'Cl'
    assert patterns['Chlorine'][1] == '[Cl]'

    # Fluorine
    primary_fluorine_molecule = RGroupMolObject(Chem.MolFromSmiles("c1cc(CCCF)ccc1"))
    patterns = primary_fluorine_molecule.find_r_groups()
    assert patterns['Fluorine'][0] == 'F'
    assert patterns['Fluorine'][1] == '[F]'

