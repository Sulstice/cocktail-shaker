# -*- coding: utf-8 -*-
#
# Testing content definitions.
#
# ------------------------------------------------


# imports
# -------
import pytest, yaml
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

    """
    molecule_success = RGroupMolObject(Chem.MolFromSmiles("c1cc(CCCO)ccc1"))

    assert molecule_success.validate_smiles(smiles="c1cc(CCCO)ccc1") is True

def test_primary_finding_r_groups():

    """

    Test finding of the R Group

    """

    # Alcohols
    # --------

    primary_alcohol_molecule = ["c1cc(CCCO)ccc1", '[OX2H]']

    assert len(Chem.MolFromSmiles(primary_alcohol_molecule[0]).GetSubstructMatches(
        Chem.MolFromSmarts(primary_alcohol_molecule[1]), uniquify=False)) == 1

    primary_alcohol_molecule = RGroupMolObject(Chem.MolFromSmiles("c1cc(CCCO)ccc1"))
    patterns = primary_alcohol_molecule.find_r_groups()
    print ("THIS IS PATTERNS" + str(patterns))
    
    # Carboxylic Acids
    # ----------------

    # included is the acid and the conjugate base
    primary_carboxylic_acid = ['C1C(CCCC1)CCC(=O)O', '[CX3](=O)[OX1H0-,OX2H1]']
    assert len(Chem.MolFromSmiles(primary_carboxylic_acid[0]).GetSubstructMatches(
        Chem.MolFromSmarts(primary_carboxylic_acid[1]), uniquify=False)) == 1

    # Halogens
    # --------

    primary_bromine_molecule = ["c1cc(CCCBr)ccc1", '[Br]']
    primary_chlorine_molecule = ["c1cc(CCCCl)ccc1", '[Cl]']
    primary_fluorine_molecule = ["c1cc(CCCF)ccc1", '[F]']

    # Bromine Compounds
    assert len(Chem.MolFromSmiles(primary_bromine_molecule[0]).GetSubstructMatches(
        Chem.MolFromSmarts(primary_bromine_molecule[1]), uniquify=False)) == 1

    # Chlorine Compounds
    assert len(Chem.MolFromSmiles(primary_chlorine_molecule[0]).GetSubstructMatches(
        Chem.MolFromSmarts(primary_chlorine_molecule[1]), uniquify=False)) == 1

    # Fluorine
    assert len(Chem.MolFromSmiles(primary_fluorine_molecule[0]).GetSubstructMatches(
        Chem.MolFromSmarts(primary_fluorine_molecule[1]), uniquify=False)) == 1


