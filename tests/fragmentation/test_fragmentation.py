# -*- coding: utf-8 -*-
#
# Testing content definitions.
#
# ------------------------------------------------


# imports
# -------
import pytest
from rdkit import Chem
from molvs import validate_smiles

def test_validation_smiles():

    """

    validate that the validator package is working correctly

    """
    molecule_success = "c1cc(CCCO)ccc1"
    molecule_fail = "c1234cc(CCCO)ccc1"

    assert len(validate_smiles(molecule_success)) == 0
    assert len(validate_smiles(molecule_fail)) >= 1

def test_primary_finding_r_groups():

    """

    Test finding of the R Group

    """

    primary_alcohol_molecule = ["c1cc(CCCO)ccc1", '[OX2H]']
    primary_bromine_molecule = ["c1cc(CCCBr)ccc1", '[Br]']
    primary_chlorine_molecule = ["c1cc(CCCCl)ccc1", '[Cl]']
    primary_fluorine_molecule = ["c1cc(CCCF)ccc1", '[F]']

    # Alcohols
    # --------
    assert len(Chem.MolFromSmiles(primary_alcohol_molecule[0]).GetSubstructMatches(
               Chem.MolFromSmarts(primary_alcohol_molecule[1]), uniquify=False)) == 1

    # Halogens
    # --------

    # Bromine Compounds
    assert len(Chem.MolFromSmiles(primary_bromine_molecule[0]).GetSubstructMatches(
               Chem.MolFromSmarts(primary_bromine_molecule[1]), uniquify=False)) == 1

    # Chlorine Compounds
    assert len(Chem.MolFromSmiles(primary_chlorine_molecule[0]).GetSubstructMatches(
               Chem.MolFromSmarts(primary_chlorine_molecule[1]), uniquify=False)) == 1

    # Fluorine
    assert len(Chem.MolFromSmiles(primary_fluorine_molecule[0]).GetSubstructMatches(
               Chem.MolFromSmarts(primary_fluorine_molecule[1]), uniquify=False)) == 1


