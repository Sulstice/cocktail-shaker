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

from fragmentation.editable_molecule import RGroupMolObject

def test_validation_smiles():

    molecule_success = "c1cc(CCCO)ccc1"
    molecule_fail = "c1234cc(CCCO)ccc1"

    assert len(validate_smiles(molecule_success)) == 0
    assert len(validate_smiles(molecule_fail)) >= 1