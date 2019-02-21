# -*- coding: utf-8 -*-
#
# Testing content definitions.
#
# ------------------------------------------------


# imports
# -------
import pytest
from sanity_check.EncodingChecker import LigandFile

def test_encoding_checker():

    success_file_object = LigandFile('resources/CHEMBL-001.mol')
    assert success_file_object.check_utf8() == True

    fail_file_object = LigandFile('resources/CHEMBL-Fail-001.mol')
    assert fail_file_object.check_utf8() == False

