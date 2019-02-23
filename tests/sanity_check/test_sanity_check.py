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

    import sys, os
    path = os.path.dirname(os.path.abspath(__file__))

    success_file_object = LigandFile(path + "/resources/CHEMBL-001.mol")
    assert success_file_object.check_utf8() == True

    fail_file_object = LigandFile(path + "/resources/CHEMBL-Fail-001.mol")
    assert fail_file_object.check_utf8() == False

