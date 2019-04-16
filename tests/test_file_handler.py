# -*- coding: utf-8 -*-
#
# Testing file content definitions
#
# ------------------------------------------------

# imports
# -------
from packages.file_handler import FileParser

def test_encoding_checker():

    import os
    path = os.path.dirname(os.path.abspath(__file__))

    success_file_object = FileParser(path + "/resources/test.sdf")
    assert success_file_object.check_utf8() == True

    fail_file_object = FileParser(path + "/resources/test_fail_encoding.sdf")
    assert fail_file_object.check_utf8() == False
