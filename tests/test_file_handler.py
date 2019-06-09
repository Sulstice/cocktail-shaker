# -*- coding: utf-8 -*-
#
# Testing file content definitions
#
# ------------------------------------------------

# imports
# -------
import os
from rdkit import Chem
from packages.r_group_enumerator import RGroupMolObject
from packages.file_handler import FileParser, FileWriter


def test_encoding_checker():

    """

    This test will handle UTF8 encoding checker, for now we don't support international encoding.
    Reason being it gets to be tricky with new line endings especially in chinese files.

    """

    import os
    path = os.path.dirname(os.path.abspath(__file__))

    success_file_object = FileParser(path + "/resources/test.sdf")
    assert success_file_object.check_utf8() == True

    fail_file_object = FileParser(path + "/resources/test_fail_encoding.sdf")
    assert fail_file_object.check_utf8() == False

def test_file_production():

    """

    This Tests the File Parser object in the ability to produce files

    """

    scaffold_molecule = RGroupMolObject([Chem.MolFromSmiles('c1cc(CCCO)ccc1'), Chem.MolFromSmiles('c1cc(CCCBr)ccc1')])
    patterns_found = scaffold_molecule.find_r_groups()
    modified_molecules = scaffold_molecule.r_group_enumerator(patterns_found=patterns_found)
    FileWriter("tests/test", modified_molecules, "sdf")
    FileWriter("tests/test", modified_molecules, "txt")

    from pathlib import Path
    dir_path = os.path.dirname(os.path.realpath(__file__))

    assert Path(dir_path + "/test.sdf").is_file() == True
    assert Path(dir_path + "/test.txt").is_file() == True


def test_file_parser():

    """

    This tests the file parsing of sdf files, txts, strings for the FileWriter object.

    """

    cwd = os.getcwd()

    if not FileParser(os.path.join(cwd, "tests/resources", "test.sdf")):
        return False

    if not FileParser(os.path.join(cwd, "tests/resources", "test.mol2")):
        return False
