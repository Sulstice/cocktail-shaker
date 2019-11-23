# -*- coding: utf-8 -*-
#
# Testing file content definitions
#
# ------------------------------------------------

# imports
# -------
import os
from cocktail_shaker.functional_group_enumerator import Cocktail
from cocktail_shaker.file_handler import FileParser, FileWriter
from cocktail_shaker.peptide_builder import PeptideMolecule

def test_encoding_checker():

    """

    This test will handle UTF8 encoding checker, for now we don't support international encoding.
    Reason being it gets to be tricky with new line endings especially in chinese files.

    """

    import os
    path = os.path.dirname(os.path.abspath(__file__))

    success_file_object = FileParser(path + "/resources/test.sdf")
    assert success_file_object.check_utf8() == True

    print ("Encoding Check Passes")

def test_file_production():

    """

    This Tests the File Parser object in the ability to produce files

    """

    peptide_backbone = PeptideMolecule(1)
    cocktail = Cocktail(
        peptide_backbone,
        ligand_library = ['Br']
    )
    modified_molecules = cocktail.shake()
    FileWriter("tests/test", modified_molecules, "sdf")
    FileWriter("tests/test", modified_molecules, "alc")
    FileWriter("tests/test", modified_molecules, "cdxml")
    FileWriter("tests/test", modified_molecules, "cerius")
    FileWriter("tests/test", modified_molecules, "charmm")
    FileWriter("tests/test", modified_molecules, "cif")
    FileWriter("tests/test", modified_molecules, "cml")
    FileWriter("tests/test", modified_molecules, "gjf")
    FileWriter("tests/test", modified_molecules, "gromacs")
    FileWriter("tests/test", modified_molecules, "hyperchem")
    FileWriter("tests/test", modified_molecules, "jme")
    FileWriter("tests/test", modified_molecules, "mol")
    FileWriter("tests/test", modified_molecules, "mol2")
    FileWriter("tests/test", modified_molecules, "pdb")
    FileWriter("tests/test", modified_molecules, "mrv")
    FileWriter("tests/test", modified_molecules, "sdf3000")
    FileWriter("tests/test", modified_molecules, "sln")
    FileWriter("tests/test", modified_molecules, "xyz")

    print ("File Writing Passes")
    from pathlib import Path
    dir_path = os.path.dirname(os.path.realpath(__file__))
    #
    assert Path(dir_path + "/test.sdf").is_file() == True
    assert Path(dir_path + "/test.alc").is_file() == True
    assert Path(dir_path + "/test.cdxml").is_file() == True
    assert Path(dir_path + "/test.cerius").is_file() == True
    assert Path(dir_path + "/test.charmm").is_file() == True
    assert Path(dir_path + "/test.cif").is_file() == True
    assert Path(dir_path + "/test.cml").is_file() == True
    assert Path(dir_path + "/test.gjf").is_file() == True
    assert Path(dir_path + "/test.gromacs").is_file() == True
    assert Path(dir_path + "/test.hyperchem").is_file() == True
    assert Path(dir_path + "/test.jme").is_file() == True
    assert Path(dir_path + "/test.mol").is_file() == True
    assert Path(dir_path + "/test.mol2").is_file() == True
    assert Path(dir_path + "/test.pdb").is_file() == True
    assert Path(dir_path + "/test.mrv").is_file() == True
    assert Path(dir_path + "/test.sdf3000").is_file() == True
    assert Path(dir_path + "/test.sln").is_file() == True
    assert Path(dir_path + "/test.xyz").is_file() == True

    print ("File writing passes")

def test_file_parser():

    """

    This tests the file parsing of sdf files, txts, strings for the FileWriter object.

    """

    cwd = os.getcwd()

    if not FileParser(os.path.join(cwd, "tests/resources", "test.sdf")):
        return False

    print ("File Parser passes")