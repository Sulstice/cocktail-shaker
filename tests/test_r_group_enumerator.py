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

def test_cocktail_shake():

    """

    Verify Cocktail Shaker is working extensively and also include the stress testing.

    """
    peptide_backbone = PeptideMolecule(1)
    cocktail = Cocktail(
        peptide_backbone,
        ligand_library = ['Br']
    )
    combinations = cocktail.shake()

    assert len(combinations) == 1
    assert combinations[0] == 'NC(Br)C(=O)NCC(=O)O'

    peptide_backbone = PeptideMolecule(2)
    cocktail = Cocktail(
        peptide_backbone,
        ligand_library = ['Br', 'I']
    )
    combinations = cocktail.shake()

    assert len(combinations) == 2
    assert 'NC(Br)C(=O)NC(I)C(=O)NCC(=O)O' in combinations
    assert 'NC(I)C(=O)NC(Br)C(=O)NCC(=O)O' in combinations

    peptide_backbone = PeptideMolecule(2)
    cocktail = Cocktail(
        peptide_backbone,
        ligand_library = ['Br', 'I', 'F', 'N', 'C', 'Cl']
    )
    combinations = cocktail.shake()
    assert len(combinations) == 30

    peptide_backbone = PeptideMolecule(6)
    cocktail = Cocktail(peptide_backbone,
                        ligand_library = ['Br', 'I', 'F', 'N', 'C', 'Cl']
                        )
    combinations = cocktail.shake()
    assert len(combinations) == 720

def test_cocktail_enumerate():

    """

    Test the enumerate function within Cocktail Class Object

    """

    peptide_backbone = PeptideMolecule(1)
    cocktail = Cocktail(
        peptide_backbone,
        ligand_library = ['Br']
    )
    cocktail.shake()
    combinations = cocktail.enumerate(enumeration_complexity='High')

    assert 'C(NCC(=O)O)(C(N)Br)=O' in combinations
    assert 'OC(CNC(C(Br)N)=O)=O' in combinations
    assert len(combinations) == 720

def test_cocktail_isomers():

    """

    Test the stereoisomer function within the Cocktail Class Object


    """

    peptide_backbone = PeptideMolecule(2)
    cocktail = Cocktail(
        peptide_backbone,
        ligand_library = ['I', 'Br'],
        enable_isomers = True
    )
    combinations = cocktail.shake()

    assert 'N[C@@H](Br)C(=O)N[C@@H](I)C(=O)NCC(=O)O' in combinations
    assert 'N[C@H](I)C(=O)N[C@H](Br)C(=O)NCC(=O)O' in combinations
    assert len(combinations) == 8

def test_amino_acids():

    """

    Test the amino acid function within the Cocktail Class Object

    """

    peptide_backbone = PeptideMolecule(2)
    cocktail = Cocktail(
        peptide_backbone,
        ligand_library = [],
        enable_isomers = True
    )
    combinations = cocktail.shake()

    assert len(combinations) == 18

    assert 'CSCCC(N)C(=O)NCC(=O)O' in combinations
    assert 'NC(Cc1ccccc1)C(=O)NCC(=O)O' in  combinations
    assert 'CC(N)C(=O)NCC(=O)O' in combinations
    assert 'CC(C)CC(N)C(=O)NCC(=O)O' in combinations
    assert 'CC(C)C(N)C(=O)NCC(=O)O' in combinations
    assert 'CCC(C)C(N)C(=O)NCC(=O)O' in combinations
    assert 'NC(CS)C(=O)NCC(=O)O' in combinations
    assert 'NC(Cc1ccc(O)cc1)C(=O)NCC(=O)O' in combinations
    assert 'NC(Cc1c[nH]cn1)C(=O)NCC(=O)O' in combinations
    assert 'N=C(N)NCCCC(N)C(=O)NCC(=O)O' in combinations
    assert 'NCCCCC(N)C(=O)NCC(=O)O' in combinations
    assert 'NC(CO)C(=O)NCC(=O)O' in combinations
    assert 'NC(=O)CCC(N)C(=O)NCC(=O)O' in combinations
    assert 'NC(CC(=O)O)C(=O)NCC(=O)O' in combinations
    assert 'NC(CCc1c[nH]c2ccccc12)C(=O)NCC(=O)O' in combinations
    assert 'NC(CCC(=O)O)C(=O)NCC(=O)O' in combinations
    assert '[H]C(N)C(=O)NCC(=O)O' in combinations
    assert 'CC(O)C(N)C(=O)NCC(=O)O' in combinations