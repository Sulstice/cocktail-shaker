# -*- coding: utf-8 -*-
#
# Testing content definitions.
#
# ------------------------------------------------


# imports
# -------
import ruamel.yaml as yaml
from cocktail_shaker.functional_group_enumerator import Cocktail
from cocktail_shaker.peptide_builder import PeptideBuilder


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
    peptide_backbone = PeptideBuilder(1)
    cocktail = Cocktail(
        peptide_backbone,
        ligand_library = ['Br']
    )
    combinations = cocktail.shake()

    assert len(combinations) == 1
    assert combinations[0] == 'NC(Br)C(=O)NCC(=O)O'

    peptide_backbone = PeptideBuilder(2)
    cocktail = Cocktail(
        peptide_backbone,
        ligand_library = ['Br', 'I']
    )
    combinations = cocktail.shake()

    assert len(combinations) == 2
    assert 'NC(Br)C(=O)NC(I)C(=O)NCC(=O)O' in combinations
    assert 'NC(I)C(=O)NC(Br)C(=O)NCC(=O)O' in combinations

    peptide_backbone = PeptideBuilder(2)
    cocktail = Cocktail(
        peptide_backbone,
        ligand_library = ['Br', 'I', 'F', 'N', 'C', 'Cl']
    )
    combinations = cocktail.shake()
    assert len(combinations) == 30

    peptide_backbone = PeptideBuilder(6)
    cocktail = Cocktail(peptide_backbone,
                        ligand_library = ['Br', 'I', 'F', 'N', 'C', 'Cl']
                        )
    combinations = cocktail.shake()
    assert len(combinations) == 720

def test_cocktail_enumerate():

    """

    Test the enumerate function within Cocktail Class Object

    """

    peptide_backbone = PeptideBuilder(1)
    cocktail = Cocktail(
        peptide_backbone,
        ligand_library = ['Br']
    )
    cocktail.shake()
    combinations = cocktail.enumerate(enumeration_complexity='High')

    assert 'C(NCC(=O)O)(C(N)Br)=O' in combinations
    assert 'OC(CNC(C(Br)N)=O)=O' in combinations

def test_cocktail_isomers():

    """

    Test the stereoisomer function within the Cocktail Class Object


    """

    peptide_backbone = PeptideBuilder(2)
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

    peptide_backbone = PeptideBuilder(1)
    cocktail = Cocktail(
        peptide_backbone,
        ligand_library = [],
        include_amino_acids = True
    )
    combinations = cocktail.shake()

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

def test_circular_peptide_production():

    """

    Test the Circular peptide production

    """
    peptide_backbone = PeptideBuilder(
        length_of_peptide = 4,
        circular = True,
    )

    # assert 'O=C1C([*:1])NC(=O)C([*:2])NC(=O)C([*:3])NC(=O)C([*:4])N1' == peptide_backbone


def test_drug_filters():

    """

    Test the drug filters

    """

    peptide_backbone = PeptideBuilder(
        length_of_peptide = 4,
        circular = True,
    )
    #
    # cocktail = Cocktail(peptide_backbone, ligand_library = ["Br", "Cl", "I", "F"])
    #
    # combinations = cocktail.shake(compound_filters=["Lipinski"])
    # assert 'O=C1NC(Br)C(=O)NC(I)C(=O)NC(Cl)C(=O)NC1F' in combinations
    #
    # combinations = cocktail.shake(compound_filters=["Ghose"])
    # assert len(combinations) == 0
    #
    # combinations = cocktail.shake(compound_filters=["Veber"])
    # assert 'O=C1NC(I)C(=O)NC(Br)C(=O)NC1F' in combinations
    #
    # combinations = cocktail.shake(compound_filters=["Rule of 3"])
    # assert len(combinations) == 0
    #
    # combinations = cocktail.shake(compound_filters=["REOS"])
    # assert 'O=C1NC(Br)C(=O)NC(I)C(=O)NC1Cl' in combinations
    #
    # combinations = cocktail.shake(compound_filters=["Drug-like"])
    # assert len(combinations) == 0
