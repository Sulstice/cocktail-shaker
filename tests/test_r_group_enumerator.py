# -*- coding: utf-8 -*-
#
# Testing content definitions.
#
# ------------------------------------------------


# imports
# -------
import yaml
from rdkit import Chem
from cocktail_shaker.functional_group_enumerator import Cocktail

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

    Come up with more a validation test for verification of molecules.

    """
    molecule_success = Cocktail(Chem.MolFromSmiles("c1cc(CCCO)ccc1"))

    assert molecule_success.validate_smiles(smiles="c1cc(CCCO)ccc1") is True

def test_primary_finding_r_groups():

    """

    This test is associated with the Cocktail being able to detect R Groups
    dependent on the datasource.

    If any test fails, it means we cannot support that R Group currently.

    TODO: Handle more sophisticated molecules and run tests on more complex strucutres.

    """

    # Alcohols
    # --------
    primary_alcohol_molecule = Cocktail(Chem.MolFromSmiles("c1cc(CCCO)ccc1"))
    patterns = primary_alcohol_molecule.detect_functional_groups()
    assert patterns['Primary Alcohol'][0] == 'O'
    assert patterns['Primary Alcohol'][1] == '[OX2H]'

    # Aldehydes
    # ---------
    # aldehyde_molecule = Cocktail(Chem.MolFromSmiles('C1=CC=C(C=C1)C=O'))
    # patterns = aldehyde_molecule.detect_functional_groups()
    # assert patterns['Aldehyde'][0] == 'C(=O)H'
    # assert patterns['Aldehyde'][1] == '[CX3H1](=O)[#6]'

    # Acetic Anhydride
    # ----------------

    acetic_anhydride_molecule = Cocktail(Chem.MolFromSmiles('c1cc(CCCCC(=O)OC(=O)C)ccc1'))
    patterns = acetic_anhydride_molecule.detect_functional_groups()
    assert patterns['Acetic Anhydride'][0] == 'CC(=O)OC(=O)C'
    assert patterns['Acetic Anhydride'][1] == '[CX3](=[OX1])[OX2][CX3](=[OX1])'

    # Acyl Halides
    # ------------
    acyl_chloride_molecule = Cocktail(Chem.MolFromSmiles('c1cc(CCC(C(=O)Cl))ccc1'))
    patterns = acyl_chloride_molecule.detect_functional_groups()
    assert patterns['Acyl Chloride'][0] == 'C(=O)Cl'
    assert patterns['Acyl Chloride'][1] == '[CX3](=[OX1])[Cl]'

    acyl_bromide_molecule = Cocktail(Chem.MolFromSmiles('c1cc(CCC(C(=O)Br))ccc1'))
    patterns = acyl_bromide_molecule.detect_functional_groups()
    assert patterns['Acyl Bromide'][0] == 'C(=O)Br'
    assert patterns['Acyl Bromide'][1] == '[CX3](=[OX1])[Br]'

    acyl_fluoride_molecule = Cocktail(Chem.MolFromSmiles('c1cc(CCC(C(=O)F))ccc1'))
    patterns = acyl_fluoride_molecule.detect_functional_groups()
    assert patterns['Acyl Fluoride'][0] == 'C(=O)F'
    assert patterns['Acyl Fluoride'][1] == '[CX3](=[OX1])[F]'

    acyl_iodide_molecule = Cocktail(Chem.MolFromSmiles('c1cc(CCC(C(=O)I))ccc1'))
    patterns = acyl_iodide_molecule.detect_functional_groups()
    assert patterns['Acyl Iodide'][0] == 'C(=O)I'
    assert patterns['Acyl Iodide'][1] == '[CX3](=[OX1])[I]'

    # Amides
    # ------
    amide_molecule = Cocktail(Chem.MolFromSmiles('C1=CC=C(C=C1)C(=O)N'))
    patterns = amide_molecule.detect_functional_groups()
    assert patterns['Amide'][0] == 'C(=O)N'
    assert patterns['Amide'][1] == '[NX3][CX3](=[OX1])[#6]'

    # Ketones
    # -------
    ketone_molecule = Cocktail(Chem.MolFromSmiles('c1cc(CCCC(=O)OC)ccc1'))
    patterns = ketone_molecule.detect_functional_groups()
    assert patterns['Ketone'][0] == 'C(=O)OC'
    assert patterns['Ketone'][1] == '[CX3]=[OX1]'

    # Carboxylic Acids
    # ----------------
    # included is the acid and the conjugate base
    primary_carboxylic_acid = Cocktail(Chem.MolFromSmiles("C1C(CCCC1)CCC(=O)O"))
    patterns = primary_carboxylic_acid.detect_functional_groups()
    assert patterns['Carboxylic Acid'][0] == 'C(=O)O'
    assert patterns['Carboxylic Acid'][1] == '[CX3](=O)[OX1H0-,OX2H1]'

    # # Halogens
    # # --------

    # Bromine Compounds
    primary_bromine_molecule = Cocktail(Chem.MolFromSmiles("c1cc(CCCBr)ccc1"))
    patterns = primary_bromine_molecule.detect_functional_groups()
    assert patterns['Bromine'][0] == 'Br'
    assert patterns['Bromine'][1] == '[Br]'

    # Chlorine Compounds
    primary_chlorine_molecule = Cocktail(Chem.MolFromSmiles("c1cc(CCCCl)ccc1"))
    patterns = primary_chlorine_molecule.detect_functional_groups()
    assert patterns['Chlorine'][0] == 'Cl'
    assert patterns['Chlorine'][1] == '[Cl]'

    # Fluorine Compounds
    primary_fluorine_molecule = Cocktail(Chem.MolFromSmiles("c1cc(CCCF)ccc1"))
    patterns = primary_fluorine_molecule.detect_functional_groups()
    assert patterns['Fluorine'][0] == 'F'
    assert patterns['Fluorine'][1] == '[F]'

    # Nitro Compounds
    nitro_molecule = Cocktail(Chem.MolFromSmiles('C1=CC=C(C=C1)[N+](=O)[O-]'))
    patterns = nitro_molecule.detect_functional_groups()
    assert patterns['Nitro'][0] == '[N+](=O)[O-]'
    assert patterns['Nitro'][1] == '[$([NX3](=O)=O),$([NX3+](=O)[O-])][!#8]'

    # Ether Compounds
    ether_molecule = Cocktail(Chem.MolFromSmiles('c1cc(CCCCOC)ccc1'))
    patterns = ether_molecule.detect_functional_groups()
    assert patterns['Ether'][0] == 'COC'
    assert patterns['Ether'][1] == '[OD2]([#6])[#6]'

    # Sulfoxide
    # sulfoxide_molecule = Cocktail(Chem.MolFromSmiles('c1cc(CCCS(=O)(=O)C)ccc1'))
    # patterns = sulfoxide_molecule.detect_functional_groups()
    # assert patterns['Sulfoxide'][0] == 'S(=O)(=O)C'
    # assert patterns['Sulfoxide'][1] == '[$([#16X3]=[OX1]),$([#16X3+][OX1-])]'

    # Azides
    # azide_molecule = Cocktail(Chem.MolFromSmiles('c1cc(CCC([N-][N+]#N))ccc1'))
    # patterns = azide_molecule.detect_functional_groups()
    # assert patterns['Azide'][0] == '([N-][N+]#N)'
    # assert patterns['Azide'][1] == '[$(*-[NX2-]-[NX2+]#[NX1]),$(*-[NX2]=[NX2+]=[NX1-])]'

    # Amines
    primary_amine_molecule = Cocktail(Chem.MolFromSmiles('c1cc(CCCN)ccc1'))
    patterns = primary_amine_molecule.detect_functional_groups()
    assert patterns['Primary Amine'][0] == 'N'
    assert patterns['Primary Amine'][1] == '[NX3;H2;!$(NC=[!#6]);!$(NC#[!#6])][#6]'

    secondary_amine_molecule = Cocktail(Chem.MolFromSmiles('c1cc(CCCNC)ccc1'))
    patterns = secondary_amine_molecule.detect_functional_groups()
    assert patterns['Secondary Amine'][0] == 'NC'
    assert patterns['Secondary Amine'][1] == '[NX3;H2,H1;!$(NC=O)]'

    enamine_amine_molecule = Cocktail(Chem.MolFromSmiles('CC(=CC(=O)OC)N'))
    patterns = enamine_amine_molecule.detect_functional_groups()
    assert patterns['Enamine'][0] == 'N'
    assert patterns['Enamine'][1] == '[NX3][CX3]=[CX3]'

