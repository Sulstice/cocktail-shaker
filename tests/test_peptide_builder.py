# -*- coding: utf-8 -*-
#
# Testing content definitions.
#
# ------------------------------------------------

from cocktail_shaker.peptide_builder import PeptideMolecule

def test_peptide_molecule_backbone():

    """

    This test will serve that the Peptide molecule string builder accurately builds the backbone with the appropiate
    number of slots

    """

    import re

    # Peptide Molecule Test 1
    peptide_molecule = PeptideMolecule(1)
    peptide_molecule = str(peptide_molecule)
    peptide_backbone_length = int(max(map(int, re.findall(r"\d+", str(peptide_molecule)))))
    assert peptide_molecule == 'NC([*:1])C(NCC(O)=O)=O'
    assert peptide_backbone_length == 1

    # Peptide Molecule Test 2
    peptide_molecule = PeptideMolecule(6)
    peptide_molecule = str(peptide_molecule)
    peptide_backbone_length = int(max(map(int, re.findall(r"\d+", str(peptide_molecule)))))
    assert peptide_molecule == 'NC([*:1])C(NC([*:2])C(NC([*:3])C(NC([*:4])C(NC([*:5])C(NC([*:6])C(NCC(O)=O)=O)=O)=O)=O)=O)=O'
    assert peptide_backbone_length == 6

    # Peptide Molecule Test 3
    peptide_molecule = PeptideMolecule(20)
    peptide_molecule = str(peptide_molecule)
    peptide_backbone_length = int(max(map(int, re.findall(r"\d+", str(peptide_molecule)))))
    assert peptide_molecule == 'NC([*:1])C(NC([*:2])C(NC([*:3])C(NC([*:4])C(NC([*:5])C(NC([*:6])C(NC([*:7])C(NC([*:8])C(NC([*:9])C(NC([*:10])C(NC([*:11])C(NC([*:12])C(NC([*:13])C(NC([*:14])C(NC([*:15])C(NC([*:16])C(NC([*:17])C(NC([*:18])C(NC([*:19])C(NC([*:20])C(NCC(O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O'
    assert peptide_backbone_length == 20

