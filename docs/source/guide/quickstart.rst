.. _gettingstarted:

Getting Started
===============

This page instructs you on how to get started with cocktail-shaker. To begin, make sure you have performed
:ref:`installing cocktail_shaker <install>`.

Basic Usage
-----------

The simplest way to use cocktail-shaker is to create a cocktail object with the ``shake`` function and create new compounds:

    >>> from cocktail_shaker import Cocktail
    >>> from cocktail_shaker import PeptideMolecule
    >>> peptide_backbone = PeptideMolecule(2)
    >>> cocktail = Cocktail(peptide_backbone,ligand_library = ['Br', 'I'])
    >>> combinations = cocktail.shake()
    >>> print (combinations)
    >>> ['NC(Br)C(=O)NC(I)C(=O)NCC(=O)O', 'NC(I)C(=O)NC(Br)C(=O)NCC(=O)O']

Write the new compounds into an SDF file:

    >>> FileWriter('new_compounds', new_compounds, 'sdf')

In this example, we have taken one SMILES string and expanded the compounds into a variety of variations into one SDF file.
