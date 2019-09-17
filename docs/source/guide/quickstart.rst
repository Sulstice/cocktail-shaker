.. _gettingstarted:

Getting Started
===============

This page instructs you on how to get started with Cocktail Shaker. To begin, make sure you have
:ref:`installing cocktail_shaker <install>`.

Basic Usage
-----------

The simplest way to use Cocktail Shaker is to create a cocktail object with the ``shake`` function and create new compounds:

    >>> cocktail = Cocktail(['c1cc(CCCO)ccc1'])
    >>> new_compounds = cocktail.shake()
    [RDKit_Mol_Object, RDKit_Mol_Object, RDKit_Mol_Object...]

Write the new compounds to an SDF file:

    >>> FileWriter('new_compounds', new_compounds, 'sdf')

In this example, we have taken one smiles string and expanded the compounds into a variety of variations in one SDF file.
