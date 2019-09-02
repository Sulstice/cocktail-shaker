.. _gettingstarted:

Getting started
===============

This page gives a introduction on how to get started with cocktail-shaker. Before we start, make sure you have
:ref:`installing cocktail_shaker <install>`.

Basic usage
-----------

The simplest way to use Cocktail Shaker is to create a Cocktail object and with the ``shake`` function and create new compounds::

    >>> cocktail = Cocktail(['c1cc(CCCO)ccc1'])
    >>> new_compounds = cocktail.shake()
    [RDKit_Mol_Object, RDKit_Mol_Object, RDKit_Mol_Object...]

Write the new compounds to an SDF file.

    >>> FileWriter('new_compounds', new_compounds, 'sdf')

Here we have taken one smiles string and expanded the compounds to a variety of variations into one SDF file.
.


