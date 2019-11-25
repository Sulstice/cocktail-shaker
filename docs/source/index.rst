.. Cocktail Shaker documentation master file, created by
   sphinx-quickstart on Wed Aug 28 23:41:51 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to Cocktail Shaker's documentation!
===========================================

.. toctree::
   :maxdepth: 2
   :caption: Contents:

.. Cocktail Shaker documentation master file, created by sphinx-quickstart on Tue Mar 24 16:12:38 2015.

Cocktail Shaker
===============

.. sectionauthor:: Suliman Sharif <sharifsuliman1@gmail.com>

Cocktail Shaker is a **high-performance drug enumeration and expansion library**.
Cocktail Shaker leverages the computational power of  **RDKit** to create and enumerate large volumes of drug compounds.

>>> from cocktail_shaker import Cocktail, FileWriter
>>> from cocktail_shaker import PeptideMolecule
>>> peptide_backbone = PeptideMolecule(2)
>>> cocktail = Cocktail(peptide_backbone,ligand_library = ['Br', 'I'])
>>> combinations = cocktail.shake()
>>> FileWriter('example', combinations, 'mol2')

Cocktail Shaker makes your drug enumeration and expansion life easy. It also generates your files for you in as many
formats needed for any cheminformatics software.

Features
--------


- File parsing of TXT, SDF, and Chemical Smiles.
- File writing in a variety of formats some of which include: cif, sdf, pdb, mol, mol2 and many others
- Ability to recognize and expand libraries of compounds some of which include halogens, acyl halides, aldehydes.
- Ability to enumerate in 1D, and 2D structures and produce those compounds.
- Supports Python versions 3.3+.
- Released under the `MIT license`_.

User guide
----------

A step-by-step guide to getting started with Cocktail Shaker.

.. toctree::
    :maxdepth: 2

    guide/installation
    guide/quickstart

API documentation
-----------------

Comprehensive API documentation with information on every function, class and method.

.. toctree::
    :maxdepth: 2

    guide/file_formats
    guide/file_handling
    guide/functional_groups
    guide/amino_acids
    guide/cocktail
    guide/contributing

Cocktail Shaker's license
-------------------------

Cocktail Shaker is released under the Mozilla Public License 2.0. This is a short, permissive software license that allows commercial use,
modifications, distribution, sublicensing and private use. Basically, you can do whatever you want with Cocktail Shaker as long as
you include the original copyright and license in any copies or derivative projects.

.. _`MIT license`: https://github.com/Sulstice/Cocktail-Shaker/blob/master/LICENSE
