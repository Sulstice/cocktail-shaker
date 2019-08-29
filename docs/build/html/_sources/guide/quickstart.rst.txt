.. _gettingstarted:

Getting started
===============

This page gives a introduction on how to get started with CIRpy. Before we start, make sure you have
:ref:`installed cocktail_shaker <install>`.

Basic usage
-----------

The simplest way to use Cocktail Shaker is to create a Cocktail object and with the ``shake`` function create new compounds!::

    >>> cocktail = Cocktail()
.

File Formats
------------

Input can additionally be parsed in a variety of file formats that are autodetected::

    >>> FileParser('name_of_file')

    For Example:

    >>> FileParser('targets.sdf')

    ['RDKit MolObject', 'RDKit MolObject']

The full list of file formats::

    txt       # Smiles Txt Format
    sdf       # 2D Structure file format
    mol2      # (Coming soon)


Output can additionally be returned in a variety of file formats that are specified your enumerated compounds and your
desired file extension::

    >>> FileWriter('name_of_file', compounds, extension of the file)

    For Example:

    >>> FileWriter('test'', compounds, '.mol2')



The full list of file formats::

    alc         # Alchemy format
    cdxml       # CambridgeSoft ChemDraw XML format
    cerius      # MSI Cerius II format
    charmm      # Chemistry at HARvard Macromolecular Mechanics file format
    cif         # Crystallographic Information File
    cml         # Chemical Markup Language
    gjf         # Gaussian input data file
    gromacs     # GROMACS file format
    hyperchem   # HyperChem file format
    jme         # Java Molecule Editor format
    maestro     # Schroedinger MacroModel structure file format
    mol         # Symyx molecule file
    mol2        # Tripos Sybyl MOL2 format
    pdb         # Protein Data Bank
    sdf         # 2D formatted structure data files
    sdf3000     # Symyx Structure Data Format 3000
    sln         # SYBYL Line Notation
    xyz         # xyz file format


Properties
----------
