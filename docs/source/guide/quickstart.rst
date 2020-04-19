.. _gettingstarted:

Getting Started
===============

This page instructs you on how to get started with cocktail-shaker. To begin, make sure you have performed
:ref:`installing cocktail_shaker <install>`.

Basic Usage
-----------

The simplest way to use cocktail-shaker is to create a Peptide Builder object to generate your "base" peptide string,
then to create cocktail object passing the peptide backbone and with the ``shake`` function you can generate your combinations:

    >>> from cocktail_shaker import Cocktail
    >>> from cocktail_shaker import PeptideBuilder
    >>> peptide_backbone = PeptideBuilder(2)
    >>> cocktail = Cocktail(
    >>>                     peptide_backbone,
    >>>                     ligand_library = ['Br', 'I']
    >>>                     )
    >>> combinations = cocktail.shake()
    >>> print (combinations)
    >>> ['NC(Br)C(=O)NC(I)C(=O)NCC(=O)O', 'NC(I)C(=O)NC(Br)C(=O)NCC(=O)O']

Write the new compounds into an SDF file:

    >>> FileWriter('new_compounds', new_compounds, 'sdf')

In this example, we have taken one SMILES string and expanded the compounds into a variety of variations into one SDF file.

Cocktail Shaker uses the PeptideBuilder class to generate base peptide backbone strings that can be then passed as into the
cocktail shaker object. Alternatively, a user can generate their own peptide string as long as it conforms with the cocktail shaker
requirements. See the peptide builder documentation for more information.

More Examples Peptide Builder
-----------------------------

Generation of a circular peptide

>>> from cocktail_shaker import PeptideBuilder
>>> peptide_molecule = PeptideBuilder(3, circular=True)
>>> print (peptide_molecule)
>>> O=C1C([*:1])NC(=O)C([*:2])NC(=O)C([*:3])N1

Using the stereoisomer function with Cocktail Shaker

>>> from cocktail_shaker import Cocktail
>>> from cocktail_shaker import PeptideBuilder
>>> peptide_backbone = PeptideBuilder(1)
>>> cocktail = Cocktail(
>>>     peptide_backbone,
>>>     ligand_library = ['Br'],
>>>     enable_isomers = True
>>> )
>>> combinations = cocktail.shake()
>>> print (combinations)
>>> ['N[C@H](Br)C(=O)NCC(=O)O', 'N[C@@H](Br)C(=O)NCC(=O)O']


More Examples of Cocktail Shaker
--------------------------------

Using the include amino acid function

>>> from cocktail_shaker import Cocktail
>>> from cocktail_shaker import PeptideBuilder
>>> peptide_backbone = PeptideBuilder(1)
>>> cocktail = Cocktail(
>>>     peptide_backbone,
>>>     include_amino_acids = True
>>> )
>>> combinations = cocktail.shake()
>>> print (combinations)
>>> ['NCCCCC(N)C(=O)NCC(=O)O', 'CC(O)C(N)C(=O)NCC(=O)O', 'NC(Cc1ccc(O)cc1)C(=O)NCC(=O)O', 'NC(=O)CCC(N)C(=O)NCC(=O)O',
>>>  'CCC(C)C(N)C(=O)NCC(=O)O', 'NC(CS)C(=O)NCC(=O)O', 'NC(CC(=O)O)C(=O)NCC(=O)O', 'N=C(N)NCCCC(N)C(=O)NCC(=O)O',
>>>  'NC(Cc1c[nH]cn1)C(=O)NCC(=O)O', 'NC(Cc1ccccc1)C(=O)NCC(=O)O', 'CC(N)C(=O)NCC(=O)O', 'NC(CCC(=O)O)C(=O)NCC(=O)O',
>>>  '[H]C(N)C(=O)NCC(=O)O', 'CSCCC(N)C(=O)NCC(=O)O', 'NC(CO)C(=O)NCC(=O)O', 'CC(C)C(N)C(=O)NCC(=O)O',
>>>  'NC(CCc1c[nH]c2ccccc12)C(=O)NCC(=O)O', 'CC(C)CC(N)C(=O)NCC(=O)O']

Using the Cocktail Shaker to generate a library of halogens & single atoms with then converting into a Pandas DataFrame
Note to use this example you will need to install an extra dependency of 'tables' for handling h4 data.

>>> from peptide_builder import PeptideBuilder
>>> from functional_group_enumerator import Cocktail
>>> import pandas as pd
>>>
>>>
>>> peptide_molecule = PeptideBuilder(length_of_peptide=7)
>>> cocktail = Cocktail(peptide_backbone=peptide_molecule,
>>>                     ligand_library=["Br", "Cl", "I", "F", "O", "N", "C"],
>>>                     include_amino_acids=False,
>>>                     enable_isomers=False)
>>> molecules = cocktail.shake()
>>>
>>> dataframe = pd.DataFrame(molecules, columns=["Smiles"])
>>> dataframe.to_hdf('data.h5', key='s', mode='w')