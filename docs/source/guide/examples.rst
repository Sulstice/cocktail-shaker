.. _examples:

Examples
========

This page instructs you on using some of the features of cocktail shaker and what you can expect.

Stereoisomers
-------------

Using the stereoisomer function with Cocktail Shaker

>>> from cocktail_shaker import Cocktail
>>> from cocktail_shaker import PeptideMolecule
>>> peptide_backbone = PeptideMolecule(1)
>>> cocktail = Cocktail(
>>>     peptide_backbone,
>>>     ligand_library = ['Br'],
>>>     enable_isomers = True
>>> )
>>> combinations = cocktail.shake()
>>> print (combinations)
>>> ['N[C@H](Br)C(=O)NCC(=O)O', 'N[C@@H](Br)C(=O)NCC(=O)O']


Amino Acids
-----------

Using the include amino acid function

>>> from cocktail_shaker import Cocktail
>>> from cocktail_shaker import PeptideMolecule
>>> peptide_backbone = PeptideMolecule(1)
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

