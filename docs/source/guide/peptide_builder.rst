.. _peptidebuilder:

Peptide Builder API
===================

This page introduces the functionality of the peptide molecule object and provides a deeper look into what it can accomplish in the future.


The PeptideBuilder Class
-------------------------

The peptide molecule class is the preliminary peptide builder work before you pass it into Cocktail Shaker. It will build
the peptides smiles for you (woo!)

You instantiate a ``PeptideBuilder`` object by parsing in the length of peptide marked by how many amino acid sides you
would like. Peptide Molecule will handle any SMARTS and validation under the hood so you will not have too. A simple example
is detailed below:

>>> from cocktail_shaker import PeptideBuilder
>>> peptide_molecule = PeptideBuilder(2)
>>> print (peptide_molecule)
>>> NCC(NC([*:1])C(NC([*:2])C(O)=O)=O)=O

Here we have two slots ready for our combinations in the cocktail. Alternatively, you can include the ``include_proline``
function on the n-terminus that will provide a SMILES string like so.

>>> from cocktail_shaker import PeptideBuilder
>>> peptide_molecule = PeptideBuilder(2, include_proline=True)
>>> print (peptide_molecule)
>>> N2CCCC2C(NC([*:1])C(NC([*:2])C(O)=O)=O)=O

Generation of a circular peptide

>>> from cocktail_shaker import PeptideBuilder
>>> peptide_molecule = PeptideBuilder(3, circular=True)
>>> print (peptide_molecule)
>>> O=C1C([*:1])NC(=O)C([*:2])NC(=O)C([*:3])N1

.. attribute:: length_of_peptide

  Length of the peptide the user would like generated

.. attribute:: include_proline (optional)

  Whether the user would like to include proline on the N-Terminus. Defaults to False.

.. attribute:: circular (optional)

  If you would like the peptide to be circular or not.