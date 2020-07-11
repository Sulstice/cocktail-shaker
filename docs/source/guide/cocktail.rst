.. _cocktail:

Cocktail API
============

This page introduces the functionality of the cocktail object and provides a deeper look into what it can accomplish in the future.


The GlobalChem Class
--------------------

.. attribute:: amino_acid_side_chains

  The peptide backbone string generated from PeptideBuilder

.. attribute:: ligand_library

  The list of the ligands you would like installed on the peptide. It can be of any order.

.. attribute:: enable_isomers (optional)

  Include stereochemistry and stereoisomers in the the results. Defaults to False.

.. attribute:: include_amino_acids (optional)

  Include the natural amino acids except for Proline (TBD), list of smiles found in :ref:`amino acids <aminoacids>`.
  Defaults to False.
