.. _filehandling:

File API Documentation
======================

This page provides an overview of how cocktail-shaker operates reading and writing files. To see
the files that cocktail-shaker supports, refer to :ref:`file formats <fileformats>`.

The FileParser Module
---------------------

You instantiate a ``FileParser``
by providing the exact path to a file with the extension included.
cocktail-shaker is smart enough to detect the file extension and allocate its specific parsing.
If the file being parsed is *not* supported, then ``FileNotSupportedError`` will be raised instead.

>>> from cocktail_shaker import FileParser
>>> molecules = FileParser('compounds.sdf')
>>> print (molecules)
>>> ['c1cc(CCCO)ccc1']

FileParser will then return a list of SMILES.

.. attribute:: file_path

   The path to a file


The FileWriter Module
---------------------

You instantiate a ``FileWriter``
by providing the path to the file, compounds to be written, and the extension you would like the files in.
If the file being written is *not* supported, then ``FileNotSupportedError`` will be raised instead.

>>> from cocktail_shaker import Cocktail, FileWriter
>>> from cocktail_shaker import PeptideBuilder
>>> peptide_backbone = PeptideBuilder(2)
>>> cocktail = Cocktail(peptide_backbone,ligand_library = ['Br', 'I'])
>>> combinations = cocktail.shake()
>>> FileWriter('new_compounds', combinations, 'sdf')
Generates an SDF File....

However, if you would like to generate the files into separate files, then you can pass in the fragmentation parameter.

>>> from cocktail_shaker import Cocktail, FileWriter
>>> from cocktail_shaker import PeptideBuilder
>>> peptide_backbone = PeptideBuilder(2)
>>> cocktail = Cocktail(peptide_backbone,ligand_library = ['Br', 'I'])
>>> combinations = cocktail.shake()
>>> FileWriter('new_compounds', combinations, 'sdf', fragmentation=2)
Generates 2 SDF Files....

.. attribute:: name

   The path to a file

.. attribute:: molecules

   List of RDKit molecules you would like to write.

.. attribute:: option

   The extension of the file you would like to write

.. attribute:: fragmentation (optional)

   How many files you would like to produce.
