.. _filehandling:

File API Documentation
======================

This page gives a introduction on how cocktail-shaker operates reading and writing files. Let's get started! To see
what files we support please head to :ref:`file formats <fileformats>`

The FileParser module
---------------------

    You instantiate a ``FileParser``
    by providing exactly the path to a file with the extension included.
    Cocktail Shaker is smart enough to detect the file extension and allocate it's specific parsing.
    If the file being parsed *is not* supported, then ``FileNotSupportedError`` will be raised instead.

    >>> from cocktail_shaker import FileParser
    >>> molecules = FileParser('compounds.sdf')
    >>> print (molecules)
    >>> ['c1cc(CCCO)ccc1']

    FileParser will then return a list of SMILES.

    .. attribute:: file_path

       The path to a file


The FileWriter module
---------------------

    You instantiate a ``FileWriter``
    by providing the path to the file, compounds to be written, and the extension you would like the files in.
    If the file being written *is not* supported, then ``FileNotSupportedError`` will be raised instead.

    >>> from cocktail_shaker import Cocktail, FileWriter
    >>> cocktail = Cocktail(['c1cc(CCCO)ccc1'])
    >>> molecules = cocktail.shake()
    >>> print (molecules)
    [RDKit_Mol_Object, RDKit_Mol_Object, RDKit_Mol_Object...]
    >>> FileWriter('new_compounds', molecules, '.sdf')
    Generates an SDF File....

    If however you would like to generate the files into separate files you can pass in the fragmentation parameter

    >>> from cocktail_shaker import Cocktail, FileWriter
    >>> cocktail = Cocktail(['c1cc(CCCO)ccc1'])
    >>> molecules = cocktail.shake()
    >>> print (molecules)
    [RDKit_Mol_Object, RDKit_Mol_Object, RDKit_Mol_Object...]
    >>> FileWriter('new_compounds', molecules, '.sdf', fragmentation=2)
    Generates 2 SDF Files....

    .. attribute:: name

       The path to a file

    .. attribute:: molecules

       List of RDKit molecules you would like to write.

    .. attribute:: option

       The extension of the file you would like to write

    .. attribute:: fragmentation (optional)

       How many files you would like to produce.
