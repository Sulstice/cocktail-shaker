.. _fileformats:

File Formats
============

Input can be parsed in a variety of file formats that are autodetected::

    >>> FileParser('compound.sdf')

    ['RDKit MolObject']

The full list of file formats that cocktail-shaker supports parsing::

    sdf       # 2D Structure file format
    mol2      # (Coming soon)


Output can be returned in a variety of file formats that are specified your enumerated compounds and your
desired file extension::

    >>> FileWriter('test'', compounds, '.mol2')

    Writes all the compounds to 1 mol2 file.


The full list of file formats::

    alc         # Alchemy format
    cdxml       # CambridgeSoft ChemDraw XML format
    cerius      # MSI Cerius II format
    charmm      # Chemistry at Harvard Macromolecular Mechanics file format
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
