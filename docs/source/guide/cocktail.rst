.. _cocktail:

cocktail-shaker API Documentation
=================================

This page introduces the functionality of the cocktail object and provides a deeper look into what it can accomplish in the future.

The cocktail Class
------------------

The cocktail object is the heart of the package and allows you to parse in a SMILES string and expand your drug
library with different functional groups as well as enumerate representations of the molecule in 1D and 2D.

    .. attribute:: molecules

      The list of SMILES you would like to be passed into the cocktail.

    You instantiate a ``Cocktail`` object by parsing in a list of SMILES.
    Cocktail-Shaker will already handle the SMILES to RDKit Mol Object for you without having to subject to create them
    yourself.
    If the SMILES fail to load and is *not* supported, then ``MoleculeError`` will be raised instead.

    >>> from cocktail_shaker import Cocktail
    >>> cocktail = Cocktail(['c1cc(CCCO)ccc1'])
    >>> cocktail = Cocktail(['c1cc(CCCO)ccc1', 'c1cc(CCCBr)ccc1'])
    >>> print (cocktail)
    >>> CocktailObject

    cocktail-shaker, under the hood, uses two validation methods to determine whether a molecule is a legitimate molecule.
    One method uses the MolVS validator and the other method is internal. 

    Using RDKit's rendering capabilities from SMILES to RDKit Mol Object, we can determine
    if a generated molecule is legitimate. By using both validation methods, cocktail-shaker can ensure validation of the SMILES.


The "Shake" Module
------------------

    The shake module detects functional groups present on a molecule, breaks their bond, and then adds a functional
    group (not itself) from the datasource library as replacement.

    Essentially, the shake compound is utilizing a chemical method known as "click" chemistry in the snapping and formation
    of new bonds linking them together.

    You instantiate a ``Cocktail`` object by parsing in a list of SMILES and then "shake" the compounds.
    cocktail-shaker shows the functional groups that can be detected and then swaps them accordingly.
    If the shake fails to work, then ``MoleculeError`` will be raised instead. In this case, contact the lead developer.

    >>> from cocktail_shaker import Cocktail
    >>> cocktail = Cocktail(['c1cc(CCCO)ccc1'])
    >>> new_compounds = cocktail.shake()
    >>> print (cocktail)
    [RDKit Mol Object, RDKit Mol Object, RDKit Mol Object]

    As mentioned before, validation occurs internally. This means that the molecules being generated are validated
    before appended the list.

    There are some minor restrictions in the first release of this package as described below.

    1. You cannot select the bond to break. This restriction is on the Roadmap for version 2.0.
    2. The limit of functional groups cocktail shaker can detect.

    Because cocktail-shaker uses SMART pattern recognition to detect functional groups, it is limited by how fast we can
    support new groups and thoroughly test them. The library is looking to expand into more classes and a variety of
    groups.

    cocktail-shaker allows you pass in a functional_groups parameter where you can selectively pick a functional
    group. Please refer to :ref:`functional groups <functionalgroups>` to see what cocktail-shaker supports.

    >>> from cocktail_shaker import Cocktail
    >>> cocktail = Cocktail(['c1cc(CCCO)ccc1'])
    >>> new_compounds = cocktail.shake(functional_groups=['Azides'])
    >>> print (cocktail)
    [RDKit Mol Object (Azide)]

   .. attribute:: functional_groups

      The list of functional groups that you would like to exchange specifically with.

The "Enumerate" Module
----------------------

    The enumerate module takes your RDKit molecule objects and generates random representations of the compounds in either
    1D, 2D, and coming soon 3D.

    Enumeration does not take into account tautomers, salts, or other configurations just yet but it's on the Roadmap.

    You instantiate a ``Cocktail`` object by parsing in a list of smiles and then "enumerate" the compounds.
    If the enumerate fails to work then ``MoleculeError`` will be raised instead. In this case, please contact the lead developer.

    >>> from cocktail_shaker import Cocktail
    >>> cocktail = Cocktail(['c1cc(CCCO)ccc1'])
    >>> new_compounds = cocktail.enumerate(enumeration_complexity='low', dimensionality='2D'))
    >>> print (cocktail)
    [RDKit Mol Object (2D Representation), RDKit Mol Object (2D Representation), RDKit Mol Object (2D Representation)]

    Alternatively, if you have just shook the compounds, cocktail-shaker is smart enough to grab the previously generated
    new compounds and apply the shake.

    >>> from cocktail_shaker import Cocktail
    >>> cocktail = Cocktail(['c1cc(CCCO)ccc1'])
    >>> cocktail.shake()
    >>> new_compounds = cocktail.enumerate(enumeration_complexity='low', dimensionality='2D'))
    >>> print (cocktail)
    [RDKit Mol Object (2D Representation), RDKit Mol Object (2D Representation), RDKit Mol Object (2D Representation)]

    The enumeration works by following the algorithm of generating random SMILES generated by RDKit. This allows
    for different representation in 1D format. Coincidentally, this algorithm works for 2D. 3D files are a little more
    complex in terms of enumeration, but is on track for version 2.0 release.

    The enumeration complexity refers to how many times cocktail-shaker will try to generate a unique random SMILES
    representation. This goes with order of magnitude of 10.

   .. attribute:: enumeration_complexity

        How many representations would you like to generate.
        'low'    = 10 Representations
        'medium' = 100 Representations
        'high    = 1000 Representations

   .. attribute:: dimensionality

        What dimensionality you would like i.e '1D', '2D', '3D (Not Supported)'




