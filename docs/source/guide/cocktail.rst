.. _cocktail:

cocktail-shaker API Documentation
=================================

This page introduces the functionality of the cocktail object and provides a deeper look into what it can accomplish in the future.


The cocktail Class
------------------

The cocktail object is the heart of the package and allows you to parse in a peptide SMILES string and expand your drug
library with different functional groups as well as enumerate representations of the molecule in 1D and 2D.

You instantiate a ``Cocktail`` object by parsing in a list of SMILES.
Cocktail-Shaker will already handle the SMILES to RDKit Mol Object for you without having to subject to create them
yourself.
If the SMILES fail to load and is *not* supported, then ``MoleculeError`` will be raised instead.

>>> from cocktail_shaker import Cocktail
>>> from cocktail_shaker import PeptideMolecule
>>> peptide_molecule = PeptideMolecule(2)
>>> cocktail = Cocktail(peptide_molecule, ligand_library=['Br', 'I'])
>>> print (cocktail)
>>> CocktailObject

cocktail-shaker, under the hood, uses two validation methods to determine whether a molecule is a legitimate molecule.
One method uses the MolVS validator and the other using RDKit's rendering capabilities from SMILES to RDKit Mol Object, we can determine
if a generated molecule is legitimate. By using both validation methods, cocktail-shaker can ensure validation of the SMILES.

.. attribute:: peptide_molecule

  The peptide backbone string generated from PeptideMolecule

.. attribute:: ligand_library

  The list of the ligands you would like installed on the peptide. It can be of any order.

.. attribute:: enable_isomers

  Include stereochemistry and stereoisomers in the the results

.. attribute:: include_amino_acids

  Include the natural amino acids except for Proline (TBD), list of smiles found in :ref:`amino acids <aminoacids>`


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
>>> from cocktail_shaker import PeptideMolecule
>>> peptide_molecule = PeptideMolecule(2)
>>> cocktail = Cocktail(peptide_molecule, ligand_library=['Br', 'I'])
>>> print (cocktail.shake())
>>> ['NC(Br)C(=O)NC(I)C(=O)NCC(=O)O', 'NC(I)C(=O)NC(Br)C(=O)NCC(=O)O']

As mentioned before, validation occurs internally. This means that the molecules being generated are validated
before appended the list.

There are some minor restrictions in the first release of this package as described below.

1. Performance is notable with higher amounts of peptides especially rendering in 3D takes ~ 2 days.
2. The limit of functional groups cocktail shaker is valid i.e Azides and Phosphorous groups
3. Circular peptides are also not supported yet.

Because cocktail-shaker uses SMART pattern recognition to detect functional groups, it is limited by how fast we can
support new groups and thoroughly test them. The library is looking to expand into more classes and a variety of
groups.

Please refer to :ref:`functional groups <functionalgroups>` to see what cocktail-shaker supports.
Please refer to :ref:`amino acids <aminoacids>` to see what cocktail-shaker supports.


The "Enumerate" Module
----------------------

The enumerate module takes your RDKit molecule objects and generates random representations of the compounds in either
1D, 2D, and coming soon 3D.

Enumeration does not take into account tautomers, salts, or other configurations just yet but it's on the Roadmap.

You instantiate a ``Cocktail`` object by parsing in a list of smiles and then "enumerate" the compounds.
If the enumerate fails to work then ``MoleculeError`` will be raised instead. In this case, please contact the lead developer.

>>> from cocktail_shaker import Cocktail
>>> from cocktail_shaker import PeptideMolecule   >>> peptide_backbone = PeptideMolecule(2)
>>> cocktail = Cocktail(peptide_backbone,ligand_library = ['Br', 'I'])
>>> combinations = cocktail.shake()
>>> print (combinations)
>>> ['NC(Br)C(=O)NC(I)C(=O)NCC(=O)O', 'NC(I)C(=O)NC(Br)C(=O)NCC(=O)O']
>>> enumerations = cocktail.enumerate()
>>> print (enumerations)
>>> ['IC(C(NCC(=O)O)=O)NC(=O)C(Br)N', 'N(CC(O)=O)C(C(I)NC(=O)C(N)Br)=O', 'NC(C(NC(I)C(NCC(=O)O)=O)=O)Br',
>>>  'OC(=O)CNC(C(NC(C(N)Br)=O)I)=O', 'IC(NC(C(N)Br)=O)C(NCC(=O)O)=O', 'N(C(=O)C(N)Br)C(C(NCC(=O)O)=O)I',
>>>  'O=C(C(N)Br)NC(I)C(NCC(=O)O)=O', 'C(C(O)=O)NC(C(NC(C(N)Br)=O)I)=O', 'OC(=O)CNC(=O)C(NC(=O)C(Br)N)I',
>>>  'N(C(=O)C(I)NC(=O)C(Br)N)CC(O)=O', 'O=C(C(Br)NC(C(I)N)=O)NCC(=O)O', 'O=C(NCC(=O)O)C(Br)NC(C(N)I)=O',
>>>  'N(CC(=O)O)C(C(Br)NC(C(N)I)=O)=O', 'N(C(C(=O)NCC(=O)O)Br)C(C(I)N)=O', 'O=C(CNC(C(Br)NC(=O)C(I)N)=O)O',
>>>  'OC(CNC(C(Br)NC(C(N)I)=O)=O)=O', 'OC(CNC(C(Br)NC(=O)C(I)N)=O)=O', 'C(NC(C(I)N)=O)(C(=O)NCC(=O)O)Br',
>>>  'BrC(C(NCC(O)=O)=O)NC(C(N)I)=O', 'O=C(C(Br)NC(C(N)I)=O)NCC(O)=O']


Alternatively, if you would like you can pass in the enumeration_complexity argument to change how many enumerations
are generated.

>>> from cocktail_shaker import Cocktail
>>> from cocktail_shaker import PeptideMolecule   >>> peptide_backbone = PeptideMolecule(2)
>>> cocktail = Cocktail(peptide_backbone,ligand_library = ['Br', 'I'])
>>> combinations = cocktail.shake()
>>> print (combinations)
>>> ['NC(Br)C(=O)NC(I)C(=O)NCC(=O)O', 'NC(I)C(=O)NC(Br)C(=O)NCC(=O)O']
>>> enumerations = cocktail.enumerate(enumeration_complexity='low')
>>> print (len(enumerations))
>>> 20
>>> enumerations = cocktail.enumerate(enumeration_complexity='med')
>>> print (len(enumerations))
>>> 186
>>> enumerations = cocktail.enumerate(enumeration_complexity='high')
>>> print (len(enumerations))
>>> 1789

Cocktail Shaker also allows you to pass in the dimensionality of the enumeration.

>>> from cocktail_shaker import Cocktail
>>> from cocktail_shaker import PeptideMolecule   >>> peptide_backbone = PeptideMolecule(2)
>>> cocktail = Cocktail(peptide_backbone,ligand_library = ['Br', 'I'])
>>> combinations = cocktail.shake()
>>> enumerations = cocktail.enumerate(dimensionality='2D')

Coming soon is mol2 3D Enumeration and on the roadmap as the big feature item for 2.0.

The enumeration works by following the algorithm of generating random SMILES generated by RDKit. This allows
for different representation in 1D format. Coincidentally, this algorithm works for 2D. 3D files are a little more
complex in terms of enumeration, but is on track for version 2.0 release.

The enumeration complexity refers to how many times cocktail-shaker will try to generate a unique random SMILES
representation. This goes with order of magnitude of 10.

.. attribute:: enumeration_complexity (optional)

     How many representations would you like to generate. Defaults to 'Low'
     'low'    = 10 Representations
     'medium' = 100 Representations
     'high    = 1000 Representations

.. attribute:: dimensionality (optional)

     What dimensionality you would like i.e '1D', '2D', '3D (Not Supported)'. Defaults to '1D'.




