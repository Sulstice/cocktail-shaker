---
title: 'Cocktail Shaker: An open source drug expansion and enumeration library for peptides'
tags:
  - Python
  - Cheminformatics
  - RDKit
authors:
  - name: Suliman Sharif
    orcid: 0000-0002-1342-9258

date: 23 November 2019
bibliography: paper.bib
---

# Statement of Need

Without expensive software, the rapid creation and design of peptide ligand libraries has been a
challenge for many drug discovery scientists. Currently, proteins and peptide-based therapeutics consist of 10% of the 
pharmaceutical market and will make up a larger proportion of the market in the future [1, 2]. With increasing interest
eludes to increasing high throughput screening of peptides computationally. There exists two platforms: Molecular Operating 
Environment (MOE) [3] and Rapid Peptides Generator (RPG) [4] but both have drawbacks. MOE works efficiently in creating peptide molecule 3D chemical 
files in one particular format (mol2) but also at the high price cost for a licence. RPG, although free of charge, does not account
for non-natural amino acids and production of multiple chemical files. In this study, I present the first open source python package,
 ```Cocktail Shaker```, developed for exploring, expanding, and synthesizing chemical peptide data.

# Methodology and Implementation

```Cocktail Shaker``` operates within the RDKit platform [5] and is designed for the chemically-oriented
computational research community. RDKit utilizes C-based functions for speed and rapid creation of molecule objects.
The toolkit offers a variety of utilities that includes: parsing and producing ready-to-use scientific files
designed for any chemical software, employing click chemistry methods for ease of exchange compounds, chemical data writing, 
and chemical representation enumeration employed for machine learning.

```Cocktail Shaker``` consists of four major class objects available to the user: PeptideMolecule, CocktailShaker, and FileWriter.
Using string manipulation Peptide Molecule can build SMILES strings with allocated slots defined by the user. The user can
then enter the produced SMILES into the CocktailShaker object with a library of ligands represented by smiles and optional arguments
of whether to include generation of stereoisomers and/or natural amino acids.  
Cocktail Shaker will generate all combinations of the library and allocate them to a slot within the peptide. This process of 
string manipulation is presented in ```Figure 1```.

<figure>
  <img src="images/figure_1.png" alt="Figure 1 String Manipulation Diagram"/>
  <figcaption><i>Figure 1: Full string manipulation diagram of Cocktail Shaker.</i></figcaption>
</figure>

```Cocktail Shaker``` also allows for File Writing of the molecules into a wide array of chemical formats (found in the documentation).
```Cocktail Shaker ``` uses RDKit to convert from 1D to 2D and the CIR Resolver built from webchem to convert 1D SMILES to 3D. At the request 
of the user store the data return into one large data file or seperate using the keyword fragmentation. This additional API
allows the user the flexibly to write a variety of files to implement in their respective chemical software.  

# Conclusion

Using  ```Cocktail Shaker```, individual research groups and companies can quickly construct private compound collections and progressively improve public
libraries with increased variations of chemical compound data.

```Cocktail Shaker``` with it's first version release provides a basis drug expansion and enumeration. For future
releases ```Cocktail Shaker``` will be expanding into specifying shapes of compounds and, recently partnered with MolPort,
vendor information on any compound generated. It was presented at the RDKit UGM conference at the University of Hamburg
to the cheminformatics community with positive feedback with it's second version 1.0.1. With incorporated feedback it will now be
released with version 1.1.0. With more contributions ```Cocktail Shaker``` will be an exciting time for drug library creation 
and drug discovery for scientists and engineers alike. 


# Citations

Citations to entries in paper.bib should be in

For a quick reference, the following citation commands can be used:
- `@author:2001`  ->  "Author et al. (2001)"
- `[@author:2001]` -> "(Author et al., 2001)"
- `[@author1:2001; @author2:2001]` -> "(Author1 et al., 2001; Author2 et al., 2002)"

# Acknowledgements

We acknowledge contributions from Ryland Forsythe as an academic advisor, Marvin Corro for quality assurance, Rose Gierth
for technical documentation, and Elena Chow for her work on the graphics. 

# References

1. Craik DJ, Fairlie DP, Liras S, Price D. The future of peptide-based drugs. Chem Biol Drug Des. 2013;81(1):136–47. Epub 2012/12/21. pmid:23253135.
2. Bruno BJ, Miller GD, Lim CS. Basics and recent advances in peptide and protein drug delivery. Ther Deliv. 2013;4(11):1443–67. Epub 2013/11/16. pmid:24228993; PubMed Central PMCID: PMCPMC3956587.
3. Reynolds CH, Merz KM, Ringe D, eds. (2010). Drug Design: Structure- and Ligand-Based Approaches (1 ed.). Cambridge, UK: Cambridge University Press. ISBN 978-0521887236.
4. Nicolas Maillet, Rapid Peptides Generator: fast and efficient in silico protein digestion, NAR Genomics and Bioinformatics, Volume 2, Issue 1, March 2020, lqz004, https://doi.org/10.1093/nargab/lqz004
5. Landrum, Greg, RDKit: Open-Source Cheminformatics Software, release 2019 
   
   
