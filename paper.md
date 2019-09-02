---
title: 'Cocktail Shaker: An open source drug expansion and enumeration library using Python and RDKit.'
tags:
  - Python
  - Cheminformatics
  - RDKit
authors:
  - name: Suliman Sharif
    orcid: 0000-0002-1342-9258

date: 01 September 2019
bibliography: paper.bib
---

# Summary

Without expensive software, the rapid creation and design of ligand libraries has been a
challenge for many drug discovery scientists. In this study, I present the first open source
python package, ```Cocktail Shaker```, developed for exploring, expanding, and synthesizing chemical
datasets.

```Cocktail Shaker``` operates within the RDKit platform and is designed for the chemically-oriented
computational research community. RDKit utilizes C-based functions for speed and rapid creation of molecule objects.
The toolkit offers a variety of utilities that includes: parsing and producing ready-to-use scientific files
designed for any chemical software, detecting common functional groups (FGs) based on SMART pattern matching, 
employing click chemistry methods for ease of exchange compounds, chemical data writing, 
and chemical representation enumeration designed for machine learning.

The toolkit is loaded with pre-set functional groups in its datasource, which is free to edit and
expand for the open source community. Using  ```Cocktail Shaker```, individual research groups and
companies can quickly construct private compound collections and progressively improve public
libraries with increased variations of chemical compound data.


```Cocktail Shaker``` with it's first version release provides a basis drug expansion and enumeration. For future
releases ```Cocktail Shaker``` will be expanding into specifying shapes of compounds and, recently partnered with MolPort,
vendor information on any compound generated. It will be presented at the RDKit UGM conference at the University of Hamburg
 to the cheminformatics community and hopefully be used in pharma, and academia alike. With more contributions ```Cocktail Shaker```
 will be an exciting time for drug library creation and drug discovery for scientists and engineers alike. 


# Citations

Citations to entries in paper.bib should be in

For a quick reference, the following citation commands can be used:
- `@author:2001`  ->  "Author et al. (2001)"
- `[@author:2001]` -> "(Author et al., 2001)"
- `[@author1:2001; @author2:2001]` -> "(Author1 et al., 2001; Author2 et al., 2002)"

# Figures

Figures can be included like this: ![Example figure.](figure.png)

# Acknowledgements

We acknowledge contributions from Ryland Forsythe, and Marvin Corro for QA Testing the application. 
And support from Elena Chow for her contribution for the graphic design. 

# References
