Cocktail Shaker: Drug Expansion and Enumeration for Peptides!
=============================================================

[![Build](https://travis-ci.org/Sulstice/Cocktail-Shaker.svg?branch=master)](https://travis-ci.org/Sulstice/Cocktail-Shaker)
[![License: MPL 2.0](https://img.shields.io/badge/License-MPL%202.0-brightgreen.svg)](https://opensource.org/licenses/MPL-2.0)
[![Coverage](https://coveralls.io/repos/github/Sulstice/Cocktail-Shaker/badge.svg?branch=master)](https://coveralls.io/github/Sulstice/Cocktail-Shaker?branch=master)
![Python](https://img.shields.io/badge/python-3.6-blue.svg)
[![Gitter](https://badges.gitter.im/Cocktail-Shaker/community.svg)](https://gitter.im/Cocktail-Shaker/community?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)
[![Zenodo](https://zenodo.org/badge/170644606.svg)](https://zenodo.org/badge/latestdoi/170644606)
[![Documentation Status](https://readthedocs.org/projects/cocktail-shaker/badge/?version=latest)](https://cocktail-shaker.readthedocs.io/en/latest/?badge=latest)
[![status](https://joss.theoj.org/papers/c2e1d3c408a5729d832b34ac680d6305/status.svg)](https://joss.theoj.org/papers/c2e1d3c408a5729d832b34ac680d6305)

<p align="center">
  <img width="200" height="300" src="images/logoshaker.png">
</p>

cocktail-shaker is a **high-performance drug enumeration and expansion
library**. cocktail-shaker leverages the computational power of
**RDKit** to create and enumerate large volumes of drug compounds. 

Announcements
=============

-   **Release!** Version 1.0.0-beta, August 26, 2019
-   **RDkit UGM 2019**: talk at September 25th at the University of
    Hamburg, Germany.

Using Cocktail Shaker
=====================

cocktail-shaker is a young library under heavy development at this time.
It targets two categories of users:

1.  **Users familiar with RDKit**, or those willing to learn RDKit, who want to
    create fast sets of data for high throughput screening or machine
    learning.
2.  **Open-Science Scientists without any knowledge of RDKit** who are
    seeking a a high-level wrapper to create chemical files for their
    software.

If you're in the first category, then you can already start using RDKit.
cocktail-shaker offers a Pythonic, easy-to-use library and you can start
channeling molecules in the expansion library. Instead of validating the
sanity of the data, cocktail-shaker takes care of that for you. With
each molecule being generated it will head into a 1D and/or 2D
validation check (3D not supported yet).

If you're in the second category, we're starting to build experimental
high-level python code to take care a lot of the underpinnings of RDKit.

Installation 
==================

Cocktail-shaker runs on Python 3.5+ and within a conda environment due to the RDKit dependency. A conda installation is coming soon (as soon as I figure it out). 

For the time being, here is the following:

- [Anaconda](https://docs.anaconda.com/anaconda/install/)
- [RDKit](https://www.rdkit.org/docs/Install.html) Version: 2019.09.1

Cocktail shaker is distributed through PyPi and can be installed via:

`your_env/bin/python -m pip install cocktail-shaker`



Development Installation
========================

As cocktail-shaker is under heavy development at this time, we highly
recommend developers to use the development version on Github (master
branch). You need to clone the repository and install cocktail-shaker with

`python setup.py install`.

As a one-liner, assuming git is installed:

    git clone https://github.com/Sulstice/Cocktail-Shaker.git

This will automatically install the latest version of cocktail-shaker.

Structure of cocktail-shaker
============================

Currently, the main subpackages are:

-   **cocktail_shaker**: Contains a lot of the high level functionality; request
    handling, file parsing/writing, enumeration, and expansion.
-   **docs**: An access point for the readthedocs implementation.
-   **cocktail_shaker/datasources**: This is where the system stores its data on
    predfined functional groups and/or shapes (coming soon).
-   **tests**: Tests that are for the file handling, requests, and
    testing molecule pattern recognition.

The API of all public interfaces are subject to change in the future,
although **datasources** are *relatively* stable at this point.

Genesis
=======

cocktail-shaker began when one developer/scientist wanted an open source
drug library.

- Lead Developer [Suliman sharif](http://sulstice.github.io/)
- Artwork [Elena Chow](http://www.chowelena.com/)
- Technical Documentation [Rose Gierth](https://www.linkedin.com/in/rose-gierth-69a4a083/)
- QA Tester [Marvin Corro](https://www.linkedin.com/in/marvincorro/)

Now cocktail-shaker looks to build on the expertise of these
developers/scientists and the broader open-science community to build an
effective drug library.

* * * * *

External links
==============

-   [Documentation](http://cocktail-shaker.readthedocs.org)

