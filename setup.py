#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# package setup
#
# ------------------------------------------------

# imports
# -------


# config
# ------
try:
    from setuptools import setup, find_packages
except ImportError:
    from distutils.core import setup, find_packages

# requirements
# ------------
with open('requirements.txt') as f:
    REQUIREMENTS = f.read().strip().split('\n')

TEST_REQUIREMENTS = [
    'pytest',
    'pytest-runner'
]

LIGAND_LIBRARY_LOADER = "Cocktail Shaker"

# files
# -----
with open('README.md') as fi:
    README = fi.read()


# exec
# ----
setup(
    name="Cocktail Shaker",
    version="0.0.1",
    description="Expand your ligand library with this amazing package!",
    long_description=README,
    author="Suliman Sharif",
    author_email="sharifsuliman1@gmail.com",
    url="www.github.com/Sulstice/Cocktail-Shaker",
    install_requires=REQUIREMENTS,
    zip_safe=False,
    keywords=['ligand-loader', 'chemistry', 'ligand-design'],
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers :: Cheminformaticians',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
    ],
    test_suite='tests',
    tests_require=TEST_REQUIREMENTS
)
