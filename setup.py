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

# # files
# # -----
# with open('README.md') as fi:
#     README = fi.read()


# exec
# ----
setup(
    name="cocktail_shaker",
    version="1.0.0-beta",
    packages = ['cocktail_shaker'],
    license='MIT',
    author="Suliman Sharif",
    author_email="sharifsuliman1@gmail.com",
    url="https://www.github.com/Sulstice/Cocktail-Shaker",
    install_requires=REQUIREMENTS,
    zip_safe=False,
    keywords=['cocktail', 'chemistry', 'ligand-design', 'shaker'],
    classifiers=[
        'Development Status :: 4 - Beta',
        'Natural Language :: English',
        'License :: OSI Approved :: MIT License',
        'Intended Audience :: Developers',
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
