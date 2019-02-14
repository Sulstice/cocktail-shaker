#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# package setup
#
# ------------------------------------------------


# config
# ------
try:
    from setuptools import setup, find_packages
except ImportError:
    from distutils.core import setup, find_packages

# requirements
# ------------
with open('requirements.txt', 'r') as reqs:
    REQUIREMENTS = map(lambda x: x.rstrip(), reqs.readlines())

TEST_REQUIREMENTS = [
    'pytest',
    'pytest-runner'
]

DEEP_CURE = "DeepCure"

# files
# -----
with open('README.md') as fi:
    README = fi.read()


# exec
# ----
setup(
    name=DEEP_CURE,
    version=DEEP_CURE.__version__,
    description=DEEP_CURE.__desc__,
    long_description=README,
    author=DEEP_CURE.__author__,
    author_email=DEEP_CURE.__email__,
    url=DEEP_CURE.__url__,
    packages=find_packages(include=[
        'mol2check',
    ]),
    include_package_data=True,
    install_requires=REQUIREMENTS,
    zip_safe=False,
    keywords=['deepcure', 'chemistry', 'ligand-design'],
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
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
