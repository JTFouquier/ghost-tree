#!/usr/bin/env python
# ----------------------------------------------------------------------------
# Copyright (c) 2015--, ghost-tree development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the LICENSE file, distributed with this software.
# ----------------------------------------------------------------------------
# setup.py copied and modified from https://github.com/biocore/scikit-bio
from setuptools import find_packages, setup
from glob import glob

__version__ = "0.0.1-dev"

classes = """
    Development Status :: 1 - Planning
    License :: OSI Approved :: BSD License
    Topic :: Scientific/Engineering
    Topic :: Scientific/Engineering :: Bio-Informatics
    Programming Language :: Python :: 3.5
    Operating System :: Unix
    Operating System :: POSIX
    Operating System :: MacOS :: MacOS X
"""
classifiers = [s.strip() for s in classes.split('\n') if s]

description = ('Tool for creating hybrid-gene phylogenetic trees')

with open('README.rst') as f:
    long_description = f.read()

setup(name='ghost-tree',
      version=__version__,
      license='BSD',
      description=description,
      long_description=long_description,
      author="ghost-tree development team",
      author_email="jennietf@gmail.com",
      maintainer="ghost-tree development team",
      maintainer_email="jennietf@gmail.com",
      test_suite='nose.collector',
      packages=find_packages(),
      scripts=glob("scripts/*"),
      install_requires=['scipy>=0.16.1', 'scikit-bio >=0.5.1', 'click >= 4.0'],
      extras_require={'test': ["nose >= 0.10.1", "pep8", "flake8"]},
      classifiers=classifiers
      )
