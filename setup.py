#!/usr/bin/env python

__version__ = "0.0.1-dev"

from setuptools import find_packages, setup
from setuptools.extension import Extension
from glob import glob

classes = """
    Development Status :: 1 - Planning
    Topic :: Scientific/Engineering
    Topic :: Scientific/Engineering :: Bio-Informatics
    Programming Language :: Python :: 2.7
    Operating System :: Unix
    Operating System :: POSIX
    Operating System :: MacOS :: MacOS X
"""
classifiers = [s.strip() for s in classes.split('\n') if s]

description = ('Tool for creating hybrid 18S + ITS fungal phylogenetic trees')

with open('README.rst') as f:
    long_description = f.read()

setup(name='ghost-tree',
      version=__version__,
      description=description,
      long_description=long_description,
      author="Jennifer Fouquier",
      author_email="jennietf@gmail.com",
      maintainer="Jennifer Fouquier",
      maintainer_email="jennietf@gmail.com",
      test_suite='nose.collector',
      packages=find_packages(),
      scripts=glob("scripts/*.py"),
      install_requires=['scikit-bio>=0.2.2', 'click'],
      extras_require={'test': ["nose >= 0.10.1", "pep8", "flake8"],
                      'doc': ["Sphinx"]},
      classifiers=classifiers,
      )
