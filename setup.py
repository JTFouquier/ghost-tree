#!/usr/bin/env python

# setup.py copied and modified from https://github.com/biocore/scikit-bio
from setuptools import find_packages, setup
from glob import glob

__version__ = "0.0.1-dev"

classes = """
    Development Status :: 1 - Planning
    License :: OSI Approved :: BSD License
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
      install_requires=['scikit-bio >= 0.2.3, < 0.3.0', 'click >= 4.0'],
      extras_require={'test': ["nose >= 0.10.1", "pep8", "flake8"]},
      classifiers=classifiers
      )
