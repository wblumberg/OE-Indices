#!/usr/bin/env python

from distutils.core import setup

setup(name='oeindices',
      version='1.0',
      description='Python script to compute convection indices from OE-files',
      author='Greg Blumberg',
      author_email='wblumberg@ou.edu',
      py_modules=['oe_indices', 'indices_helper']
)


