#from setuptools import setup, find_packages
from distutils.core import setup
import os, sys

__version__ = '0.1'

setup(name = 'mNuCosmoMap',
      version = __version__,
      description = 'TBD',
      author='ChangHoon Hahn',
      author_email='changhoonhahn@lbl.gov',
      url='',
      platforms=['*nix'],
      license='GPL',
      requires = ['numpy', 'matplotlib', 'scipy'],
      provides = ['mNuCosmoMap'],
      packages = ['mnucosmomap'],
      scripts=['mnucosmomap/catalogs.py', 'mnucosmomap/util.py']
)
