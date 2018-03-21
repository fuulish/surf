#!/usr/bin/env python

import numpy
from setuptools import setup
from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize


ext_modules = [Extension("surf/_surf",
                         sources=["surf/_surf.pyx", "surf/c_surf.c" ],
                         include_dirs=[numpy.get_include()],
                         extra_compile_args = ['-fopenmp'],
                         extra_link_args=['-fopenmp'],
                         libraries=['gsl', 'gslcblas'],),]

setup(name='surf',
      version='0.1',
      description='determine instantaneous liquid interfaces',
      url='http://NOYB.com',
      author='Frank Uhlig',
      author_email='uhlig.frank@gmail.com',
      license='GPLv3',
      packages=['surf'],
      zip_safe=False,
      ext_modules = cythonize(ext_modules),
      )
