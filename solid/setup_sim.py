#!/usr/bin/python
import pyximport; pyximport.install()
from distutils.core import setup
from Cython.Build import cythonize
from Cython.Distutils import build_ext
import sys
import numpy as np

setup(name = 'SetupTb',
      ext_modules=cythonize(['poly_chain.pyx','poly_spike.pyx','simulate_prewetting.pyx',]),
      include_dirs = [np.get_include()],
      cmdclass={'build_ext': build_ext}
)

#simulate_kinetic_events.pyx'
