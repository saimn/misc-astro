import numpy as np
from setuptools import setup, Extension
from Cython.Build import cythonize

setup(name="hackinghour",
      ext_modules=cythonize([Extension('fastloop', ['fastloop.pyx'],
                                      include_dirs=[np.get_include()])])
      )
