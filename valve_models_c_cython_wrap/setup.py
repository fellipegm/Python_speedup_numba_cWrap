from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
import numpy


sources = ['valve_wrap.pyx', 'valve_models.c']

setup(
    name='First Order Filter Pure C',
    ext_modules=[
            Extension('_valve_models_c_wrap_cython', 
                      sources, 
                      include_dirs=[numpy.get_include()],
                      extra_compile_args=["-O2", '-march=native', '-mtune=native']),
                ],
    cmdclass={'build_ext': build_ext},
)