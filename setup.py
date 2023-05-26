from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
import numpy

#Adding multiple files..
#ext		=		Extension('_header_read_cy', '_header_Fring_cy', '_header_geometric_cy', '_header_gps_cy', '_header_delay_cy', ['call_to_read_dev_dev_dev.pyx, header_Fringe_Search_dev_dev.pyx', 'geometric.pyx', 'header_gps_dev.pyx', 'delay.pyx'],language ="c++",)

setup(
  name='Cython modules',
  ext_modules=[ Extension('CYTHON._header_Fring_cy',        ['CYTHON/header_Fringe_Search_proc.pyx'],language ="c++",), \
                Extension('CYTHON._header_read_cy',         ['CYTHON/call_to_read_dev_dev_dev.pyx'],language ="c++",), \
                Extension('CYTHON._header_geometric_cy',    ['CYTHON/geometric.pyx'],language ="c++",), \
                Extension('CYTHON._header_gps_cy',          ['CYTHON/header_gps_dev.pyx'],language ="c++",), \
                Extension('CYTHON._header_delay_cy',        ['CYTHON/delay.pyx'],language ="c++",)], \
                include_dirs=[numpy.get_include()],\
                cmdclass={'build_ext': build_ext, },
)
