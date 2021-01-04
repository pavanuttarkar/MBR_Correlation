# DS
from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

ext_modules = [
    Extension("_header_Fring_cy", ["header_Fringe_Search.pyx"], ), 
    Extension('_header_read_cy', ['call_to_read_dev.pyx'], language ="c++"),
    Extension('_header_geometric_cy', ['geometric.pyx'],)
    ]

for e in ext_modules:
    e.cython_directives = {"embedsignature": True}

setup(
    name='All header for MBR Data Format',
    cmdclass={'build_ext': build_ext},
    ext_modules=ext_modules
)


# libraries = ["/home/devansh/Desktop/Schwartz/MBR_Corr/test/spdocs_test/source/internal_loop2"]
# libraries = ["/home/devansh/Desktop/Schwartz/MBR_Corr/test/spdocs_test/source/internal_loop2"]