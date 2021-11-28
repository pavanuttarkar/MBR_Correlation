from numpy.distutils.core import Extension as FExtension
from numpy.distutils.core import setup as fsetup

from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

f_modules = [
    FExtension(
        name="api.internal_loop_RFI",
        sources=["api/internal_standalone_WORKING_with_RFI.f90"],
        libraries=["fftw3", "m"],
        library_dirs=["/usr/lib/x86_64-linux-gnu/"],
        extra_compile_args=[
            "-O3",
            "-L /usr/lib/x86_64-linux-gnu/",
            "-lfftw3",
            "-lm",
            "-DF2PY_REPORT_ON_ARRAY_COPY=1",
        ],
    ),
    FExtension(
        name="api.internal_loop2",
        sources=["api/internal_standalone_WORKING_with_RFI.f90"],
        libraries=["fftw3", "m"],
        library_dirs=["/usr/lib/x86_64-linux-gnu/"],
        extra_compile_args=[
            "-O3",
            "-L /usr/lib/x86_64-linux-gnu/",
            "-lfftw3",
            "-lm",
            "-DF2PY_REPORT_ON_ARRAY_COPY=1",
        ],
    ),
    FExtension(
        name="api.single_file_combined_RFInoRFI",
        sources=["api/internal_standalone_WORKING_with_RFI_for_single_file_dev.f90"],
        libraries=["fftw3", "m"],
        library_dirs=["/usr/lib/x86_64-linux-gnu/"],
        extra_compile_args=[
            "-O3",
            "-L /usr/lib/x86_64-linux-gnu/",
            "-lfftw3",
            "-lm",
            "-DF2PY_REPORT_ON_ARRAY_COPY=1",
        ],
    ),
]

fsetup(
    name="python-fortran",
    ext_modules=f_modules,
)

ext_modules = [
    Extension(
        "api._header_Fring_cy",
        ["api/header_Fringe_Search_proc.pyx"],
        language="c++",
    ),
    Extension(
        "api._header_read_cy",
        ["api/call_to_read_dev_dev_dev.pyx"],
        language="c++",
    ),
    Extension(
        "api._header_geometric_cy",
        ["api/geometric.pyx"],
        language="c++",
    ),
    Extension(
        "api._header_gps_cy",
        ["api/header_gps_dev.pyx"],
        language="c++",
    ),
    Extension(
        "api._header_delay_cy",
        ["api/delay.pyx"],
        language="c++",
    ),
]

for e in ext_modules:
    e.cython_directives = {"language_level": 2}

setup(
    name="Cython modules",
    ext_modules=ext_modules,
    cmdclass={"build_ext": build_ext, "embedsignature": True},
)
