import glob
import os
import shutil
import tarfile
from contextlib import contextmanager
from subprocess import check_call
from tempfile import TemporaryDirectory, NamedTemporaryFile

import numpy as np
from Cython.Build import cythonize
from setuptools import setup, Extension, Command


MESASDK_ROOT = os.path.abspath(
    os.environ.get('MESASDK_ROOT', '/Applications/mesasdk')
)
MESASDK_INCLUDE = os.path.join(MESASDK_ROOT, 'include')
MESASDK_LIBS = [
    os.path.join(MESASDK_ROOT, 'lib'),
    os.path.join(MESASDK_ROOT, 'math-slots/default/lib')
]
MESA_DIR = os.path.abspath(
    os.environ.get('MESA_DIR', './mesa/')
)
MESA_INCLUDE = os.path.join(MESA_DIR, 'include')
MESA_LIB = os.path.join(MESA_DIR, 'lib')
os.environ['PATH'] = os.path.join(MESASDK_ROOT, 'bin') + ':' + os.environ['PATH']
os.environ['CC'] = os.path.join(MESASDK_ROOT, 'bin/gcc')


def build_f90(sources, *, shell=None):
    if shell is None:
        base_popenargs = []
    else:
        base_popenargs = [shell]
    objects = []
    for source in sources:
        obj = os.path.splitext(source)[0] + '.o'
        objects.append(obj)
        check_call(base_popenargs
                   + ['gfortran',
                      '-O3',
                      '-I' + MESA_INCLUDE,
                      '-fopenmp',
                      '-fPIC',
                      '-c',
                      source,
                      '-o' + obj])
    return objects


def fortranize(extensions):
    for ext in extensions:
        f90_sources = [f for f in ext.sources if f.endswith('.f90')]
        objects = build_f90(f90_sources)
        for source in f90_sources:
            ext.sources.remove(source)
        ext.extra_objects += objects
    return extensions


def main():
    from sys import argv

    extensions = [Extension(
        name='opacity',
        sources=['src/opacity.f90', 'src/opacity.pyx'],
        include_dirs=[MESASDK_INCLUDE, MESA_INCLUDE, np.get_include()],
        library_dirs=MESASDK_LIBS + [MESA_LIB],
        extra_compile_args=['-fopenmp'],
        extra_link_args=[
            # MESA SDK:
            '-lcrlibm', '-lcrmath', '-lhdf5', '-lhdf5_fortran', '-lgomp',
            '-llapack', '-lblas',
            # MESA:
            '-leos', '-lkap', '-lchem', '-linterp_2d',
            '-linterp_1d', '-lnum', '-lmath', '-lmtx', '-lconst', '-lutils',
            '-lhdf5io', '-lauto_diff',
        ],
    )]

    extensions = fortranize(extensions)
    extensions = cythonize(extensions, annotate=True, force=True)

    setup(
        ext_modules=extensions,
    )


if __name__ == '__main__':
    main()

