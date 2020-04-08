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


@contextmanager
def cd(path):
    old_path = os.getcwd()
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(old_path)


def compile_util(name, tar, dest):
    with TemporaryDirectory() as src:
        with tarfile.open(tar) as tar:
            tar.extractall(src)
        with cd(os.path.join(src, os.listdir(src)[0])):
            check_call(['./configure'])
            check_call(['make', name])
            shutil.move(name, dest)


def get_build_env(path=None):
    env = {'MESA_DIR': MESA_DIR,
            'MESASDK_ROOT': MESASDK_ROOT}
    if path is not None:
        env['PATH'] = path + ':' + os.environ['PATH']
    return env


def build_mesa(clean=True):
    call_mesa_script_sh = os.path.abspath('./call_mesa_script.sh')
    with cd(MESA_DIR), TemporaryDirectory() as tmp_path:
        compile_util('ndiff', './utils/ndiff-2.00.tar.gz', tmp_path)
        compile_util('makedepf90', './utils/makedepf90-2.8.8.tar.gz', tmp_path)
        env = get_build_env(tmp_path)
        if clean:
            check_call([call_mesa_script_sh, './clean'], env=env)
        check_call([call_mesa_script_sh, './install'], env=env)


class BuildMesaCommand(Command):
    """Build Mesa Setup Command"""
    description = "Build mesa"
    user_options = [('clean', None, 'clean before build',),]

    def initialize_options(self):
        self.clean = None

    def finalize_options(self):
        pass

    def run(self):
        build_mesa(clean=self.clean)


def build_f90(sources=('opacity.f90',)):
    env = get_build_env()
    objects = []
    for source in sources:
        obj = os.path.splitext(source)[0] + '.o'
        objects.append(obj)
        with NamedTemporaryFile(suffix='.f90') as tf:
            check_call(['cpp', source, tf.name])
            check_call(['gfortran',
                        '-O3',
                        '-I' + MESA_INCLUDE,
                        '-fopenmp',
                        '-fPIC',
                        '-c',
                        tf.name,
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

    BUILDMESA = 'buildmesa' in argv

    extensions = [Extension(
        name='opacity',
        sources=['opacity.f90', 'opacity.pyx'],
        include_dirs=[MESASDK_INCLUDE, MESA_INCLUDE, np.get_include()],
        library_dirs=MESASDK_LIBS + [MESA_LIB],
        extra_compile_args=['-fopenmp'],
        extra_link_args=[
            # MESA SDK:
            '-lcrlibm', '-lcrmath', '-lhdf5', '-lhdf5_fortran', '-lgomp',
            '-llapack', '-lblas',
            # MESA:
            '-lnet', '-leos', '-lkap', '-lrates', '-lchem', '-linterp_2d',
            '-linterp_1d', '-lnum', '-lmath', '-lmtx', '-lconst', '-lutils',
        ],
    )]

    if not BUILDMESA:
        extensions = fortranize(extensions)
        extensions = cythonize(extensions, annotate=True, force=True)

    setup(
        ext_modules=extensions,
        cmdclass={'buildmesa': BuildMesaCommand,},
    )


if __name__ == '__main__':
    main()

