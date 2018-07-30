import glob
import os
import shutil
import tarfile
from contextlib import contextmanager
from subprocess import check_call
from tempfile import TemporaryDirectory 

from numpy.distutils.core import Extension, setup


MESASDK_ROOT = '/Applications/mesasdk'


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


def build_mesa(mesa_dir='mesa', clean=True):
    call_mesa_script_sh = os.path.abspath('./call_mesa_script.sh')
    with cd(mesa_dir), TemporaryDirectory() as tmp_path:
        compile_util('ndiff', './utils/ndiff-2.00.tar.gz', tmp_path)
        compile_util('makedepf90', './utils/makedepf90-2.8.8.tar.gz', tmp_path)
        env = {'MESA_DIR': mesa_dir,
               'MESASDK_ROOT': MESASDK_ROOT,
               'PATH': tmp_path + ':' + os.environ['PATH']}
        if clean:
            check_call([call_mesa_script_sh, './clean'], env=env)
        check_call([call_mesa_script_sh, './install'], env=env)


os.environ['PATH'] = MESASDK_ROOT + '/bin/' + ':' + os.environ['PATH']
# build_mesa(clean=True)

extensions = [Extension(
    name='opacity',
    sources=['opacity.f90'],
    include_dirs=['mesa/include'],
    extra_objects=glob.glob('mesa/lib/lib*.a'),
    extra_compile_args=['-fopenmp'],
    extra_link_args=['-lgomp']
)]
setup(ext_modules=extensions)

