from os import path, popen, uname
from setuptools import setup
from setuptools.command.install import install
from distutils.command.build import build
from subprocess import call

USE_OPENMP=False
BASEPATH = path.abspath(path.join(path.dirname(__file__),'..'))
LONG_DESC= "TODO: DESCRIPTION OF LIBCMPB-PYTHON HERE"



class CMPB_Build(build):
  def run(self):
    EXT = ".dylib" if uname()[0] == "Darwin" else ".so"

    # run original build code
    build.run(self)

    # build libcmpb
    cmd = [ 'make', 'all' ]
    message = ""

    # run Make for each condition (make CPU/GPU, 32/64)
    def compile():
          call(cmd, cwd=BASEPATH)

    self.execute(compile, [], message)

    # copy Make output to library build folder
    self.mkpath(self.build_lib)
    for target in target_files:
          self.copy_file(path.join(BASEPATH,"libcmpb" + EXT), self.build_lib)


class CMPB_Install(install):
    def initialize_options(self):
        install.initialize_options(self)
        self.build_scripts = None

    def finalize_options(self):
        install.finalize_options(self)
        self.set_undefined_options('build', ('build_scripts', 'build_scripts'))

    def run(self):
        # run original install code
        install.run(self)

        # install CMPB executable
        self.copy_tree(self.build_lib, self.install_lib)

setup(
    name='mathprogbase-python',
    version='0.0.0',
    author='Miles Lubin, Steven Diamond, Baris Ungun',
    author_email='mlubin@mit.edu',
    url='https://github.com/mlubin/cmpb',
    package_dir={'optkit': 'optkit'},
    packages=['cmpb'],
    license='TODO: LICENSE',
    zip_safe=False,
    description='TODO: SHORT DESCRIPTION',
    long_description=LONG_DESC,
    install_requires=["numpy >= 1.8",
                      "scipy >= 0.13"],
    cmdclass={'build' : CMPB_Build, 'install' : CMPB_Install}
)
