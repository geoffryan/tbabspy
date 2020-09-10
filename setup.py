import setuptools
import subprocess
from pathlib import Path
from setuptools import setup, Extension
# from numpy.distutils.core import setup as npsetup
# from numpy.distutils.core import Extension as npExtension
from numpy.distutils.fcompiler import new_fcompiler
from setuptools.command.build_ext import build_ext
from distutils.errors import DistutilsSetupError
from distutils import log as distutils_logger

import numpy as np

version = {}
with open("tbabspy/version.py", "r") as f:
    exec(f.read(), version)

with open("README.md", "r") as f:
    long_description = f.read()

phfit2Src = "tbabspy/phfit2Src/phfit2.f"
phfit2Obj = "tbabspy/phfit2Src/phfit2.o"

coreInc = [np.get_include()]

coreSrc = ["tbabspy/tbabscoremodule.c", "tbabspy/util.c",
           "tbabspy/tbabsSrc/tbvabs.c"]
coreDep = ["tbabspy/util.h",
           phfit2Src,
           "tbabspy/tbabsSrc/tbvabs.h",
           "tbabspy/tbabsSrc/fephoto.h",
           "tbabspy/tbabsSrc/nephoto.h",
           "tbabspy/tbabsSrc/ophoto.h"]
coreObj = [phfit2Obj]



# phfit2Src = ["tbabspy/phfit2Src/phfit2.f"]

# phfit2Module = Extension('tbabspy.phfit2', sources=phfit2Src,
                         # extra_compile_args=['-Wall'])

coreModule = Extension('tbabspy.tbabscore',
                        sources=coreSrc,
                        depends=coreDep,
                        # libraries=coreLib, library_dirs=['tbabspy'],
                        include_dirs=coreInc,
                        extra_objects=coreObj,
                        extra_compile_args=['-Wall'])

# grabbed from lucianopaz's answer to this question:
# QUESTION:   https://stackoverflow.com/questions/41169711/
#  python-setuptools-distutils-custom-build-for-the-extra-package-with-makefile
# 
# ANSWER:  https://stackoverflow.com/a/48641638/2418959
#
# This may not be the whole story, if this breaks on other platforms something
# like https://stackoverflow.com/a/49765622/2418959 may be necessary to ensure
# all necessary libraries are packaged in the wheel.
class custom_build_ext(build_ext, object):

    my_ext = coreModule.name

    def build_extension(self, ext):

        # If you're not building tbabscore, do the usual thing
        if ext.name != self.my_ext:
            super(custom_build_ext, self).build_extension(ext)

        else:
            # First we need to build the phfit2 Fortran library

            # Hack to compile phfit2.f manually
            # distutils_logger.info("Compiling " + phfit2Src
            #                       + " into " + phfit2Lib)
            # subprocess.call(["gfortran", "-O3", "-Wall",
            #                  "-o", phfit2Obj,
            #                  "-c", phfit2Src])
            # distutils_logger.info("Compiled " + phfit2Lib)

            distutils_logger.info("Compiling " + phfit2Src + " with numpy")

            # Use numpy's fortan enabled distutils
            # These lines are copied from numpy.distutil.command.build_ext's
            # run() method.
            f77_compiler = new_fcompiler(compiler=None,
                                         verbose=self.verbose,
                                         dry_run=self.dry_run,
                                         force=self.force,
                                         requiref90=False,
                                         c_compiler=self.compiler)

            if f77_compiler:
                ctype = f77_compiler.compiler_type
                f77_compiler.customize(self.distribution)
            if f77_compiler and f77_compiler.get_version():
               f77_compiler.customize_cmd(self)
               f77_compiler.show_customization()
            else:
                self.warn('f77_compiler=%s is not available.' % (ctype))
                f77_compiler = None

            if f77_compiler is None:
                raise DistutilsSetupError("no Fortran compiler found")

            f77_compiler.extra_f77_compile_args = ["-Wall"]

            # We should now have a setup compiler compliant with the
            # python installation.  Compile the .f source!
            f_objects = f77_compiler.compile([phfit2Src],
                                             output_dir=".",
                                             macros=[],
                                             include_dirs=[],
                                             debug=False,
                                             extra_postargs=[])

            f_obj = f_objects[0]

            # I *think* with output_dir="." the object file will be put
            # alongside the source file (ie. phfit2Obj).  But I don't know if
            # this is true on all systems, so check and adjust the
            # extra_objects parameter as necessary.
            default_path = Path(phfit2Obj)
            if not default_path.exists()\
                    or not default_path.samefile(Path(f_obj)):
                distutils_logger.info("phfit2.f compiled into unexpected"
                                      + " location")
                ext.extra_objects += f_objects
            
            distutils_logger.info("Compiled "
                                  + " ".join([str(x) for x in f_objects]))
            
            # Now that the fortran file has been compiled we can proceed
            # with a normal build!
            super(custom_build_ext, self).build_extension(ext)



setup(
    name='tbabspy',
    version=version['__version__'],
    author="Geoffrey Ryan",
    author_email="gsryan@umd.edu",
    description='tbabs X-ray absorption models',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://github.com/geoffryan/tbabspy',
    packages=['tbabspy'],
    ext_modules=[coreModule],
    cmdclass={'build_ext': custom_build_ext},
    classifiers=[
        "Programming Language :: Python :: 3",
        "Programming Language :: C",
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Astronomy"],
    install_requires=['numpy>=1.10'],
    extras_require={
        'docs': ['numpydoc']
        },
    project_urls={
        "Source Code": "https://github.com/geoffryan/tbabspy",
        "tbabs XSPEC manual": ("https://heasarc.gsfc.nasa.gov/xanadu/"
                               +"xspec/manual/node265.html"),
        "tbabs beta development": ("https://pulsar.sternwarte.uni-erlangen.de"
                                   +"/wilms/research/tbabs/"),
        "phfit2 Code": "http://www.pa.uky.edu/~verner/photo.html"}
    )
